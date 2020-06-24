#!bin/bash

### Preparation
proj="LocalWBP"

gh_dir="${GIT_PATH}/${proj}"

# for local MacOS
local_dir="$HOME/Box/ZamanianLab/Data/WBP"

# for server
# local_dir="~/WBP"

### Define wormbase source links and download genomes

wbp_prefix="ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS14/species"

# Species list based on genome quality
species="${gh_dir}"/aux/species_selected.txt

# Download data ----------------------------------------------------------------

# -N: only download newer versions
# -nc: no clobber; ignore server files that aren't newer than the local version
# -r: recursive
# -nH: don't mimick the server's directory structure
# -cut-dirs=7: ignore everything from pub to species in the recursive search
# --no-parent: don't ascend to the parent directory during a recursive search
# -A: comma-separated list of names to accept
# -P:
while IFS= read -r line
do
  species_dl="$wbp_prefix/$line/"
  printf ${species_dl}"\n"
  wget -nc -r -nH --cut-dirs=7 --no-parent --reject="index.html*" -A 'canonical_geneset.gtf.gz','genomic_masked.fa.gz','protein.fa.gz','CDS_transcripts.fa.gz' $species_dl -P $local_dir
done <"$species"

# Make BLAST databases ---------------------------------------------------------

while IFS= read -r line
do
  # split the line into an array with species + BioProject
  array=($(echo "$line" | sed 's/\// /g'))
  species_folder=$local_dir/${array[0]}/${array[1]}
  gunzip -k $species_folder/*.genomic*.gz
  makeblastdb -in $species_folder/*.genomic*.fa -dbtype nucl
  rm $species_folder/*.genomic*.fa

  gunzip -k $species_folder/*.protein*.gz
  makeblastdb -in $species_folder/*.protein.fa -dbtype prot
  rm $species_folder/*.protein.fa
done <"$species"

# make BLAST databases of concatenated FASTA files
mkdir $local_dir/all
find $local_dir -name '*.genomic*.gz' -exec cat {} + > $local_dir/all/all.genomic_masked.fa.gz
gunzip -k $local_dir/all/all.genomic_masked.fa.gz
makeblastdb -in $local_dir/all/all.genomic_masked.fa -dbtype nucl

find $local_dir -name '*.protein.fa.gz' -exec cat {} + > $local_dir/all/all.protein.fa.gz
gunzip -k $local_dir/all/all.protein.fa.gz
makeblastdb -in $local_dir/all/all.protein.fa -dbtype prot

# Create GTF and exon RDS files ------------------------------------------------

while IFS= read -r line
do
  # split the line into an array with species + BioProject
  array=($(echo "$line" | sed 's/\// /g'))
  species_folder=$local_dir/${array[0]}/${array[1]}
  gunzip -k $species_folder/*.gtf.gz
  Rscript scripts/parse_GTF.R $species_folder/*.gtf ${array[0]}
  rm $species_folder/*.gtf
done <"$species"
