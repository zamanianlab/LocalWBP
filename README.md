# LocalWBP
A collection of scripts for organizing a local version of WormBase ParaSite

## Scripts
**regen_species.sh**
  1. Regenerates `aux/species_all.txt`, which is a list of all available species on the WormBase ParaSite FTP

**WBP_download.sh**
  1. Downloads the GTF, genomic_masked.fasta, protein.fasta, and CDS_transcripts.fasta files for each species/BioProject in `aux/species_selected.txt`  
      - The user must create this file which includes the preferred BioProject (many species have more than one assembly)  
      - This module will only overwrite current local files if the server files have been updated, otherwise those files will be skipped
  2. Makes BLAST databases for the protein and genomic_masked sequences from each species; deletes gunzipped FASTA files after creating the database, as they are unnecessary for BLAST searches
  3. Makes BLAST databases for the concatenated protein and genomic_masked sequences; saves these databases in `all/`
  4. Runs `scripts/parse_GTF.R` on the GTF for each species (see below for description)

**parse_GTF.R**
  1. Reads GTF and removes frame, score, and source columns
  2. Creates new columns for gene_id, transcript_id, and biotype
  3. Creates two RDS files for faster loading into R:  
      - `genus_species.gtf.rds`: the parsed GTF with only the transcript features  
      - `genus_species.exon.rds`: the parsed GTF with only the exon features (useful for plotting gene models with `ggplot2`)

**WBP_biomart.R**
  - A short template for querying the WormBase ParaSite BioMart within R, which a focus on pulling orthologs of a given list of genes

**GO_template.R**
  - A template for performing GO enrichment tests using GO terms pulled from WormBase ParaSite or VectorBase
  - File requirements:  
    - GO mappings: a file with a strict organizational definition - the first column should contain gene_ids and the second column should contain a comma-separated list of GO terms associated with a given gene_id; examples can be found in `aux/GO_mappings`) and can be regenerated with this script  
    - a list of genes of interest; an example can be found at `aux/gene.list.sus.csv`
