# data wrangling/plotting
library(readr)
library(dplyr)
library(stringr)

# misc
library(conflicted)

conflict_prefer("filter", "dplyr")

# get arguments
args <- commandArgs(trailingOnly = TRUE)

# get the base directory
pattern <- str_glue(".*", args[2], ".*/")
save_location <- str_extract(args[1], pattern)

# read GTF
gtf <- read_delim(args[1],
                  delim = '\t',
                  comment = '#',
                  col_names = c('scaffold', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'))

parsed_gtf <- select(gtf, -frame, -score, -source) %>%
  mutate(gene_id = str_remove(attribute, ';.*'),
         gene_id = str_remove_all(gene_id, '\''),
         gene_id = str_remove(gene_id, 'gene_id '),
         gene_id = str_remove_all(gene_id, '\"'),
         transcript_id = str_remove(attribute, '.*transcript_id '),
         transcript_id = str_remove(transcript_id, ';.*'),
         transcript_id = str_remove_all(transcript_id, '\''),
         transcript_id = str_remove(transcript_id, '\\.[0-9]'),
         transcript_id = str_remove_all(transcript_id, '\"'),
         biotype = str_extract(attribute, 'transcript_biotype \"[a-z|_]*\"'),
         biotype = str_remove(biotype, 'transcript_biotype '),
         biotype = str_remove_all(biotype, '\"'))

transcripts <- filter(parsed_gtf, feature == 'transcript') %>%
  select(chr = scaffold, gene_id, transcript_id, start, end, strand, biotype, attribute)

# Write transcripts -------------------------------------------------------
# create the file to write
transcript_outfile <- str_glue(save_location, args[2], ".gtf.rds")

# save the parsed GTF as an RDS
saveRDS(transcripts, transcript_outfile)

# Write exons -------------------------------------------------------------
exons <- filter(parsed_gtf, feature == 'exon') %>%
  select(chr = scaffold, gene_id, transcript_id, start, end, strand, biotype, attribute)

# create the file to write
exon_outfile <- str_glue(save_location, args[2], ".exons.rds")

# save the parsed GTF as an RDS
saveRDS(exons, transcript_outfile)

# save the parsed GTF as an RDS
saveRDS(gtf, exon_outfile)

