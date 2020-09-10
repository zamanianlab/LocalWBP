library(tidyverse)
library(here)
library(conflicted)
library(biomaRt)

conflict_prefer("filter", "dplyr")

# An example of pulling ortholog data from WormBase ParaSite --------------

# establish the BioMart that will be used
mart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)

# example for a single gene

# pull the desired data
# Note: It can be hard to know the parameters that you will want to use
#       here. The best way to find your parameters is to use the
#       WormBase ParaSite BioMart and click the XML button. The below
#       example will pull orthologs of Bm1711 (Bma-osm-9).
osm9_orthlogs <- getBM(mart = mart, 
                       filters = c("wbps_transcript_id"),
                       value = list("Bm1711a.1"),
                       attributes = c("wbps_transcript_id", "brpahaprjeb497_gene", "caelegprjna13758_gene", "assuumprjna80881_gene", "aslumbprjeb4950_gene", "strattprjeb125_gene")) %>%
  janitor::clean_names() %>%
  pivot_longer(cols = 1:6, names_to = 'name', values_to = 'transcript_id')


# example for a list of genes
seeds <- read_csv(here('aux', 'celegans_aminer.txt'), col_names = 'seed')

# set parameters
filters = c('wbps_gene_id')
value = list(seeds$seed)
attributes = c('wbps_peptide_id',
               'acviteprjeb1697_homolog_ensembl_peptide',
               'anceylprjna72583_homolog_ensembl_peptide',
               'assuumprjna62057_homolog_ensembl_peptide',
               'brmalaprjna10729_homolog_ensembl_peptide',
               'brpahaprjeb497_homolog_ensembl_peptide',
               'brtimoprjeb4663_homolog_ensembl_peptide',
               'buxyloprjea64437_homolog_ensembl_peptide',
               'caelegprjna13758_homolog_ensembl_peptide',
               'diimmiprjeb1797_homolog_ensembl_peptide',
               'drmediprjeb500_homolog_ensembl_peptide',
               'elelapprjeb502_homolog_ensembl_peptide',
               'hacontprjeb506_homolog_ensembl_peptide',
               'hepolyprjeb15396_homolog_ensembl_peptide',
               'lisigmprjeb3075_homolog_ensembl_peptide',
               'loloaprjna246086_homolog_ensembl_peptide',
               'neamerprjna72135_homolog_ensembl_peptide',
               'onflexprjna230512_homolog_ensembl_peptide',
               'onocheprjeb1204_homolog_ensembl_peptide',
               'onvolvprjeb513_homolog_ensembl_peptide',
               'parediprjna186477_homolog_ensembl_peptide',
               'paunivprjna386823_homolog_ensembl_peptide',
               'patricprjeb515_homolog_ensembl_peptide',
               'prpaciprjna12644_homolog_ensembl_peptide',
               'stcarpprjna202318_homolog_ensembl_peptide',
               'strattprjeb125_homolog_ensembl_peptide',
               'tocaniprjna248777_homolog_ensembl_peptide',
               'trspirprjna12603_homolog_ensembl_peptide',
               'trmuriprjeb126_homolog_ensembl_peptide',
               'wubancprjna275548_homolog_ensembl_peptide'
)

# pull the desired data
aminer_orthologs <- getBM(mart = mart, 
                   filters = filters,
                   value = value,
                   attributes = attributes) %>%
  janitor::clean_names() %>%
  pivot_longer(cols = 1:30) %>%
  dplyr::select(name, transcript_id = value) %>%
  filter(transcript_id != '') %>%
  distinct() %>%
  arrange(transcript_id)
