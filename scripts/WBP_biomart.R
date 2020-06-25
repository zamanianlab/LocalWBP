library(tidyverse)
library(here)
library(conflicted)
library(biomaRt)

conflict_prefer("filter", "dplyr")

# An example of pulling ortholog data from WormBase ParaSite --------------

# establish the BioMart that will be used
mart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)

# pull the desired data
# Note: It can be hard to know the parameters that you will want to use
#       here. The best way to find your parameters is to use the
#       WormBase ParaSite BioMart and click the XML button. The below
#       example will pull orthologs of Bm1711 (Bma-osm-9).
osm9_orthlogs <- getBM(mart = mart, 
                       filters = c("species_id_1010", "wbps_transcript_id"),
                       value = list("brmalaprjna10729", "Bm1711a.1"),
                       attributes = c("production_name_1010", "wbps_gene_id", "brpahaprjeb497_gene", "caelegprjna13758_gene", "assuumprjna80881_gene", "aslumbprjeb4950_gene", "strattprjeb125_gene")) %>%
  janitor::clean_names()

