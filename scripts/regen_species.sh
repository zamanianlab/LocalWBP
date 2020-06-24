#!bin/bash

### Use this script to regenerate the entire list of species hosted at WormBase ParaSite

### Preparation
proj="LocalWBP"

gh_dir="${GIT_PATH}/${proj}"

curl -l ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/current/species/ > "${gh_dir}"/aux/species_all.txt
