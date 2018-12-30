#!/bin/sh

## adding conda channels
conda config --add channels defaults 2> /dev/null
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels au-eoed

## creating GToTree environment and installing dependencies
conda create -n gtotree biopython hmmer muscle trimal fasttree prodigal taxonkit gnu-parallel --yes

## activating environment
source activate gtotree

## creating directory for conda-env-specific source files
mkdir -p ${CONDA_PREFIX}/etc/conda/activate.d

## adding GToTree bin path and GToTree_HMM_dir variable:
echo '#!/bin/sh'" \


export PATH=\"$(pwd)/bin:"'$PATH'\"" \

export GToTree_HMM_dir=\"$(pwd)/hmm_sets\"" >> ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh

## downloading ncbi tax database for taxonkit and setting variable
mkdir ncbi_tax_info
cd ncbi_tax_info

curl --silent --retry 10 -O ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzvf taxdump.tar.gz 2> /dev/null

echo "export TAXONKIT_DB=$(pwd)" >> ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh

cd ../

# re-activating environment so variable and PATH changes take effect
source activate gtotree
