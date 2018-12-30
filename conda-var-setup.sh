#!/bin/sh

## creating directory for conda-env-specific source files
mkdir -p ${CONDA_PREFIX}/etc/conda/activate.d

## adding GToTree bin path and GToTree_HMM_dir variable:
echo '#!/bin/sh'" \


export PATH=\"$(pwd)/bin:"'$PATH'\"" \

export GToTree_HMM_dir=\"$(pwd)/hmm_sets\"" >> ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh

# re-activating environment
source activate gtotree
