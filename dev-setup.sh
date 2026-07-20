#!/usr/bin/env bash

# this is meant to be done after a full conda install (or setting up all dependencies)
## making sure we are in the gtotree-dev conda environment
if [ "${CONDA_DEFAULT_ENV}" != "gtotree-dev" ]; then
    printf "\n    This should be run in the 'gtotree-dev' conda environment..\n"
    printf "    You know this, Mike...\n\n"
    return 1 2>/dev/null || exit 1
fi

rm -rf build/ gtotree.egg-info/

pip uninstall -y gtotree 2>/dev/null || true

pip install --force-reinstall --no-build-isolation -e .

## if changing conda versions and wanting to install locally entirely (rather than using a prior official conda install of gtotree)
# conda build -c conda-forge -c bioconda conda-recipe/ ; conda create -n gtotree-dev -c conda-forge -c bioconda --use-local gtotree -y
