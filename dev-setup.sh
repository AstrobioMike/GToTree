#!/usr/bin/env bash

if [ "${CONDA_DEFAULT_ENV}" != "gtotree-dev" ]; then
    printf "\n    This should be run in the 'gtotree-dev' conda environment..\n"
    printf "    You know this, Mike...\n\n"
    return 1 2>/dev/null || exit 1
fi

if [[ "${1:-}" == "conda" && "${BASH_SOURCE[0]}" == "$0" ]]; then
    printf "\n    The 'conda' option must be used with sourcing this script:\n"
    printf "        source dev-setup.sh conda\n\n"
    exit 1
fi


# new conda build if wanted
if (($# > 0)); then
    if [[ "$1" == "conda" ]]; then
        eval "$(conda shell.bash hook)"
        conda deactivate
        conda build -c conda-forge -c bioconda conda-recipe/
        conda create -n gtotree-dev -c conda-forge -c bioconda --use-local gtotree -y
        conda activate gtotree-dev
    elif [[ "$1" != "conda" ]]; then
        printf "\n    Unrecognized argument: $1\n"
        printf "    This script only accepts 'conda' as an argument.\n\n"
        return 1 2>/dev/null || exit 1
    fi
fi

# normal pip update
rm -rf build/ gtotree.egg-info/
pip uninstall -y gtotree 2>/dev/null || true
pip install --force-reinstall --no-build-isolation -e .
