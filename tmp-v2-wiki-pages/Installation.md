[GToTree](https://github.com/AstrobioMike/GToTree/wiki) runs in a Unix-like command-line environment. This means it will work on Mac and Linux computers in the standard terminal programs available with them. And to use GToTree on a Windows computer, I would recommend installing the [Windows Subsystem for Linux (WSL)](https://docs.microsoft.com/en-us/windows/wsl/install), then when in the WSL terminal, install a Linux version of [miniconda](https://astrobiomike.github.io/unix/conda-intro#getting-and-installing-conda). Then installing with `conda` as shown below will work in the WSL environment 👍 

## Conda quickstart!
If you don't already have the glorious package manager [conda](https://conda.io/docs/), I **highly** recommend you get it. This really isn't the venue to go into why it's so helpful, but it really is, I promise 🙂 

To get conda up and running (which is very quick), you can follow the instructions to install miniconda (a light-weight version) for your appropriate system starting from [here](https://conda.io/en/latest/miniconda.html). You will want a python 3.X version, and more than likely a 64-bit version. And if you'd like to learn more about conda sometime, I have [an introduction page here](https://astrobiomike.github.io/unix/conda-intro) 🙂 

---

The following line will create a gtotree conda environment and install GToTree, you want to run these in the base conda environment:

```bash
conda create -y -n gtotree -c astrobiomike -c conda-forge -c bioconda gtotree

# IF you are on mac, and aren't sure if you installed an x86_64 miniconda, then do this
CONDA_SUBDIR=osx-64 conda create -y -n gtotree -c astrobiomike -c conda-forge -c bioconda gtotree
```

**DONE!**

Now you should be able to enter and exit the environment with `conda activate gtotree` and `conda deactivate`. If you enter the environment and run the following:

```bash
gtt hmms
```

It will print out where the GToTree default HMMs directory is located, and list the available pre-built SCG-sets. Running `gtotree -h` gives the condensed help menu, and `gtotree -s` the detailed one.

Running `gtt` by itself gives an overview of all the helper subcommands available.

---

## Reference data

GToTree uses a few reference databases (GTDB taxonomy, NCBI assembly summary, and Pfam and KOFamScan profiles if you use them). These are **downloaded automatically the first time they're needed**, so for most people there is nothing to do here.

But if you'd like to find them, or fetch or update them at any time, you can access them with the following `gtt data` helpers.

To see where GToTree is storing them:

```bash
gtt data locations check
```

To set where they should be stored:

```bash
gtt data locations set
```

To download any (again, will be done automatically when needed anyway, but if you need to ahead of time, here's how you can):

```bash
gtt data get gtdb-data
gtt data get ncbi-assembly-data
gtt data get pfam-data
gtt data get kofamscan-data
```

And adding `-f`/`--force-update` will re-download something already present, which can be good if you want to make sure you have the latest available.

---

## Test run
You can run a quick end-to-end test of the installed environment like so:

```bash
gtt test
```

If that works, it will print a message saying the test completed successfully.

The test cleans up after itself. Add `-k`/`--keep` if you want to poke around at the output directory it makes (`test-gtotree-output/`) afterward.

---

## Updating to a newer version
If wanting to update to the latest GToTree version, it is best to remove the previous conda environment and install fresh. This can be done as follows:

```bash
# from outside the gtotree conda environment (assuming that's what it was named like the install above)

# then re-installing in a new environment same as above (this will remove and overwrite the prior one)
conda create -y -n gtotree -c astrobiomike -c conda-forge -c bioconda gtotree

# or to be sure to grab the right stuff if on a mac, you can add this variable up front
CONDA_SUBDIR=osx-64 conda create -y -n gtotree -c astrobiomike -c conda-forge -c bioconda gtotree
```

Then the new environment can be activated with `conda activate gtotree`.


---

## Versions of utilized programs
Every GTOtree run produces a citations.txt file that holds the versions and citation info of all primary programs used. That should be the main source of truth for version tracking and citation info for GToTree and the wonderful tools it depends on. 

You can also check specific versions in your active conda gtotree environment with `conda list`.
