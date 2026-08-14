[GToTree](https://github.com/AstrobioMike/GToTree/wiki) runs in a Unix-like command-line environment. This means it will work on Mac and Linux computers in the standard terminal programs available with them. And to use GToTree on a Windows computer, I would recommend installing the [Windows Subsystem for Linux (WSL)](https://docs.microsoft.com/en-us/windows/wsl/install), then when in the WSL terminal, install a Linux version of [miniconda](https://astrobiomike.github.io/unix/conda-intro#getting-and-installing-conda). Then installing with `conda` as shown below will work in the WSL environment 👍 

## Conda quickstart!
If you don't already have the glorious package manager [conda](https://conda.io/docs/), I **highly** recommend you get it. This really isn't the venue to go into why it's so helpful, but it really is, I promise 🙂 

To get conda up and running (which is very quick), you can follow the instructions to install miniconda (a light-weight version) for your appropriate system starting from [here](https://conda.io/en/latest/miniconda.html). You will want a python 3.X version, and more than likely a 64-bit version. And if you'd like to learn more about conda sometime, I have [an introduction page here](https://astrobiomike.github.io/unix/conda-intro) 🙂 

---

The following line will create a gtotree conda environment and install GToTree, you want to run these in the base conda environment:

```bash
conda create -y -n gtotree -c astrobiomike -c conda-forge -c bioconda gtotree
```

**DONE!**

Now you should be able to enter and exit the environment with `conda activate gtotree` and `conda deactivate`. If you enter the environment and run the following:

```bash
gtt hmms
```

It will print out where the GToTree default HMMs directory is located, and list the available pre-built SCG-sets. Running `gtotree -h` gives the condensed help menu, and `gtotree -s` the detailed one.

Running `gtt` by itself gives an overview of all the packaged helper subcommands.

---

## Reference data

GToTree uses a few reference databases (GTDB taxonomy, NCBI assembly summary, and – if you use those features – Pfam and KOFamScan profiles). These are **downloaded automatically the first time they're needed**, so for most people there is nothing to do here.

If you'd like to fetch or update them ahead of time, e.g. on a compute node without internet access at runtime:

```bash
gtt data get gtdb-data
gtt data get ncbi-assembly-data
gtt data get pfam-data
gtt data get kofamscan-data
```

Add `-f`/`--force-update` to re-download something already present.

To see where GToTree is storing/looking for these:

```bash
gtt data locations check
```

And to change those locations (useful on shared systems, or where a home directory has a quota):

```bash
gtt data locations set
```

---

## Test run
You can run a quick end-to-end test of the installed environment like so:

```bash
gtt test
```

<!-- TODO(Mike): paste the current `gtt test` output here -->

The test cleans up after itself. Add `-k`/`--keep` if you want to poke around at the output directory it makes (`test-gtotree-output/`) afterward.

---

## Updating to a newer version
If wanting to update to the latest GToTree version, it is best to remove the previous conda environment and install fresh. This can be done as follows:

```bash
# from outside the gtotree conda environment (assuming that's what it was named like the install above)
conda env remove -n gtotree

# then re-installing in a new environment same as above
conda create -y -n gtotree -c astrobiomike -c conda-forge -c bioconda gtotree
```

Then the new environment can be activated with `conda activate gtotree`.

> **Updating from v1 to v2?** v2 is a full rewrite and several commands and outputs changed names – the helper programs are now subcommands of `gtt` (e.g. `gtt-hmms` is now `gtt hmms`), and the packaged SCG-sets are rebuilt against GTDB. See the [CHANGELOG](https://github.com/AstrobioMike/GToTree/blob/master/CHANGELOG.md) before re-running old commands. Installing into a fresh environment as above is the way to go.

---

## Installation without conda (not recommended)

Again, the conda installation is **highly** recommended as it is more robust across different systems. But to try installing without conda, GToTree v2 is a python package, so it can be installed with `pip` from a release tarball (be sure to change the version below to the latest found [here](https://github.com/AstrobioMike/GToTree/releases/latest)):

```bash
curl -L https://github.com/AstrobioMike/GToTree/archive/v2.0.0.tar.gz -o GToTree-v2.0.0.tar.gz
tar -xzvf GToTree-v2.0.0.tar.gz
cd GToTree-2.0.0
pip install .
```

That installs the `gtotree`/`GToTree` and `gtt` commands and their python dependencies. You will still need to install the non-python dependencies yourself – see [below](installation#installing-dependencies-without-conda).

### Add path to included HMM files
The packaged SCG-set HMMs live in the `hmm_sets` directory of the source. So that you don't need to provide a full path to them, set an environment variable pointing at it:

```bash
cd hmm_sets/ # from the unpacked source directory
echo "export GToTree_HMM_dir=\"$(pwd)/\"" >> ~/.bash_profile
source ~/.bash_profile
```

You can run `gtt hmms` with no arguments to make sure the default HMM directory is set, and see what taxa the currently available SCG-sets target.

---

## Installing dependencies without conda
By **far**, the easiest way to get all the dependencies up and running is with [conda](https://conda.io/docs/) as done above. But if you don't want to use conda, here's what's needed.

The python dependencies (biopython, dendropy, pandas, pyarrow, pyhmmer, requests, rich, rich-argparse, tqdm, argcomplete) are handled by `pip install .` above. The rest are external programs that need to be on your `PATH`.

### Essential dependencies
If you use GToTree, please be sure to cite these folks – a `citations.txt` file including used programs is produced with each run to help 🙂

- **[pyhmmer](https://pyhmmer.readthedocs.io/)** - [citation](https://academic.oup.com/bioinformatics/article/39/5/btad214/7131068) – HMM searching is done in-process through pyhmmer in v2, so a separate `hmmsearch` binary is not needed for the main workflow. Please also cite [HMMER3](http://hmmer.org/), which pyhmmer wraps.
- **[Muscle](https://www.drive5.com/muscle/downloads.htm)** v5.1 - [citation](https://academic.oup.com/nar/article/32/5/1792/2380623)  
- **[Trimal](http://trimal.cgenomics.org/downloads)** v1.4.1 - [citation](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2712344/)  
- **[FastTree](http://www.microbesonline.org/fasttree/)** v2.1.11 - [citation](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0009490)  
- **[Prodigal](https://github.com/hyattpd/Prodigal)** v2.6.3 - [citation](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2848648/)  
  - *needed whenever genes have to be called – i.e. input genomes in fasta format, GenBank format with no CDS annotations, or NCBI accessions to genomes with no gene calls*

### Optional dependencies depending on use
If you use GToTree in a manner that uses these tools, please cite these folks – a `citations.txt` file including used programs is produced with each run to help 🙂  

- **[IQ-TREE](http://www.iqtree.org/)** >=2.2.2 - [citation](https://academic.oup.com/mbe/article/37/5/1530/5721363)
  - *if using `-T IQTREE`*
- **[VeryFastTree](https://github.com/citiususc/veryfasttree)** - [citation](https://academic.oup.com/bioinformatics/article/36/17/4658/5843801)
  - *if using `-T VeryFastTree`*
- **[KOFamScan](https://github.com/takaram/kofam_scan)** v1.3.0 - [citation](https://academic.oup.com/bioinformatics/article/36/7/2251/5631907)
  - *if searching for target KOs (`-K`, or `gtt search-kos`)*
- **[HMMER3](http://hmmer.org/download.html)** - [citation](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002195)
  - *KOFamScan calls `hmmsearch` itself, so a HMMER3 binary is needed if you use the KO-searching features*
- **[Genome Taxonomy Database](https://www.nature.com/articles/s41587-020-0501-8)** - [citation](https://academic.oup.com/nar/article/50/D1/D785/6370255)  
  - *if adding GTDB taxonomy information (`-D`), or gathering genomes by GTDB taxonomy (`-w`)*
- **[Pfam](https://www.ebi.ac.uk/interpro/entry/pfam/)** - [citation](https://academic.oup.com/nar/article/51/D1/D418/6825349)
  - *if searching for target Pfams (`-p`, or `gtt search-pfams`), and cited for the packaged SCG-sets, which are Pfam-derived*

> **Note on versions**  
> The versions listed above are what the conda recipe pins or what was used at a point in GToTree's history, and are left here as a reference for anyone installing without `conda`. With the conda installation it can sometimes be better to be more flexible with regard to versions. You can check specific versions in your conda installation with `conda list`, and the `citations.txt` file produced by a GToTree run lists the versions of programs used for that specific run.

<!-- TODO(Mike): worth confirming this dependency split before publishing. It's derived
     from conda-recipe/meta.yaml plus which tools are actually invoked in the v2 code.
     Note the recipe currently also lists entrez-direct, dos2unix, parallel, bc,
     coreutils, sed, file, gzip, and curl — I couldn't find remaining uses of the first
     three in the v2 codebase, and the rest only back the small shell pipelines in
     citations.py. If those get pruned from the recipe, this section stays correct as-is.
     The recipe also lists pytest/pytest-cov under run requirements, which means every
     user installs the test stack. -->
