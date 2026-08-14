This page serves as the general user guide for GToTree, which may be helpful if you're looking for something specific. To jump right into practical ways GToTree can be helpful it may be more useful to start with the [Example-usage page](https://github.com/AstrobioMike/GToTree/wiki/example-usage) :) 

--- 

# User-Guide Contents

* [**Required inputs**](user-guide#required-inputs)
  * [**Input genomes**](user-guide#input-genomes)
  * [**Specifying which single-copy gene-set to use**](user-guide#specifying-which-single-copy-gene-set-to-use)
* [**Outputs**](user-guide#outputs)
  * [**Primary output files**](user-guide#primary-output-files)
  * [**Report output files**](user-guide#report-output-files)
* [**Optional arguments and parameters**](user-guide#optional-arguments-and-parameters)
* [**Resuming an interrupted run**](user-guide#resuming-an-interrupted-run)
* [**Options set for programs run**](user-guide#options-set-for-programs-run)
* [**All programs used citation info**](user-guide#citation-information)

---

> **NOTE:** Running `gtotree -h` will provide the condensed help menu, and `gtotree -s` the detailed one with every parameter.

> **NOTE for v1 users:** several things changed in v2 — helper programs are now subcommands of `gtt`, some output files were renamed, and a few parameters were dropped. See the [CHANGELOG](https://github.com/AstrobioMike/GToTree/blob/master/CHANGELOG.md) for the full list.

---

# Required Inputs
The minimum required inputs to GToTree are specifying the genomes you want to incorporate and specifying which single-copy gene-set to use.

Genomes can be provided via any combination of a target taxon (`-w`), NCBI accessions (`-a`), GenBank files (`-g`), and/or nucleotide (`-f`) or amino-acid (`-A`) fasta files.

## Input Genomes

### Reference genomes by taxonomy
New in v2, you can ask GToTree to gather reference genomes for a taxon rather than collecting accessions yourself. Pass a taxon name to the **`-w`** argument:

```bash
gtotree -w Nitrospirota -o nitrospirota-tree
```

* The taxon can be at any rank (e.g. `-w Bacteria`, `-w Alteromonas`, `-w "Escherichia coli"`).
* **`-w` may be given more than once** to pool several taxa (e.g. `-w Bacteria -w Archaea`). Each is resolved and dereplicated on its own, then merged.
* **`--source`** controls which taxonomy the genomes are drawn from, `gtdb` (default) or `ncbi`.
* **`--target-rank`** disambiguates a name that exists at more than one rank.
* **`--derep-rank`** dereplicates down to one genome per unique value of a rank. The default, `auto`, uses two ranks finer than the target; `off` keeps every genome under the taxon.

This can be combined freely with the other input types, so you can add your own new genomes as fasta files alongside references pulled by taxonomy.

### NCBI Accessions
You can specify which NCBI-archived genomes you'd like to incorporate by providing a single-column file holding NCBI assembly accessions to the **`-a`** argument. This file can be created "manually" by searching NCBI's website and downloading a results table, or generated at the command line with `gtt get-accs-from-ncbi` or `gtt get-accs-from-gtdb` – examples are on the [examples page here](https://github.com/AstrobioMike/GToTree/wiki/example-usage#accessions).

* Those provided can have version numbers (what comes after the "." in the accession, e.g. GCF_000153765.1), or they can be version-less (e.g. GCF_000153765). In the case where no version is provided, GToTree will automatically take the newest released version of that accession.
* If any of the provided accessions cannot be found at NCBI, they will be reported in `run-files/removed-genomes.tsv`.

### GenBank files
To specify which GenBank files to include, you need to provide a single-column file that holds the file names (or [paths](https://astrobiomike.github.io/unix/basics#absolute-vs-relative-path)) to each of the GenBank files you'd like incorporated. This is passed to the **`-g`** argument.

### Fasta files
Nucleotide fasta files are provided similarly to the GenBank files, but passed to the **`-f`** argument. Genes will be called with [Prodigal](https://github.com/hyattpd/Prodigal).

### Amino-acid files
Amino-acid files can be provided similarly, but passed to the **`-A`** argument. Each file should hold the proteins for one genome.

## Specifying which single-copy gene-set to use
GToTree needs to know which SCG-set to use – passed with the **`-H`** flag. There are 46 packaged with the program (discussed in more detail [here](https://github.com/AstrobioMike/GToTree/wiki/scg-sets)). You can view which are available by running **`gtt hmms`**, and you don't need to specify a full path, just the name as printed there (e.g. `-H Bacteria`). You can also pass a path to your own HMM file.

> **`-H` is optional when `-w` is used.** If you provide a target taxon and don't specify `-H`, GToTree selects an appropriate packaged set automatically: the lowest-rank set covering all the taxa you asked for. Taxa under one phylum get that phylum's set, taxa spanning both domains get `Bacteria-and-Archaea`, and anything reaching outside Bacteria and Archaea falls back to `Universal-Hug-et-al`. The set chosen and the reason for it are printed in the run's opening banner. Passing `-H` always overrides this.

# Outputs
Each GToTree run creates an output directory to hold all of the output files. This defaults to "gtotree-output", but can be specified with the **`-o`** argument.

### Primary output files

#### Tree
* **&lt;output-dir-name&gt;.tre**
  * The final tree file in newick format, midpoint-rooted. Named after the output directory, so `-o alteromonas-tree` produces `alteromonas-tree/alteromonas-tree.tre`.
  * The pre-midpoint-rooted tree is kept in `run-files/`.
  * FastTree reports "local support values" that appear as labels on internal nodes to estimate the reliability of each split in the tree. You can find more information about this at their user page [here](http://www.microbesonline.org/fasttree/#Support).
  * IQ-TREE reports ultrafast bootstrap (UFBoot) support values. Their [help pages](http://www.iqtree.org/doc/Frequently-Asked-Questions) state that values of 95% indicate a 95% probability that clade is true.
  * If run with the `-N` option, no tree will be produced, and only the alignment will be generated.

#### Alignment files
* **aligned-SCGs.faa** (or `aligned-SCGs.fasta` in nucleotide mode)
  * Alignment file in fasta format.
* **aligned-SCGs-mod-names.faa**
  * Produced if labels were modified – i.e. if a mapping file was given (`-m`) or taxonomy was added (`-t` or `-D`).

#### Genomes summary info
* **genomes-summary-info.tsv**
  * A tab-delimited table of summary information for each genome:

|Column|Contents|
|:---|---|
| genome_id | the input genome ID (either the accession or base file name depending on input source) |
| input | the input the genome came from (e.g. the accessions file, or the `-w` taxon) |
| source | which input type the genome came from |
| label | the label assigned to the genome in the output tree file |
| taxid | the NCBI taxid, where one is available |
| num_SCG_hits | number of gene hits to the target HMMs |
| num_uniq_SCG_hits | number of unique gene hits to the target HMMs |
| num_SCG_hits_after_filtering | number of gene hits remaining after length filtering |
| num_total_genes | total genes for the genome |
| in_final_tree | Yes or No, did this genome end up in the final tree |
| reason_removed | if it didn't, why |

> If target Pfams (`-p`) or KOs (`-K`) were searched, `pfam_search_completed` / `ko_search_completed` columns are added. If taxonomy was added (`-t` or `-D`), a column per taxonomic rank is added.

#### Citations
* **citations.txt**
  * Citation information for every program actually used in that specific run, with versions. See [below](user-guide#citation-information).

#### Log
* **gtotree-runlog.txt**
  * A record of what was printed to the screen during the run.

### Report output files
These live in the output sub-directory `run-files/`.

**removed-genomes.tsv**
  * Every input genome that didn't make it into the final tree, at whatever stage it was dropped, in one table. Columns are `genome_id`, `input`, `source`, `stage_removed`, and `reason_removed`.
  * In v1 this information was spread across several files (`NCBI_accessions_not_found.txt`, `Genomes_removed_for_too_few_hits.tsv`, and others); it is all consolidated here now, with the `stage_removed` column telling you which case each row is.

**SCG-info.tsv**
  * One row per target SCG, with how many genomes had hits to it at each stage, whether it was retained, and if not, why.

**partitions.txt** and **partitions.nex**
  * Partition files compatible with treeing programs capable of using different models for each gene. See, for example, [iqtree's info here](http://www.iqtree.org/doc/Advanced-Tutorial).

**run-data.json**
  * The run's internal state. This is what `-R`/`--resume` reads to pick a run back up.

---

## Optional arguments and parameters
The condensed help menu can be viewed by running `gtotree -h`, and the detailed one, with full descriptions of every parameter, with `gtotree -s`.

---

### Output directory
* **[-o \<str\>] default: gtotree-output**

GToTree writes all output files to an output directory. By default this is "gtotree-output", but you can specify it by passing an argument to the **`-o`** flag. (E.g.: `-o alteromonas-output`)

---

### Overwrite an existing output directory
* **[-F ] default: false**

By default GToTree won't write into an output directory that already exists. Provide **`-F`**/`--force-overwrite` to delete and recreate it.

---

### Specify desired genome labels
* **[-m \<file\>] specify desired genome labels**

Often it is helpful to have specific labels for specific genomes in a tree (as exemplified in the [Alteromonas example](https://github.com/AstrobioMike/GToTree/wiki/example-usage#mapping-file-for-labeling-specific-genomes)). GToTree adds lineage information to genomes that have it if we specify we want NCBI taxonomies added (with the `-t` flag) or GTDB taxonomy (with the `-D` flag). We can also swap labels of specific genomes we know we care about and want to be able to find more easily, or just append certain information to the label.

Either or both of these can be done by providing a mapping file to the **`-m`** argument. It should be a 2- or 3-column tab-delimited file that has the initial genome ID in the first column (this will be either the NCBI accession or the file name, depending on how the input genome was provided). The second column may or may not be empty. If we want to specify the complete label ourselves for that genome, then we should put that new label in column 2. If we don't want to specify the complete label, leave column 2 empty. Column 3 may or may not be empty. If we'd like to append something to the label (whether that's the initial label, the modified lineage label, or the label we may have specified in column 2), then add that text to column 3. If there is nothing we want to append, we should leave column 3 empty.

>**NOTE:** Not all input genomes need to be provided in the file being passed to `-m`.

---

### Specify to add NCBI lineage info to genome labels
* **[-t ] default: false**

By setting the **`-t`** flag, GToTree will look up NCBI lineage information for any genomes that have an associated taxid and add it to the genome labels – making the output tree much more useful than just a collection of odd identifiers. Which specific taxonomic ranks get added can be specified with the **`-L`** argument.

> **Changed in v2:** taxids are no longer pulled out of input GenBank files; they are used only for genomes provided as NCBI accessions.

---

### Specify to add [GTDB](https://gtdb.ecogenomic.org/) lineage info to genome labels
* **[-D ] default: false**

By setting the **`-D`** flag, GToTree will add GTDB taxonomy information to the labels that appear in the final alignment and in the tree. Which specific taxonomic ranks get added can be specified with the **`-L`** argument. If a given accession isn't present in GTDB, the NCBI lineage info will be used (and "_NCBI" will be appended to it).

> `-D` and `-t` can't both be given; pick one.

---

### Specify which taxonomic ranks to add to genome labels
* **[-L \<str\>] default: domain,phylum,class,genus,species**

Provide the **`-t`** or **`-D`** flag in order to add lineage info to the genome labels. By default this adds domain, phylum, class, genus, and species, where available. This may be suitable when making a tree across multiple domains, but may be unnecessarily cumbersome when just making a tree of one genus. You can specify which ranks you'd like added with the **`-L`** argument as a comma-separated list, e.g. `-L domain,phylum,class,order,family,genus,species,strain`.

> **Changed in v2:** the default is now `domain,phylum,class,genus,species` — genus is included and strain is not. Rank names are case-insensitive, so v1-style `-L Domain,Phylum,Class` still works.

---

### Filtering gene-hits by length
* **[-c \<float\>] default: 0.2**

When scanning many genomes for many genes, it becomes harder or completely impractical to visually inspect alignments of everything. One way to try to filter out potential spurious gene hits is to filter by some expected length. The **`-c`** parameter uses the median length of each particular gene-set to calculate an upper- and lower-length threshold to filter out potentially spurious genes. It takes a float between 0-1 specifying the range about the median of sequences to be retained. For example, under the default setting, if the median length of a set of sequences is 100 AAs, those genes with sequences longer than 120 or shorter than 80 will be filtered out before alignment of that gene set. This becomes less useful when using very few genomes however (see [note here](https://github.com/AstrobioMike/GToTree/wiki/things-to-consider#filtering-hits-by-gene-length)).

---

### Filtering genes by how many genomes have them
* **[-r \<float\>] default: 0.1**

New in v2. The **`-r`** parameter sets the minimum proportion of genomes that must have a hit to a target gene for that gene to be retained. It takes a float between 0-1. For example, with 100 input genomes and the default of 0.1, a target gene with hits in only 9 of them would be dropped from the analysis.

This is calculated on the genomes remaining after processing but *before* genome-level filtering with `-G`, and it is not revisited afterward — that avoids `-r` and `-G` iteratively pruning each other.

---

### Filtering genomes based on hits to target genes
* **[-G \<float\>] default: 0.5**

The **`-G`** parameter allows you to filter out genomes that have too few hits to the target genes. It takes a float between 0-1 specifying the minimum fraction of hits a genome must have of the SCG-set. For example, under the default setting, if 100 target genes will make up the final tree, and genome X only has hits to 49 of them, it will be removed from analysis. How you want this set may depend on the breadth of diversity of the tree you are making (see [note here](https://github.com/AstrobioMike/GToTree/wiki/things-to-consider#filtering-genomes-by-fraction-of-hits-to-targets)).

> **Changed in v2:** this is now calculated against the number of target genes that will actually contribute to the final tree (i.e. after `-r` filtering), rather than against the starting total. So if you start with 100 target genes and 20 are dropped by `-r`, the fraction is taken against 80.

---

### Best-hit mode
* **[-B ] default: false**

Provide the **`-B`** flag with no arguments if you'd like to run GToTree in "best-hit" mode. By default, if a target gene has more than one hit in a given genome, GToTree won't include a sequence for that target gene from that genome in the final alignment. With this flag provided, GToTree will take the best hit and incorporate it into the alignment, even if that genome has more than one hit to the target gene. See [here](https://github.com/AstrobioMike/GToTree/wiki/Things-to-consider#best-hit-mode) for more discussion on this.

---

### Searching for additional target genes
* **[-p \<file\>] target Pfams**
* **[-K \<file\>] target KOs**

Single-column files of Pfam accessions (e.g. `PF00789`) or KO IDs (e.g. `K01601`) to search the input genomes for alongside the SCGs. Results are written to `pfam-search-results/` and `ko-search-results/` in the output directory, including hit-count tables and files for visualization in [iToL](https://itol.embl.de/).

> These same searches are available as standalone subcommands, `gtt search-pfams` and `gtt search-kos`, if you want them without building a tree.

---

### Number of jobs to run in parallel where possible
* **[-j \<int\>] default: 4**

This determines how many jobs to run in parallel during steps that are parallelizable – such as the processing/searching of each individual genome, the filtering of genes and genomes, and alignment of each individual gene-set.

> **Changed in v2:** the default is now 4 (was 1). The separate `-n` parameter for HMM-search cpus has been dropped.

---

### Number of threads passed to muscle
* **[-M \<int\>] default: 5**

How many threads each individual muscle alignment job is given.

---

### Overriding the super5 alignment algorithm
* **[-X ] default: false**

With more than 1,000 input genomes, GToTree uses muscle's [super5](https://drive5.com/muscle5/manual/super5_algo.html) algorithm to keep alignment times reasonable. Provide **`-X`** to use the regular [PPP](https://drive5.com/muscle5/manual/ppp_algo.html) algorithm anyway – see [here](https://github.com/AstrobioMike/GToTree/wiki/things-to-consider#working-with-many-genomes) for more.

---

### Choosing the treeing program
* **[-T \<str\>] default: FastTreeMP**

Which program builds the tree. Options are `FastTree`, `FastTreeMP`, `VeryFastTree`, and `IQTREE`.

---

### Skipping the tree
* **[-N ] default: false**

Provide **`-N`** to stop after producing the alignment, without building a tree.

---

### Keeping individual gene alignments
* **[-k ] default: false**

Provide **`-k`** to keep the individual per-gene alignments rather than only the concatenated one.

---

### Keep working directory
* **[-d ] default: false**

Provide the **`-d`** flag with no arguments if you'd like to keep the working directory that is used during the run. This is mostly useful for debugging purposes.

---

### Temporary directory location
* **[--tmp-dir \<path\>] default: &lt;output-dir&gt;/gtt-tmp-\***

Where the temporary working directory is created, if you need it somewhere other than inside the output directory (e.g. on a faster or larger disk).

---

### Generate nucleotide alignment and tree
* **[-z ] default: false**

Provide the **`-z`** flag with no arguments if you'd like to make the alignment and tree based on nucleotide sequences instead of amino-acids (which can provide greater resolution on closely related input genomes). Note this mode can only accept NCBI accessions (passed to `-a`) and genome fasta files (passed to `-f`) as input sources (since we can't confidently reverse-translate amino-acid seqs). GToTree still finds the target genes based on amino-acid HMM searches.

---

# Resuming an interrupted run

New in v2. If a run is interrupted – a walltime limit, a dropped connection partway through downloading genomes, a killed process – you can pick it back up by re-running the same command with **`-R`**/`--resume` added:

```bash
gtotree -w Nitrospirota -o nitrospirota-tree -R
```

GToTree fingerprints the inputs and settings of each run. On a resume, if anything that would change the result has been altered (a different SCG-set, a modified input file, a different cutoff), the run stops and tells you which parameter changed rather than silently mixing results from two different analyses. Settings that only affect how the run executes, like `-j`, can be changed freely.

`-R` and `-F` can't be used together, since one reuses the previous run and the other deletes it.

---

# Options set for programs run
> **NOTE**  
> In order to give conda more freedom in managing its environments, specific versions have been removed from the conda installation. If you installed with conda, and want to know the specific versions, you can check in the environment with `conda list`.

## prodigal
[Prodigal](https://github.com/hyattpd/Prodigal) is run with default settings other than setting the `-c` flag, which means only include complete genes (and `-q` to quiet its output).

## HMM searching
HMM searches are performed in-process with [pyhmmer](https://pyhmmer.readthedocs.io/), rather than by shelling out to `hmmsearch`. The searches use the gathering thresholds (`--cut_ga`-equivalent) stored in the HMM profiles being used for cutoff values, as before. GToTree will exit up front if a provided HMM file has profiles lacking gathering thresholds.

## muscle
[Muscle](https://www.drive5.com/muscle/downloads.htm) is run with default settings using the [`-align` PPP algorithm](https://drive5.com/muscle5/manual/ppp_algo.html) when working with fewer than 1,000 target genomes. When run with greater than 1,000 input genomes, the [`-super5` algorithm](https://drive5.com/muscle5/manual/super5_algo.html) is used. This can be overridden by adding `-X` to the GToTree call – see [here](https://github.com/AstrobioMike/GToTree/wiki/things-to-consider#working-with-many-genomes) for more on this. The number of threads given to each alignment job is set with `-M` (default 5).

## trimal
[Trimal](http://trimal.cgenomics.org/downloads) is run with default settings other than setting the `-automated1` flag, which performs "a heuristic selection of the automatic method based on similarity statistics. (Optimized for Maximum Likelihood phylogenetic tree reconstruction)."

## FastTree
[FastTree](http://www.microbesonline.org/fasttree/) and FastTreeMP are run with default settings for amino-acid alignments, and with the flags `-nt` and `-gtr` for nucleotide alignments. FastTreeMP is given `OMP_NUM_THREADS` matching `-j`.

## VeryFastTree
[VeryFastTree](https://github.com/citiususc/veryfasttree) is run with `-threads` matching `-j`, and with `-nt -gtr` for nucleotide alignments.

## IQ-TREE
[IQ-TREE](http://www.iqtree.org/) is run with `-nt` matching `-j`, `-m MFP` for model selection, and `-B 1000` for ultrafast bootstraps. Bootstrapping is skipped with fewer than 4 genomes, as IQ-TREE requires at least that many.

---

# Citation information

GToTree relies on many great programs. Along with all other outputs, it will generate a `citations.txt` file with citation information specific for every run that accounts for all programs it relies upon. Please be sure to cite the developers appropriately :)

Here is an example output `citations.txt` file from a run, and how I'd cite it in the methods:

```
GToTree v2.0.0
Lee MD. GToTree: a user-friendly workflow for phylogenomics. Bioinformatics. 2019; (March):1-3. doi:10.1093/bioinformatics/btz188

Prodigal v2.6.3
Hyatt, D. et al. Gene and translation initiation site prediction in metagenomic sequences. Bioinformatics. 2010; 28, 2223–2230. doi.org/10.1186/1471-2105-11-119

HMMER (via pyhmmer v0.11.0)
Eddy SR. Accelerated profile HMM searches. PLoS Comput. Biol. 2011; (7)10. doi:10.1371/journal.pcbi.1002195

Muscle v5.1
Edgar RC. MUSCLE v5 enables improved estimates of phylogenetic tree confidence by ensemble bootstrapping. bioRxiv. 2021. doi.org/10.1101/2021.06.20.449169

TrimAl v1.4.rev15
Gutierrez SC. et al. TrimAl: a Tool for automatic alignment trimming. Bioinformatics. 2009; 25, 1972–1973. doi:10.1093/bioinformatics/btp348

FastTree 2 v2.1.11
Price MN et al. FastTree 2 - approximately maximum-likelihood trees for large alignments. PLoS One. 2010; 5. doi:10.1371/journal.pone.0009490
```

**Example methods text based on above citation output (be sure to modify as appropriate for your run)**
> *The archaeal phylogenomic tree was produced with GToTree v2.0.0 (Lee 2019), using the prepackaged single-copy gene-set for archaea (51 target genes). Briefly, prodigal v2.6.3 (Hyatt et al. 2010) was used to predict genes on input genomes provided as fasta files. Target genes were identified with HMMER3 (Eddy 2011) via pyhmmer v0.11.0, individually aligned with muscle v5.1 (Edgar 2021), trimmed with trimal v1.4.rev15 (Capella-Gutiérrez et al. 2009), and concatenated prior to phylogenetic estimation with FastTree2 v2.1.11 (Price et al. 2010).*
