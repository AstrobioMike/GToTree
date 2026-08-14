There are 46 single-copy gene (SCG)-set HMMs packaged with GToTree. 45 of them were
built from the [GTDB](https://gtdb.ecogenomic.org/) taxonomy (release r232) as described
below, and one covers all 3 domains (from
[Hug et al. 2016](https://www.nature.com/articles/nmicrobiol201648)).

We can view which are available by running **`gtt hmms`**.

> **Note for GToTree v1 users:** the packaged sets changed substantially in v2. The v1
> sets are still available for download — see
> [The original (v1) SCG-sets](#the-original-v1-scg-sets) below.

# Available SCG-sets

### Bacteria

|SCG-set|Rank|Target genes|Genomes used|Dereplicated to|
|----|:----:|:----:|:----:|:----:|
|Bacteria|domain|52|5146|family|
|Acidobacteriota|phylum|97|985|genus|
|Actinomycetota|phylum|92|1993|genus|
|Aquificota|phylum|194|114|species|
|Armatimonadota|phylum|119|201|genus|
|Babelota|phylum|108|112|genus|
|Bacillota|phylum|79|3257|genus|
|Bacillota_I|phylum|71|513|genus|
|Bacteroidota|phylum|118|2462|genus|
|Bdellovibrionota|phylum|93|245|genus|
|Campylobacterota|phylum|182|120|genus|
|Chlamydiota|phylum|119|189|genus|
|Chloroflexota|phylum|69|1061|genus|
|Cyanobacteriota|phylum|119|525|genus|
|Deinococcota|phylum|132|55|genus|
|Desulfobacterota|phylum|101|918|genus|
|Elusimicrobiota|phylum|107|112|genus|
|Fibrobacterota|phylum|114|52|genus|
|Fidelibacterota|phylum|116|89|genus|
|Fusobacteriota|phylum|129|142|species|
|Gemmatimonadota|phylum|129|267|genus|
|Myxococcota|phylum|86|129|genus|
|Myxococcota_A|phylum|83|406|genus|
|Nitrospinota|phylum|217|57|species|
|Nitrospirota|phylum|134|244|genus|
|Omnitrophota|phylum|130|366|genus|
|Patescibacteriota|phylum|53|3454|genus|
|Planctomycetota|phylum|100|1186|genus|
|Pseudomonadota|phylum|97|4932|genus|
|Spirochaetota|phylum|74|338|genus|
|Synergistota|phylum|126|62|genus|
|Thermotogota|phylum|161|117|species|
|Verrucomicrobiota|phylum|95|694|genus|
|Zixibacteria|phylum|149|107|genus|
|Alphaproteobacteria|class|106|2343|genus|
|Gammaproteobacteria|class|108|2550|genus|

### Archaea

|SCG-set|Rank|Target genes|Genomes used|Dereplicated to|
|----|:----:|:----:|:----:|:----:|
|Archaea|domain|51|492|family|
|Asgardarchaeota|phylum|87|90|species|
|Halobacteriota|phylum|112|225|genus|
|Methanobacteriota|phylum|221|224|species|
|Micrarchaeota|phylum|76|93|genus|
|Nanobdellota|phylum|53|421|genus|
|Thermoplasmatota|phylum|93|180|genus|
|Thermoproteota|phylum|74|338|genus|

### Spanning more than one domain

|SCG-set|Scope|Target genes|Source|
|----|:----:|:----:|----|
|Bacteria-and-Archaea|multi-domain|30|https://doi.org/10.1093/bioinformatics/btz188|
|Universal-Hug-et-al|universal|16|www.nature.com/articles/nmicrobiol201648|

More information on each — including the specific Pfam accessions in every set — is in the
`hmm-sources-and-info.tsv` file
[here](https://github.com/AstrobioMike/GToTree/blob/master/hmm_sets/hmm-sources-and-info.tsv),
and in the GToTree installed location (see `gtt hmms`).

# How these were generated

All 45 GTDB-based sets were built with `gtt gen-scg-hmms`, the same helper program that
ships with GToTree for building your own (run `gtt gen-scg-hmms -h`). In outline:

1. All Pfams from release 38.2 were downloaded. Those that didn't cover more than 50% (on
   average) of the underlying protein sequences that built that Pfam's HMM were filtered
   out, so no two Pfams from the same source protein could both be included.

2. For each target taxon, the genomes classified to it in GTDB r232 were gathered and
   dereplicated (to one per genus for most sets; see the "Dereplicated to" column above),
   which keeps a handful of heavily sequenced genera from dominating the counts.

3. The filtered Pfams were searched against all coding sequences of those genomes using
   the gathering thresholds curated in the Pfam models, and hits per Pfam per genome
   were counted.

4. Only Pfams with exactly 1 hit in at least 90% of the genomes for that taxon were kept
   as that taxon's SCG-set.

Sets were built for each domain, and for taxa with enough genomes to support the above.

## Building your own

`gtt gen-scg-hmms` does all of the above for any target group. The simplest form takes a
taxon directly:

```bash
gtt gen-scg-hmms -w Nitrospirota -o my-nitrospirota-scgs
```

It accepts the same genome inputs the main program does (`-a`, `-f`, `-A`, `-g`), so you
can also build a set from your own genomes. See `gtt gen-scg-hmms -h` for the filtering
options (`--min-completeness`, `--max-contamination`, `--percent-single-copy`,
`--min-pfam-coverage`).

# The original (v1) SCG-sets

The sets packaged with GToTree v1 were built in 2018 against NCBI taxonomy and Pfam 32.0.
They are no longer packaged with GToTree, but they are all still available as assets on
[this release](https://github.com/AstrobioMike/GToTree/releases/tag/original-SCG-HMMs):

Actinobacteria, Alphaproteobacteria, Archaea, Bacteria, Bacteria_and_Archaea,
Bacteroidetes, Betaproteobacteria, Chlamydiae, Cyanobacteria, Epsilonproteobacteria,
Firmicutes, Gammaproteobacteria, Proteobacteria, Tenericutes, and Universal-Hug-et-al.

To use one, download it and pass the file path to `-H`:

```bash
curl -LO https://github.com/AstrobioMike/GToTree/releases/download/original-SCG-HMMs/Firmicutes.hmm

gtotree -a accessions.txt -H Firmicutes.hmm -o my-output
```

> **Note:** `Bacteria` and `Archaea` exist in both the old and new sets, but they are not
> the same sets. The packaged v2 `-H Bacteria` has 52 target genes where v1's had 74, and
> `-H Archaea` has 51 where v1's had 76. If you are extending or reproducing an earlier
> v1 analysis, download the original set from the release above rather than relying on
> the name.
