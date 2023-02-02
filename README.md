
<img align="right" width="600" src="https://github.com/AstrobioMike/AstrobioMike.github.io/blob/master/images/GToTree-logo-1200px.png">  

<br>
<br>
<br>
<br>

<a href="https://scholar.google.com/citations?view_op=view_citation&citation_for_view=-ONw6lsAAAAJ:_FxGoFyzp5QC"><img align="right" alt="Citations" src="https://img.shields.io/badge/Citations-190+-blue" height="22"></a>
<br>
<a href="https://github.com/AstrobioMike/GToTree/wiki/installation#conda-quickstart"><img align="right" alt="Conda installs" src="https://img.shields.io/badge/Conda%20installs-2.4k-blue" height="22"></a>
<br>
<a href="https://doi.org/10.1093/bioinformatics/btz188"><img align="right" alt="DOI" src="https://img.shields.io/badge/DOI-10.1093/bioinformatics/btz188-blue" height="22"></a>
<br>
<a href="https://twitter.com/AstrobioMike"><img align="right" alt="Twitter Follow" src="https://img.shields.io/twitter/follow/AstrobioMike?color=blue&style=social"></a>
<br>

---

# GToTree: a user-friendly workflow for phylogenomics
[GToTree](https://github.com/AstrobioMike/GToTree/wiki) is a user-friendly workflow for phylogenomics intended to give more researchers the capability to easily create phylogenomic trees. Documentation and examples can be found [at the wiki here](https://github.com/AstrobioMike/GToTree/wiki), and the open-access Bioinformatics Journal publication is available [here](https://doi.org/10.1093/bioinformatics/btz188).

---

**See the [conda quickstart](https://github.com/AstrobioMike/GToTree/wiki/installation#conda-quickstart) installation page to have things up and running in just a single step!**

```
conda create -y -n gtotree -c conda-forge -c bioconda -c astrobiomike gtotree
```

---

GToTree is a more structured implementation of a workflow I would put together everytime I wanted to make a large-scale phylogenomic tree. What do I mean by large-scale? Anything from a full-blown Tree of Life with all 3 domains, down to, for example, all available genomes of *Staphylococcus* alongside new isolate genomes. At its heart it just takes in genomes and outputs an alignment and phylogenomic tree based on the specified HMM profiles. But I think its value comes from three main things: 1) its flexibility with regard to input format - taking fasta files, GenBank files, and/or NCBI accessions (So if you just recovered a bunch of new genomes and you want to see where they fit in with references, you can provide references by accession and your new genomes as fasta files.); 2) its automation of required between-tool tasks such as filtering hits by gene-length, filtering out genomes with too few hits to the target genes, and swapping genome labels for something more useful; and 3) its scalability – GToTree can turn ~1,700 input genomes into a tree in ~60 minutes on a standard laptop.

Also included are several newly generated single-copy gene-sets for 13 different taxonomical groupings. These are presented in the [wiki](https://github.com/AstrobioMike/GToTree/wiki/SCG-sets), along with an explanation and example code/steps used in the generation of them. 

GToTree utilizes helper scripts written in python, but is primarily implemented in bash. Every attempt is being made to make it portable across all variations of GNU/Unix, including on Macs, so if you run into any issues, it'd be appreciated if you could [report them](https://github.com/AstrobioMike/GToTree/issues) so the problems can be found and fixed!  

<p align="center">
<a href="https://github.com/AstrobioMike/AstrobioMike.github.io/blob/master/images/GToTree-Overview-main.png"><img src="https://github.com/AstrobioMike/AstrobioMike.github.io/blob/master/images/GToTree-Overview.png"></a>
</p>

See the ["What is GToTree?" wiki page](https://github.com/AstrobioMike/GToTree/wiki/what-is-gtotree%3F) for some more detail on the processing steps pictured above. For practical ways GToTree can be helpful, check out the [Example usage page](https://github.com/AstrobioMike/GToTree/wiki/example-usage). And for detailed information on using GToTree, see the [User guide](https://github.com/AstrobioMike/GToTree/wiki/user-guide).

---

**See the [conda quickstart](https://github.com/AstrobioMike/GToTree/wiki/installation#conda-quickstart) installation page to get GToTree up and running in just a single step!**

```
conda create -y -n gtotree -c conda-forge -c bioconda -c astrobiomike gtotree
```

---

## Citation information

GToTree will print out a `citations.txt` file with citation information specific for every run that accounts for all programs it relies upon. Please be sure to cite the developers appropriately :)

Here is an example output `citations.txt` file from a run, and how I'd cite it in the methods:

```
GToTree v1.6.31
Lee MD. GToTree: a user-friendly workflow for phylogenomics. Bioinformatics. 2019; (March):1-3. doi:10.1093/bioinformatics/btz188

Prodigal v2.6.3
Hyatt, D. et al. Gene and translation initiation site prediction in metagenomic sequences. Bioinformatics. 2010; 28, 2223–2230. doi.org/10.1186/1471-2105-11-119

HMMER3 v3.3.2
Eddy SR. Accelerated profile HMM searches. PLoS Comput. Biol. 2011; (7)10. doi:10.1371/journal.pcbi.1002195

Muscle v5.1
Edgar RC. MUSCLE v5 enables improved estimates of phylogenetic tree confidence by ensemble bootstrapping. bioRxiv. 2021. doi.org/10.1101/2021.06.20.449169

TrimAl v1.4.rev15
Gutierrez SC. et al. TrimAl: a Tool for automatic alignment trimming. Bioinformatics. 2009; 25, 1972–1973. doi:10.1093/bioinformatics/btp348

TaxonKit v0.9.0
Shen W and Ren H. TaxonKit: a practical and efficient NCBI Taxonomy toolkit. Journal of Genetics and Genomics. 2021. doi.org/10.1016/j.jgg.2021.03.006

FastTree 2 v2.1.11
Price MN et al. FastTree 2 - approximately maximum-likelihood trees for large alignments. PLoS One. 2010; 5. doi:10.1371/journal.pone.0009490
```

**Example methods text based on above citation output (be sure to modify as appropriate for your run)**
> *The archaeal phylogenomic tree was produced with GToTree v1.6.31 (Lee 2019), using the prepackaged single-copy gene-set for archaea (76 target genes). Briefly, prodigal v2.6.3 (Hyatt et al. 2010) was used to predict genes on input genomes provided as fasta files. Target genes were identified with HMMER3 v3.2.2 (Eddy 2011), individually aligned with muscle v5.1 (Edgar 2021), trimmed with trimal v1.4.rev15 (Capella-Gutiérrez et al. 2009), and concatenated prior to phylogenetic estimation with FastTree2 v2.1.11 (Price et al. 2010). TaxonKit (Shen and Ren 2021) was used to connect full lineages to taxonomic IDs.*
