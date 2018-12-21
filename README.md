# GToTree: a user-friendly workflow for phylogenomics
The purpose of GToTree is to give more researchers the capability to create a full-scale Tree of Life, or any subset of it, to aid in their work.

---

## Currently under initial development - not yet released

GToTree is a more structured implementation of a workflow I would put together everytime I wanted to make a large-scale phylogenomic tree. What do I mean by large-scale? Anything from a full-blown Tree of Life with all 3 domains, down to, for example, all available genomes of *Staphylococcus*. At its heart it just takes in genomes and outputs an alignment and phylogenomic tree based on the specified HMM profiles. But I think its value comes from three main things: 1) its flexibility with regard to input format - taking fasta files, GenBank files, and/or NCBI accessions (So if you just recovered a bunch of new genomes and you want to see where they fit in with references, you can provide references by accession and your new genomes as fasta files.); 2) its automation of required in-between-tool tasks such as filtering hits by gene-length, filtering out genomes with too few hits to the target genes, and swapping genome labels for something more useful; and 3) its scalability – GToTree can turn 1,000 input genomes into a tree in ~90 minutes on a standard laptop.

Also included are several newly generated single-copy gene-sets for 13 different taxonomical groupings. These are presented in the [wiki](https://github.com/AstrobioMike/GToTree/wiki/SCG-sets), along with all code/steps used in the generation of them.  

GToTree utilizes helper scripts written in python, but is primarily implemented in bash. Every attempt is being made to make it portable across all variations, including Darwin on OSX, so please report an [issue](https://github.com/AstrobioMike/GToTree/issues) for any hiccups if you encounter them.  

---

See the [installation wiki](https://github.com/AstrobioMike/GToTree/wiki/installation) for help/suggestions with installing GToTree or any of the dependencies. *Spoiler alert:* [conda](https://conda.io/docs/) makes everything a cinch.  

---

## Dependencies
### Required to use at all:
If you use GToTree, please cite these folks :)  

- **[Biopython](https://biopython.org/wiki/Download)** - [citation](https://academic.oup.com/bioinformatics/article/25/11/1422/330687)
- **[HMMER3](http://hmmer.org/download.html)** - citation: they note in the [user manual](http://eddylab.org/software/hmmer/Userguide.pdf) to cite the website, but there is also [this paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002195)  
- **[Muscle](https://www.drive5.com/muscle/downloads.htm)** - [citation](https://academic.oup.com/nar/article/32/5/1792/2380623)  
- **[Trimal](http://trimal.cgenomics.org/downloads)** - [citation](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2712344/)  
- **[FastTree](http://www.microbesonline.org/fasttree/)** - [citation](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0009490)  

### Required based on inputs and options:
If you use GToTree in a manner that uses these tools, please cite these folks :)  

- **[Prodigal](https://github.com/hyattpd/Prodigal)** - [citation](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2848648/)  
  - *if providing input genomes in fasta format, or GenBank format with no CDS annotations, or NCBI accessions to genomes with no gene calls*
- **[TaxonKit](https://bioinf.shenwei.me/taxonkit/)** - citation: I'm not sure how to cite them yet
  - *if changing genome labels based on lineage information for input genomes with associated NCBI taxids*
- **[GNU Parallel](https://www.gnu.org/software/parallel/)** - [citation info](https://www.gnu.org/software/parallel/parallel_design.html#Citation-notice)
  - *if running in parallel (can make a big difference even on a laptop)*

## Testing instructions for those of you kind enough to be helping me squash da bugs
Make sure all dependencies above are installed and functioning properly. Clone the repo, add the "bin" to your PATH, and in the "test\_data" directory there are example files and an "example\_run\_log" code that should work. If it does not, please let me know what happened. If it does work, please start changing the inputs and use as you might to help find da bugs :)

If you have trouble getting E-Direct to work, it would still be a big help if you test things using just GenBank and/or fasta files as inputs (no accessions, no E-Direct needed). Thanks!

Running `GToTree` with no arguments gives the help menu. 

