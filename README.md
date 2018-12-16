# GToTree: a user-friendly workflow for large-scale phylogenomics

## Currently under initial development - not yet ready for use

GToTree is a more structured implementation of a workflow I would put together everytime I wanted to make a large-scale phylogenomic tree. What do I mean by large-scale? Anything from a full-blown Tree of Life with all 3 domains, down to, for example, all available genomes of Staphylococcus. At its heart it just takes in genomes and outputs a phylogenomic tree based on the specified HMM profiles. But its value comes from its flexibility with regard to input format, enabling the automation of required in-between-tool tasks (such as filtering hits by gene-length, filtering out genomes with too few hits to the target genes, and swapping genome labels for something more useful), and its scalability – GToTree can turn 1,000 input genomes into a tree in about 90 minutes on a standard laptop. 

Also included is a newly generated bacterial single-copy gene-set based on Pfams from release 32.0 searched against 11,405 complete bacterial genomes available from NCBI (accessed on 09-DEC-18). All code/steps used in the generation of the gene set are reported in the [wiki](https://github.com/AstrobioMike/GToTree/wiki/Bacterial-SCG-set-generation).

GToTree utilizes helper scripts written in python and R, but is primarily implemented in bash. Every attempt is being made to make it portable across all variations, including Darwin on OSX, so please report an issue for any hiccups if you encounter them.  

## Dependencies
### Required to use at all:
If you use GToTree, please cite these folks :)  

- **HMMER3** - [http://hmmer.org/download.html](http://hmmer.org/download.html) (citation: they note to cite the website)  
- **Muscle** - [https://www.drive5.com/muscle/downloads.htm](https://www.drive5.com/muscle/downloads.htm) (citation: [https://academic.oup.com/nar/article/32/5/1792/2380623](https://academic.oup.com/nar/article/32/5/1792/2380623))
- **Trimal** - [http://trimal.cgenomics.org/downloads](http://trimal.cgenomics.org/downloads) (citation: [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2712344/](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2712344/))
- **FastTree** - [http://www.microbesonline.org/fasttree/](http://www.microbesonline.org/fasttree/) (citation: [https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0009490](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0009490))  

### Required based on inputs and options:
If you use GToTree in a manner that uses these tools, please cite these folks :)  

- **prodigal** (*if providing input genomes in fasta format, or GenBank format with no CDS annotations, or NCBI accessions to genomes with no gene calls*) - [https://github.com/hyattpd/Prodigal](https://github.com/hyattpd/) (citation: [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2848648/](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2848648/))  
- **NCBI E-Direct** (*if providing input genomes as NCBI assembly accessions*) - [https://www.ncbi.nlm.nih.gov/books/NBK179288/](https://www.ncbi.nlm.nih.gov/books/NBK179288/) (citation: cite that page)
- **TaxonKit** (*if changing genome labels based on lineage information for input genomes with associated NCBI taxids*) - [https://bioinf.shenwei.me/taxonkit/](https://bioinf.shenwei.me/taxonkit/) (citation: I'm not sure how to cite them yet)

## Testing instructions for those of you kind enough to be helping me squash da bugs
Make sure all dependencies above are installed and functioning properly. Clone the repo, add the "bin" to your PATH, and in the "test\_data" directory there are example files and an "example\_run\_log" code that should work. If it does not, please let me know what happened. If it does work, please start changing the inputs and use as you might to help try to find da bugs :)

Running `GToTree` with no arguments gives the help menu. 
