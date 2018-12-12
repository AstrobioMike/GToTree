# GToTree

## Currently under development

GToTree is a command-line tool that can take any combination of flat fasta files, GenBank files, and NCBI assembly accessions as input and primarily outputs an alignment file, estimates of genome completeness and redundancy, and a phylogenomic tree based on the specified single-copy gene (SCG) set. GToTree enables automation of tasks such as: filtering out genes of spurious length and genomes with too few hits to the target SCG set via adjustable parameters; adding in needed gap-sequences for genes not detected in genomes being retained in the analysis; and swapping genome identifiers so resulting trees and alignments can be explored more easily.

GToTree utilizes helper scripts written in python and R, but is primarily implemented in bash. Every attempt is being made to make it portable across all variations, including Darwin on OSX, so please send in any hiccups if you encounter them.  

Dependencies include - if you use GToTree, please cite these folks :)  

- **HMMER3** - [http://hmmer.org/download.html](http://hmmer.org/download.html){:target="_blank"} (they note to cite the website)  
- **FastTree** - [http://www.microbesonline.org/fasttree/](http://www.microbesonline.org/fasttree/){:target="_blank"} (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0009490)  
- **prodigal** (*if providing input genomes in fasta format, or GenBank format with no CDS annotations*) - [https://github.com/hyattpd/Prodigal](https://github.com/hyattpd/){:target="_blank"} (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2848648/)  
- **NCBI E-Direct** (*if providing input genomes as NCBI assembly accessions*) - [https://www.ncbi.nlm.nih.gov/books/NBK179288/](https://www.ncbi.nlm.nih.gov/books/NBK179288/) (cite that page)
- **TaxonKit** (*if changing genome labels based on lineage information for input genomes with associated NCBI taxids*) - [https://bioinf.shenwei.me/taxonkit/](https://bioinf.shenwei.me/taxonkit/){:target="_blank"} (I'm not sure how to cite them yet)

Installation instructions - for now, as this is probably only going to be seen by Meren :)  
Clone the repo, and in the "test\_data" directory there are example files and an "example\_run\_log" code that works if the "bin" is in your PATH
