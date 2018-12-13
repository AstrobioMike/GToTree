# GToTree

## Currently under development

GToTree is a command-line tool that can take any combination of flat fasta files, GenBank files, and NCBI assembly accessions as input and primarily outputs an alignment file, estimates of genome completeness and redundancy, and a phylogenomic tree based on the specified single-copy gene (SCG) set. GToTree enables automation of tasks such as: filtering out genes of spurious length and genomes with too few hits to the target SCG set via adjustable parameters; adding in needed gap-sequences for genes not detected in genomes being retained in the analysis; and swapping genome identifiers so resulting trees and alignments can be explored more easily.

Also included is a newly generated bacterial single-copy gene-set based on Pfams from release 32.0 searched against 11,405 complete bacterial genomes available from NCBI (accessed on 09-DEC-18). All code/steps used in the generation of the gene set are reported in the [wiki](https://github.com/AstrobioMike/GToTree/wiki).

GToTree utilizes helper scripts written in python and R, but is primarily implemented in bash. Every attempt is being made to make it portable across all variations, including Darwin on OSX, so please send in any hiccups if you encounter them.  

Dependencies include - if you use GToTree, please cite these folks :)  

- **HMMER3** - [http://hmmer.org/download.html](http://hmmer.org/download.html) (citation: they note to cite the website)  
- **Muscle** - [https://www.drive5.com/muscle/downloads.htm](https://www.drive5.com/muscle/downloads.htm) (citation: [https://academic.oup.com/nar/article/32/5/1792/2380623](https://academic.oup.com/nar/article/32/5/1792/2380623))
- **FastTree** - [http://www.microbesonline.org/fasttree/](http://www.microbesonline.org/fasttree/) (citation: [https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0009490](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0009490))  
- **prodigal** (*if providing input genomes in fasta format, or GenBank format with no CDS annotations*) - [https://github.com/hyattpd/Prodigal](https://github.com/hyattpd/) (citation: [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2848648/](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2848648/))  
- **NCBI E-Direct** (*if providing input genomes as NCBI assembly accessions*) - [https://www.ncbi.nlm.nih.gov/books/NBK179288/](https://www.ncbi.nlm.nih.gov/books/NBK179288/) (citation: cite that page)
- **TaxonKit** (*if changing genome labels based on lineage information for input genomes with associated NCBI taxids*) - [https://bioinf.shenwei.me/taxonkit/](https://bioinf.shenwei.me/taxonkit/) (citation: I'm not sure how to cite them yet)

Testing instructions  
Clone the repo, and in the "test\_data" directory there are example files and an "example\_run\_log" code that works if the "bin" is in your PATH
