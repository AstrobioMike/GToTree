
<a href="https://github.com/AstrobioMike/AstrobioMike.github.io/blob/master/images/GToTree-logo-1200px.png"><img align="right" width="600" src="https://github.com/AstrobioMike/AstrobioMike.github.io/blob/master/images/GToTree-logo-1200px.png"></a>  


<br>  
<br>  
<br>
<br>
<br>

<a href="https://scholar.google.com/citations?view_op=view_citation&citation_for_view=-ONw6lsAAAAJ:_FxGoFyzp5QC"><img align="right" alt="Citations" src="https://img.shields.io/badge/Citations-650+-blue" height="22"></a>
<br>
<a href="https://github.com/AstrobioMike/GToTree/wiki/installation#conda-quickstart"><img align="right" alt="Conda installs" src="https://img.shields.io/badge/Conda%20installs-30k+-blue" height="22"></a>
<br>
<a href="https://doi.org/10.1093/bioinformatics/btz188"><img align="right" alt="DOI" src="https://img.shields.io/badge/DOI-10.1093/bioinformatics/btz188-blue" height="22"></a>
<br>
<a href="https://twitter.com/AstrobioMike"><img align="right" alt="Twitter Follow" src="https://img.shields.io/twitter/follow/AstrobioMike?color=blue&style=social"></a>
<br>

---

# Welcome to the GToTree wiki!

[GToTree](https://github.com/AstrobioMike/GToTree) is a user-friendly workflow for phylogenomics intended to give more researchers the capability to create phylogenomic trees and to make the process of iterating phylogenomic trees much easier. The open-access Bioinformatics Journal publication is available [here](https://doi.org/10.1093/bioinformatics/btz188). GToTree can be installed and run on a Mac or Linux machine, as well as on Windows within a Windows Subsystem for Linux environment 👍 

**For an overview of what this is all about, see the ["What is GToTree?" page](what-is-gtotree%3F)**. Or, to see some practical ways GToTree can be helpful, scan through the [Example usage page](example-usage).

---

**A quick [conda installation](https://github.com/AstrobioMike/GToTree/wiki/installation#conda-quickstart) can be run like so:**

```
conda create -y -n gtotree -c astrobiomike -c conda-forge -c bioconda gtotree
```

**And a tree can be as quick as naming the group you want:**

```
gtotree -w Nitrospirota -D -j 4 -o nitrospirota-tree
```

---

# Wiki contents

* [**What is GToTree?**](what-is-gtotree%3F)
* [**Installation**](installation)
  * [Conda Quickstart!](installation#conda-quickstart)
  * [Updating to a newer version](installation#updating-to-a-newer-version)
  * [Installation without Conda](installation#installation-without-conda)
* [**Example usage**](example-usage)
  * [Alteromonas example](example-usage#alteromonas-example)
  * [Building a tree from a taxon name](example-usage#building-a-tree-from-a-taxon-name)
  * [Using GToTree with the Genome Taxonomy Database (GTDB)](example-usage#using-gtotree-with-the-genome-taxonomy-database-gtdb)
  * [Visualization of gene-presence/absence across Bacteria](example-usage#visualization-of-gene-presenceabsence-across-the-bacterial-domain)
  * [Using the alignment and partitions file with another program](example-usage#using-the-alignment-and-partitions-file-with-another-program)
  * [Tree-of-Life example](example-usage#tol-example)
* [**User Guide**](user-guide)
  * [Required inputs](user-guide#required-inputs)
  * [Outputs](user-guide#outputs)
  * [Optional arguments and parameters](user-guide#optional-arguments-and-parameters)
  * [Resuming an interrupted run](user-guide#resuming-an-interrupted-run)
  * [Options set for programs run](user-guide#options-set-for-programs-run)
  * [All programs used citation info](user-guide#citation-information)
* [**Single-copy gene-sets**](scg-sets)
  * [The original (v1) SCG-sets](scg-sets#the-original-v1-scg-sets)
* [**Things to consider**](things-to-consider)
  * [An important caveat on the idea of a workflow for phylogenomics](things-to-consider#an-important-caveat-on-the-idea-of-a-workflow-for-phylogenomics)
  * [When to use GToTree and when not](things-to-consider#when-to-use-gtotree-and-when-not)
    * [Is GToTree useful for assigning taxonomy?](things-to-consider#is-gtotree-useful-for-assigning-taxonomy)
    * [Should GToTree be used for estimating genome/MAG/bin completeness?](things-to-consider#should-gtotree-be-used-for-estimating-genomemagbin-quality)
  * [Consider using "representative" genomes](things-to-consider#consider-using-representative-genomes)
  * [Working with many genomes](things-to-consider#working-with-many-genomes)
  * [Filtering genes by length](things-to-consider#filtering-hits-by-gene-length)
  * [Filtering genomes by hits to targets](things-to-consider#filtering-genomes-by-fraction-of-hits-to-targets)
  * [Best-hit mode](things-to-consider#best-hit-mode)
  * [Eukaryotes with GToTree?](things-to-consider#are-eukaryotic-genomes-appropriate-for-use-with-gtotree)

---

<p align="center">
<a href="https://github.com/AstrobioMike/AstrobioMike.github.io/blob/master/images/GToTree-Overview-main.png"><img src="https://github.com/AstrobioMike/AstrobioMike.github.io/blob/master/images/GToTree-Overview-main.png"></a>
</p>