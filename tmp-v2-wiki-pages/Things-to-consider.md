# An important caveat on the idea of a “workflow for phylogenomics”
Phylogenetics is an incredibly complicated and well-researched field, and things become even more complicated when working with many concatenated genes as is the case with phylogenomics. GToTree is meant to be a relatively high-throughput, user-friendly, and reproducible workflow, something I believe is useful due to the high volumes of sequencing data and genomes we are often working with these days. But anything designed this way needs to inherently sacrifice something in terms of flexibility, options, precision, etc. It is important that users new to this arena understand that many things impact the outcome of a phylogenetic/genomic analysis, particularly including the alignment algorithm used, and the model and program used for tree construction. Currently, GToTree employs only one alignment tool, and four options for tree construction. Users can also take the concatenated alignment output by GToTree (and the [partitions](example-usage#using-the-alignment-and-partitions-file-with-another-program) file if they'd like) and use that with many other tree construction tools. But please keep in mind that phylogenetic analysis is complicated, and no one program or tool is an "absolute answer" or the "truth" – another way to think of this is with the old adage "all models are wrong".

---

# When to use GToTree and when not?
GToTree is very useful in many situations, but not all. For example, it is useful if you want to make a large-scale Tree of Life spanning all 3 domains including a lot of genomes ([like demonstrated here](https://github.com/AstrobioMike/GToTree/wiki/example-usage#tol-example)). And it is useful if you want to infer evolutionary relationships between some newly recovered genomes and references on a smaller scale ([like the *Alteromonas* example here](https://github.com/AstrobioMike/GToTree/wiki/example-usage#alteromonas-example)). But even if you use a specific marker-gene set to make a tree of all the organisms of interest (like the Gammaproteobacteria set in the Alteromonas example), this is only useful at the level of resolution those marker genes provide. Often that may be enough for your purposes, but sometimes you might need or want to go deeper. In cases like this, you may want to use GToTree to figure out where your new genomes fit in with, say, 500 reference genomes, and then you could use that tree to identify which reference genomes you actually want to include in a pangenomic analysis with your new genomes.

## Is GToTree useful for assigning taxonomy?
No, GToTree is not a tool for taxonomic assignment. For assigning taxonomy to genomes there are dedicated programs that work very well like the **G**enome **T**axonomy **D**ata**b**ase **T**ool**k**it ([GTDB-Tk](https://github.com/Ecogenomics/GTDBTk)). I would typically assign taxonomy to my new genomes with GTDB-Tk, and then I would use that information to figure out which reference genomes I'd want to include in a de novo phylogenomic tree I'd build with GToTree. You can however make a de novo tree, look at where your new genomes fall on that tree, and see which references they are more closely related to based on that tree. For example, with the [*Alteromonas* example here](https://github.com/AstrobioMike/GToTree/wiki/example-usage#alteromonas-example), we start there "knowing" the new genome is an *Alteromonas*, and we build a de novo tree with all RefSeq reference *Alteromonas* genomes and our new one. Ahead of where that example starts, I may have figured out that my new genome was an *Altermonas* by using [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk).

## Should GToTree be used for estimating genome/MAG/bin quality?
No. Earlier versions of GToTree reported a rough estimated completion/redundancy based on the SCG-set being used, since it was basically free information given what we're already looking for. That was removed in v2 – it was a rudimentary view of something that dedicated tools do properly, and it invited more misreading than it was worth. For estimating genome/MAG/bin quality, I'd recommend [CheckM2](https://github.com/chklovski/CheckM2) 👍 

The `num_SCG_hits` and `num_uniq_SCG_hits` columns in `genomes-summary-info.tsv` still tell you how many target genes were found in each genome, which is what the `-G` filter acts on, but treat those as a description of the run rather than a quality metric.


---

# Consider using "representative" genomes
It is not often that we want to cover a large breadth of diversity while **also** needing multiple very closely related genomes in our tree. If we are trying to cover a large breadth of diversity, it can be very helpful to use "representative" reference genomes – a slimmed set of manually and computationally derived genomes that are designed to capture the breadth of microbial diversity by choosing representatives in genomic lineage space. I regularly make use of both [NCBI's "reference" genomes](https://www.ncbi.nlm.nih.gov/refseq/about/prokaryotes/#representative_genomes) and [GTDB's species representatives](https://gtdb.ecogenomic.org/faq#how-are-gtdb-species-clusters-formed). 

For example, if we wanted to build a tree of the Staphylococcaceae family, we can first check how many genomes that is:

```bash
gtt get-accs-from-gtdb -t Staphylococcaceae --get-taxon-counts
```

<!-- TODO(Mike): paste current output here -->

That is generally way too many to practically make a phylogenomic tree out of. Limiting to [GTDB species representatives](https://gtdb.ecogenomic.org/faq#gtdb_species_clusters) with `-G` cuts it down substantially:

```bash
gtt get-accs-from-gtdb -t Staphylococcaceae --get-taxon-counts -G
```

<!-- TODO(Mike): paste current output here -->

Which is much more manageable, while still covering the breadth of diversity within the family, just with fewer closely related genomes.

`-R`/`--refseq-reference-genomes-only` does the analogous thing with [NCBI's reference genomes](https://www.ncbi.nlm.nih.gov/refseq/about/prokaryotes/#representative_genomes).

## Dereplicating by rank

Representative genomes get us a long way, but sometimes they are still more than we want – especially when making a tree that spans an entire domain. For that, GToTree can dereplicate a selection down to one genome per unique value of a taxonomic rank, with `--derep-rank`.

This is available directly in the main call:

```bash
  # one genome per order across the whole archaeal domain
gtotree -w Archaea --derep-rank order -D -j 4 -o archaea-by-order
```

and in the accession-gathering helper, if you want the list rather than the tree:

```bash
gtt get-accs-from-gtdb -t Archaea --derep-rank order
```

One per order is generally fine if the goal is showing where new MAGs fit into a domain-level tree, in my experience.

When `-w` is used, `--derep-rank` defaults to `auto`, which dereplicates to two ranks finer than whatever you targeted – so `-w Archaea` (a domain) gives one genome per class. Set it explicitly to go coarser or finer, or use `--derep-rank off` to keep every genome under the taxon.

> **Note for v1 users:** this replaces the old `gtt-subset-GTDB-accessions` helper, which no longer exists. `--derep-rank` does the same job, and does it inside the main run rather than as a separate step on the metadata table.

---

# Working with many genomes
When working with many thousands of genomes, the time the individual gene alignments take can quickly become prohibitive. First, I'd suggest considering working with [representative genomes](#consider-using-representative-genomes) as discussed above, if not doing so already. But if we are targeting greater than 1,000 genomes, GToTree by default will use the [super5](https://drive5.com/muscle5/manual/super5_algo.html) algorithm of muscle to help speed up the alignments. A message will pop up notifying you when the run starts, but if you'd like to use the regular [PPP](https://drive5.com/muscle5/manual/ppp_algo.html) muscle algorithm still despite having many sequences to align, then pass the `-X` option with no arguments to the GToTree call.

---

# Filtering hits by gene-length
The default setting for this value (set with **`-c`**) is 0.2. This value will filter out genes that are 20% larger or smaller in length than the median of all hits to a target single-copy gene. For example, this means if the median length of all genes selected as best hits to marker-gene "A" is 100, genes that were hits to marker-gene "A" that are greater in length than 120 or shorter in length than 80 will be removed from the analysis. This seems to work well in my experience, but only when there are enough genes in the gene set to give somewhat of a representative distribution of the lengths of genes that exist within that target gene set. Meaning, at the extreme end, if we only had 3 genes to consider, and their lengths were 100, 100, and 121, the 121-length gene would be filtered out, but maybe it shouldn't be. If running GToTree with very few genomes, you might consider increasing this threshold and/or visually inspecting some of the alignments.

---

# Filtering genomes by fraction of hits to targets
The default for this value (set with **`-G`**) is set to 0.5, meaning that if 100 target genes are going into the final tree, genomes with hits to fewer than 50 of them will be dropped from the analysis. This seems to me to be reasonable in most cases. But if you are seeing a lot of genomes dropped due to this filter, you could consider lowering it (the [Tree of Life example](example-usage#tol-example) uses `-G 0.4`, since a tree spanning both domains will naturally have more genomes missing more targets).

> **Changed in v2:** the fraction is now taken against the number of target genes that will actually contribute to the final tree – i.e. after genes have been filtered by `-r`/`--gene-rep-cutoff` – rather than against the starting number in the SCG-set. So if you start with 100 targets and 20 are dropped for having too few hits across genomes, `-G 0.5` means needing hits to 40 of the remaining 80, not 50 of the original 100.

---

# Best-hit mode
By default, if a given genome has more than one hit to a specific HMM profile (target gene), GToTree **won't** include a sequence for that target gene from that genome in the final concatenated alignment (it will insert a gap-sequence just as would be the case if that genome had 0 hits to the target). This is a conservative way to go, because if there are multiple copies of a target SCG present within a genome, the copies may not all be under the same evolutionary pressures, and which one we choose may impact the alignment and tree in ways we do not want it to. So I figure in general, being conservative is better for default settings. But if you'd like, you can specify the **`-B`** flag with no arguments to tell GToTree to run in "best-hit" mode. In this case, when a given genome has more than one hit to a specific target gene, GToTree will take the best hit and add it to the alignment.

---

# Are eukaryotic genomes appropriate for use with GToTree?
If only using highly conserved ribosomal proteins (like those in the `Universal-Hug-et-al` SCG-set), and/or if all genes are already identified (e.g. the input source is an NCBI accession with gene calls or a GenBank file with gene calls), then GToTree is suitable for working with Eukaryotes in addition to Bacteria and Archaea. If no gene-calls are available, then GToTree is likely not suitable for eukaryotic genomes as the only gene-caller currently implemented is prodigal.

Note also that the 45 pre-built SCG-sets other than `Universal-Hug-et-al` are built from GTDB, which covers only Bacteria and Archaea. If you gather genomes by taxonomy (`-w`) and the target reaches into eukaryotes, auto-selection will fall back to the universal set for that reason.

