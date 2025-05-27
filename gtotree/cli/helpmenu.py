from gtotree.utils.messaging import color_text, gtotree_header

helpmenu = gtotree_header()

helpmenu += f"""

 ----------------------------------  {color_text("HELP INFO", "yellow")}  ----------------------------------

  This program takes input genomes from various sources and ultimately produces
  a phylogenomic tree. You can find detailed usage information at:
                                  github.com/AstrobioMike/GToTree/wiki


 -------------------------------  {color_text("REQUIRED INPUTS", "yellow")}  -------------------------------

      1) Input genomes in one or any combination of the following formats:
        - [-a <file>] single-column file of NCBI assembly accessions
        - [-g <file>] single-column file with the paths to each GenBank file
        - [-f <file>] single-column file with the paths to each fasta file
        - [-A <file>] single-column file with the paths to each amino acid file,
                      each file should hold the coding sequences for just one
                      genome

      2)  [-H <file>] location of the uncompressed target SCGs HMM file being
                      used, or just the HMM name if the 'GToTree_HMM_dir'
                      environment variable is set (which is typically handled
                      by the conda installation), run 'gtt-hmms' to view
                      available gene-sets


 ------------------------------  {color_text("OPTIONAL SETTINGS", "yellow")}  ------------------------------

      {color_text("Output directory specification:", "yellow")}

        - [-o <str>] default: gtotree-output


      {color_text("User-specified modification of genome labels:", "yellow")}

        - [-m <file>] mapping file specifying desired genome labels
                  A two- or three-column tab-delimited file where column 1 holds
                  either the input NCBI accession or input filename basename of
                  the genome to re-label (depending on the input source), column
                  2 holds the desired new genome label, and column 3 holds
                  something to be appended to either the initial or modified
                  labels (e.g. useful for "tagging" genomes in the tree based
                  on some characteristic). Columns 2 or 3 can be empty, and the
                  file does not need to include all input genomes.


      {color_text("Options for adding taxonomy information:", "yellow")}

        - [-D ] add GTDB taxonomy; default: false
                  Provide this flag with no arguments if you'd like to add
                  taxonomy from the Genome Taxonomy Database (GTDB;
                  gtdb.ecogenomic.org) to the sequence headers. This will
                  only be effective for input genomes provided as NCBI accessions
                  (provided to the `-a` argument). See the `-L` argument for
                  specifying desired ranks.

                  Note: See helper script `gtt-get-accessions-from-GTDB` for help
                  getting input NCBI accessions based on GTDB taxonomy searches.

        - [-t ] add NCBI taxonomy; default: false
                  Provide this flag with no arguments if you'd like to add NCBI
                  taxonomy info to the sequence headers. This will only be
                  effective for input genomes provided as NCBI accessions (provided
                  to the `-a` argument). See the `-L` argument for specifying
                  desired ranks.

        - [-L <str>] wanted lineage ranks; default: Domain,Phylum,Class,Genus,Species
                  A comma-separated list of the taxonomic ranks you'd like added
                  to the labels if adding taxonomic information. E.g., all would
                  be "-L Domain,Phylum,Class,Order,Family,Genus,Species,Strain".

                  Note: Strain-level information is available through NCBI, but
                  not GTDB.


      {color_text("Filtering settings:", "yellow")}

        - [-c <float>] sequence length cutoff; default: 0.2
                  A float between 0-1 specifying the range about the median of
                  sequences to be retained. For example, if the median length of
                  a set of sequences is 100 AAs, those seqs longer than 120 or
                  shorter than 80 will be filtered out before alignment of that
                  gene set with the default 0.2 setting.

        - [-G <float>] genome hits cutoff; default: 0.5
                  A float between 0-1 specifying the minimum fraction of hits a
                  genome must have of the target SCGs that end up in the final
                  tree. For example, if 100 target genes are going to comprise
                  the final tree, and Genome X only has hits to 49 of them, that
                  genome would be removed from the analysis with the default
                  value of 0.5 for this parameter.

                  Note: This is calculated on the total amount of SCGs that
                  contribute to the final tree, not necessarily the total
                  starting number of target SCGs.

        - [-B ] best-hit mode; default: false
                  Provide this flag with no arguments if you'd like to run
                  GToTree in "best-hit" mode. By default, if a SCG has more than
                  one hit in a given genome, GToTree won't include a sequence
                  for that target from that genome in the final alignment. With
                  this flag provided, GToTree will use the best hit. See here
                  for more discussion:
                      github.com/AstrobioMike/GToTree/wiki/things-to-consider


      {color_text("KO searching:", "yellow")}

        - [-K <file>] single-column file of KO targets to search for in each
                  genome. A table of hit counts, fastas of hit sequences, and
                  files compatible with the iToL web-based tree-viewer will be
                  generated for each target. You can see a visualization of
                  gene presence/absence example at:
                      github.com/AstrobioMike/GToTree/wiki/example-usage


      {color_text("Pfam searching:", "yellow")}

        - [-p <file>] single-column file of Pfam targets to search for in each
                  genome. A table of hit counts, fastas of hit sequences, and
                  files compatible with the iToL web-based tree-viewer will be
                  generated for each target. You can see a visualization of
                  gene presence/absence example at:
                      github.com/AstrobioMike/GToTree/wiki/example-usage


      {color_text("General run settings:", "yellow")}

        - [-R | --resume] resume mode; default: false
                  Provide this flag with no arguments if you'd like to try to resume a
                  previous run. This cannot be used if any inputs or options have changed.

        - [-z ] nucleotide mode; default: false
                  Make alignment and/or tree with nucleotide sequences instead
                  of amino-acid sequences.

                  Note: This mode can only accept NCBI accessions (passed to
                  `-a`) and genome fasta files (passed to `-f`) as input
                  sources. (GToTree still finds target genes based on
                  amino-acid HMM searches.)

        - [-N ] do not make a tree; default: false
                  No tree produced. Stop after producing the concatenated
                  alignment.

        - [-k ] keep individual target-gene alignments; default: false
                  Keep individual alignment files.

        - [-T <str>] tree program to use; default: FastTreeMP
                  Which program to use for tree generation. Currently supported
                  are "FastTree", "FastTreeMP", "VeryFastTree", and "IQTREE".
                  These run with default settings only (and IQTREE includes
                  "-m MFP" and "-B 1000"). To run any with more specific options
                  you can use the output alignment file from GToTree (and the
                  partitions file if wanted for mixed-model specification) as
                  input into a dedicated treeing program (the GToTree `-N`
                  option will generate the alignment only and skip internal
                  treeing if wanted).

        - [-j ] num jobs; default: 1
                  The number of jobs you'd like to run in parallel during steps
                  that are parallelizable. This includes things like downloading
                  input accession genomes and running parallel alignments, and
                  portions of the treeing step if using FastTreeMP or
                  VeryFastTree.

                  Note: I've occassionally noticed NCBI not being happy with
                  over ~50 downloads being attempted concurrently. So if using a
                  `-j` setting around there or higher, and GToTree is saying a
                  lot of input accessions were not successfully downloaded,
                  consider trying with a lower `-j` setting.

        - [-n <int> ] num HMM cpus; default: 2
                  The number of cpus you'd like to use during the HMM search.
                  Given these are individual small searches on single genomes,
                  2 is probably always sufficient. (Keep in mind this will be
                  multiplied by the number of jobs running concurrently if also
                  modifying the `-j` parameter.)

        - [-M <int> ] num muscle threads; default: 5
                  The number of threads muscle will use during alignment. (Keep
                  in mind this will be multiplied by the number of jobs running
                  concurrently if also modifying the `-j` parameter.)

        - [-X ] override super5 alignment; default: false
                  If working with greater than 1,000 target genomes, GToTree
                  will by default use the 'super5' muscle alignment algorithm
                  to increase the speed of the alignments. Provide this flag
                  with no arguments if you want to use the standard muscle
                  alignment.

                  Note: See sections in the link below for more information on
                  working with many genomes:
                      github.com/AstrobioMike/GToTree/wiki/things-to-consider

        - [-P (DEPRECATED)] use http instead of ftp (DEPRECATED)
                  Previously, GToTree used ftp where possible by default. Now it
                  uses http by default wherever it can, so this flag is no longer needed.

        - [-F ] force overwrite; default: false
                  Provide this flag with no arguments if you'd like to force
                  overwriting the output directory if it exists.

        - [--tmp-dir <path>] temporary directory location; default: <output-dir>/gtotree-tmp-*
                  If you want to specify where the temporary working directory will
                  be created, you can provide the path to this parameter.

        - [-d ] debug mode; default: false
                  Provide this flag with no arguments if you'd like to keep the
                  temporary directory.


 --------------------------------  {color_text("EXAMPLE USAGE", "yellow")}  --------------------------------

	GToTree -a ncbi-accessions.txt -f fasta-files.txt -H Bacteria -D -j 4

"""



