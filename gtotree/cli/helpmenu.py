from gtotree.utils.messaging import color_text, gtotree_header

helpmenu = gtotree_header()

helpmenu += f"""

 ----------------------------------  {color_text("HELP INFO", "yellow")}  ----------------------------------

  This program takes input genomes from various sources and ultimately produces
  a phylogenomic tree. You can find detailed usage information at:
                                  github.com/AstrobioMike/GToTree/wiki


 -------------------------------  {color_text("REQUIRED INPUTS", "yellow")}  -------------------------------

      1) Input genomes in one or any combination of the following formats:
        - [{color_text("-a <file>", "teal", bold = True)}] single-column file of NCBI assembly accessions
        - [{color_text("-g <file>", "teal", bold = True)}] single-column file with the paths to each GenBank file
        - [{color_text("-f <file>", "teal", bold = True)}] single-column file with the paths to each fasta file
        - [{color_text("-A <file>", "teal", bold = True)}] single-column file with the paths to each amino-acid file,
                      each file should hold the coding sequences for just one
                      genome

      2)  [{color_text("-H <file>", "teal", bold = True)}] location of the uncompressed target SCGs HMM file being
                      used, or just the HMM name if the 'GToTree_HMM_dir'
                      environment variable is set (which is typically handled
                      by the conda installation), run 'gtt-hmms' to view
                      available gene-sets


 ------------------------------  {color_text("OPTIONAL SETTINGS", "yellow")}  ------------------------------

      {color_text("Output directory specification:", "orange")}

        - [{color_text("-o <dir>", "teal", bold = True)}] default: gtotree-output


      {color_text("User-specified modification of genome labels:", "orange")}

        - [{color_text("-m <file>", "teal", bold = True)}] mapping file specifying desired genome labels
                  A two- or three-column tab-delimited file where column 1 holds
                  either the input NCBI accession or input filename basename of
                  the genome to re-label (depending on the input source), column
                  2 holds the desired new genome label, and column 3 holds
                  something to be appended to either the initial or modified
                  labels (e.g. useful for "tagging" genomes in the tree based
                  on some characteristic). Columns 2 or 3 can be empty, and the
                  file does not need to include all input genomes.


      {color_text("Options for adding taxonomy information:", "orange")}

        - [{color_text("-D | --add-gtdb-tax", "teal", bold = True)}] add GTDB taxonomy; default: false
                  Provide this flag with no arguments if you'd like to add
                  taxonomy from the Genome Taxonomy Database (GTDB;
                  gtdb.ecogenomic.org) to the sequence headers. This will
                  only be effective for input genomes provided as NCBI accessions
                  (passed to the `-a` parameter) that are in the GTDB. See the
                  `-L` argument for specifying desired ranks.

                  Note: You can use `gtt-get-accessions-from-GTDB` to get
                  input NCBI accessions based on GTDB taxonomy searches.

        - [{color_text("-t | --add-ncbi-tax", "teal", bold = True)}] add NCBI taxonomy; default: false
                  Provide this flag with no arguments if you'd like to add NCBI
                  taxonomy info to the sequence headers. This will only be
                  effective for input genomes provided as NCBI accessions (passed
                  to the `-a` parameter). See the `-L` argument for specifying
                  desired ranks.

        - [{color_text("-L <str>", "teal", bold = True)}] wanted lineage ranks; default: Domain,Phylum,Class,Genus,Species
                  A comma-separated list of the taxonomic ranks you'd like added
                  to the labels if adding taxonomic information. E.g., all would
                  be "-L Domain,Phylum,Class,Order,Family,Genus,Species,Strain".

                  Note: Strain-level information is available through NCBI, but
                  not GTDB.


      {color_text("Filtering settings:", "orange")}

        - [{color_text("-c <float>", "teal", bold = True)}] sequence-length cutoff; default: 0.2
                  A float between 0-1 (inclusive) specifying the range about the
                  median of sequences to be retained. For example, if the median
                  length of a target gene-set of sequences is 100 AAs, those seqs
                  longer than 120 or shorter than 80 will be filtered out before
                  alignment of that gene-set with the default 0.2 setting.

        - [{color_text("-r <float>", "teal", bold = True)}] gene-representation cutoff; default: 0.1
                  A float between 0-1 (inclusive) specifying the minimum proportion
                  of genomes that must have hits to a target gene for that gene to
                  be retained and used in the final tree. For example, if 100 input
                  genomes are provided, and target-gene X only has hits to 19 of them,
                  that target gene would be removed from the analysis with the default
                  value of 0.2 for this parameter.

                  Note: This is calculated based on the total number of genomes remaining
                  after preprocessing steps, but before genome filtering based on the '-G'
                  parameter, and it is not revisited after genome-level filtering. This is
                  to avoid iterative pruning effects between the '-g' and '-G' parameters.

        - [{color_text("-G <float>", "teal", bold = True)}] genome-hits cutoff; default: 0.5
                  A float between 0-1 (inclusive) specifying the minimum proportion
                  of target-gene hits a genome must have of the target SCGs that
                  end up in the final tree. For example, if 100 target genes are
                  going to comprise the final tree, and Genome X only has hits to
                  49 of them, that genome would be removed from the analysis with
                  the default value of 0.5 for this parameter.

                  Note: This is calculated based on the total amount of SCGs that will
                  contribute to the final tree (after any filtering), not necessarily
                  the total starting number of target SCGs.

        - [{color_text("-B | --best-hit-mode", "teal", bold = True)}] best-hit mode; default: false
                  Provide this flag with no arguments if you'd like to run
                  GToTree in "best-hit" mode. By default, if a SCG has more than
                  one hit in a given genome, GToTree won't include a sequence
                  for that target from that genome in the final alignment. With
                  this flag provided, GToTree will use the best hit. See here
                  for more discussion:
                      github.com/AstrobioMike/GToTree/wiki/things-to-consider


      {color_text("KO searching:", "orange")}

        - [{color_text("-K <file>", "teal", bold = True)}] single-column file of KO targets to search for in each
                  genome. A table of hit counts, fastas of hit sequences, and
                  files compatible with the iToL web-based tree-viewer will be
                  generated for each target. You can see a visualization of
                  gene presence/absence example at:
                      github.com/AstrobioMike/GToTree/wiki/example-usage


      {color_text("Pfam searching:", "orange")}

        - [{color_text("-p <file>", "teal", bold = True)}] single-column file of Pfam targets to search for in each
                  genome. A table of hit counts, fastas of hit sequences, and
                  files compatible with the iToL web-based tree-viewer will be
                  generated for each target. You can see a visualization of
                  gene presence/absence example at:
                      github.com/AstrobioMike/GToTree/wiki/example-usage


      {color_text("General run settings:", "orange")}

        - [{color_text("-j <int>", "teal", bold = True)}] num jobs; default: 4
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

        - [{color_text("-M <int>", "teal", bold = True)}] num muscle threads; default: 5
                  The number of threads muscle will use during alignment. (Keep
                  in mind this will be multiplied by the number of jobs running
                  concurrently if also modifying the `-j` parameter.)

        - [{color_text("-X | --no-super5", "teal", bold = True)}] override super5 alignment; default: false
                  If working with greater than 1,000 target genomes, GToTree
                  will by default use the 'super5' muscle alignment algorithm
                  to increase the speed of the alignments. Provide this flag
                  with no arguments if you want to use the standard muscle
                  alignment.

                  Note: See sections in the link below for more information on
                  working with many genomes:
                      github.com/AstrobioMike/GToTree/wiki/things-to-consider

        - [{color_text("-N | --no-tree", "teal", bold = True)}] do not make a tree; default: false
                  No tree produced. Stop after producing the concatenated
                  alignment.

        - [{color_text("-T <str>", "teal", bold = True)}] tree program to use; default: FastTreeMP
                  Which program to use for tree generation. Currently supported
                  are "{color_text("FastTree", "teal", bold = True)}", "{color_text("FastTreeMP", "teal", bold = True)}", "{color_text("VeryFastTree", "teal", bold = True)}", and "{color_text("IQTREE", "teal", bold = True)}".
                  These run with default settings only (and IQTREE includes
                  "-m MFP" and "-B 1000"). To run any with more specific options
                  you can use the output alignment file from GToTree (and the
                  partitions file if wanted for mixed-model specification) as
                  input into a dedicated treeing program (the GToTree `-N`
                  option will generate the alignment only and skip internal
                  treeing if wanted).

        - [{color_text("-z | --nucleotide-mode", "teal", bold = True)}] nucleotide mode; default: false
                  Make alignment and/or tree with nucleotide sequences instead
                  of amino-acid sequences. (GToTree still finds target genes
                  based on amino-acid HMM searches.)

                  Note: This mode can only accept NCBI accessions (passed to
                  `-a`) and genome fasta files (passed to `-f`) as input
                  sources.

        - [{color_text("-k | --keep-gene-alignments", "teal", bold = True)}] keep individual target-gene alignments; default: false
                  Keep individual alignment files.

        - [{color_text("-R | --resume", "teal", bold = True)}] resume mode; default: false
                  Provide this flag with no arguments if you'd like to try to resume a
                  previous run. This cannot be used if any inputs or options have changed.

        - [{color_text("-F | --force-overwrite", "teal", bold = True)}] force overwrite; default: false
                  Provide this flag with no arguments if you'd like to force
                  overwriting the output directory if it exists.

        - [{color_text("--tmp-dir <path>", "teal", bold = True)}] temporary directory location; default: <output-dir>/gtt-tmp-*
                  If you want to specify where the temporary working directory will
                  be created, you can provide the path to this parameter.

        - [{color_text("-d | --debug", "teal", bold = True)}] debug mode; default: false
                  Provide this flag with no arguments if you'd like to keep the
                  temporary directory.


 --------------------------------  {color_text("EXAMPLE USAGE", "yellow")}  --------------------------------

	GToTree -a ncbi-accessions.txt -f fasta-files.txt -H Bacteria -D

"""



