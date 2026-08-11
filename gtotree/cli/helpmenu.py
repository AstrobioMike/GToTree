import sys
import textwrap

from gtotree.utils.misc.messaging import color_text, gtotree_header


WIDTH = 82

SUBSECTION_INDENT = 6
ENTRY_INDENT = 8
DETAIL_INDENT = 18


################################################################################
# rendering helpers
################################################################################

def banner(title):
    """The ' ------  TITLE  ------ ' section rules from the original menu."""
    inner = len(title) + 4          # two spaces of padding either side
    dashes = WIDTH - 1 - inner
    left = dashes // 2 + dashes % 2
    right = dashes // 2
    return f"\n\n {'-' * left}  {color_text(title, 'yellow')}  {'-' * right}\n"


def subsection(title, extra_gap=False):
    lead = "\n" if extra_gap else ""
    return f"{lead}\n{' ' * SUBSECTION_INDENT}{color_text(title, 'orange')}\n"


def flag(text):
    return color_text(text, "teal", bold=True)


def bold(text):
    if not text or not sys.stdout.isatty():
        return text
    return f"\033[1m{text}\033[0m"


CONDENSED_NOTE = [
    "This is the condensed help menu. Run with [-S | --show-detailed-help] to",
    "see all available parameters and their full descriptions.",
]


def condensed_note():
    out = "\n"
    for line in CONDENSED_NOTE:
        pad = " " * max(0, (WIDTH - len(line)) // 2)
        if "[-S | --show-detailed-help]" in line:
            before, after = line.split("[-S | --show-detailed-help]")
            rendered = bold(before + "[") + flag("-S | --show-detailed-help") + bold("]" + after)
        else:
            rendered = bold(line)
        out += pad + rendered + "\n"
    return out


def pick(text, detailed):
    """
    A summary line can be either a plain string (same wording in both menus) or a
    (condensed, detailed) pair when the two tiers should read differently
    """
    if isinstance(text, tuple):
        condensed, full = text
        return full if detailed else condensed
    return text


def entry(flags, brief, detailed=False):
    """A '        - [-o <dir>] one-line summary' row, wrapped to the detail column."""
    brief = pick(brief, detailed)
    lead = f"{' ' * ENTRY_INDENT}- [{flag(flags)}] "
    # the color codes don't occupy visible columns, so measure the plain version
    plain_lead_len = ENTRY_INDENT + 3 + len(flags) + 2
    wrapped = textwrap.fill(
        " ".join(brief.split()),
        width=WIDTH,
        initial_indent=" " * plain_lead_len,
        subsequent_indent=" " * DETAIL_INDENT,
    )
    return lead + wrapped[plain_lead_len:] + "\n"


def detail_block(text):
    text = textwrap.dedent(text.strip("\n"))
    out = []
    for para in text.split("\n\n"):
        lines = []
        buf = []
        for ln in para.split("\n"):
            if ln.startswith("    "):        # preserving pre-indented lines (like URLs)
                if buf:
                    lines.append((" ".join(buf), 0))
                    buf = []
                lines.append((ln.strip(), len(ln) - len(ln.lstrip())))
                continue
            buf.append(ln.strip())
        if buf:
            lines.append((" ".join(buf), 0))

        rendered = []
        for ln, extra in lines:
            rendered.append(textwrap.fill(
                ln,
                width=WIDTH,
                initial_indent=" " * (DETAIL_INDENT + extra),
                subsequent_indent=" " * (DETAIL_INDENT + extra),
                break_on_hyphens=False,
            ))
        out.append("\n".join(rendered))
    return "\n" + "\n\n".join(out) + "\n"


################################################################################
# the one table both menus render from
################################################################################
# Each entry: (flags, brief, detail, brief_tier)
#   brief_tier True  -> shown in both menus
#   brief_tier False -> only appears under -S

BRIEF, FULL = True, False

ENTRIES = [
    ("Output directory settings:", [
        ("-o | --output-dir <dir>",
         'default: gtotree-output',
         None,
         BRIEF),

        ("-R | --resume",
         "resume mode",
         """
         Provide this flag if you'd like to try to resume a previous
         run. This cannot be used if any inputs or options have changed.
         """,
         BRIEF),

        ("-F | --force-overwrite",
         "force overwrite output dir",
         """
         Provide this flag if you'd like to force overwriting the
         output directory if it exists.
         """,
         BRIEF),
    ]),

    ("Adding reference genomes by taxonomy:", [
        ("-w | --wanted-ref-tax <str>",
         "ref tax to add genomes for; default: none",
         """
         A taxon name (e.g., "Bacteria", "Nitrospirota", "Escherichia coli") whose
         reference genomes will be added to the tree in addition to any genomes
         provided through the other input parameters.
         """,
         FULL),

        ("--source <str>",
         "where to pull reference genomes from; default: gtdb",
         """
         Which database the `-w` taxon and its genomes are drawn from. One of
         "gtdb" or "ncbi".
         """,
         FULL),

        ("--target-rank <str>",
         "rank of the `-w` taxon; default: auto-detected",
         """
         Only needed to disambiguate if the `-w` target taxon name occurs at more
         than one rank.
         """,
         FULL),

        ("--derep-rank <str>",
         "dereplicate to one genome per rank; default: auto",
         """
         Keeps a single best genome per unique value of this rank within the `-w`
         taxon's rank (to control tree size/complexity). For example, "-w Bacteria
         --derep-rank class" keeps one genome per bacterial class. The default
         "auto" uses a rank two levels finer than the taxon's own rank. Pass "none"
         to disable dereplication and include all genomes under the requested taxon.

         Note: a `--derep-rank` set equal to the `-w` taxon's rank returns a single
         best genome (e.g., could be useful for adding an outgroup).
         """,
         FULL),
    ]),

    ("User-specified modification of genome labels:", [
        ("-m | --mapping-file <file>",
         "mapping file specifying desired genome labels",
         """
         A two- or three-column tab-delimited file where column 1 holds either the
         input NCBI accession or input filename basename of the genome to re-label
         (depending on the input source), column 2 holds the desired new genome
         label, and column 3 holds something to be appended to either the initial or
         modified labels (e.g., useful for tagging genomes in the tree based on some
         characteristic). Columns 2 or 3 can be empty, and the file does not need to
         include all input genomes.
         """,
         FULL),
    ]),

    ("Options for adding taxonomy information to tree/alignment:", [
        ("-D | --add-gtdb-tax",
         "add GTDB taxonomy",
         """
         Provide this flag if you'd like to add GTDB taxonomy info
         to the sequence headers for any genomes included based on
         `-w` or `-a` inputs. See the `--lineage-ranks` argument
         for specifying desired ranks.

         Note: You can use the `-w` parameter here to include reference genomes
         based on GTDB taxonomy searches.
         """,
         BRIEF),

        ("-t | --add-ncbi-tax",
         "add NCBI taxonomy",
         """
         Provide this flag if you'd like to add NCBI taxonomy info
         to the sequence headers for any genomes included based on
         `-w` or `-a` inputs. See the `--lineage-ranks` argument
         for specifying desired ranks.

         Note: You can use the `-w` and `--source` parameters here to include reference
         genomes based on NCBI taxonomy searches.
         """,
         BRIEF),

        ("-L | --lineage-ranks <str>",
         "wanted lineage ranks; default: domain,phylum,class,genus,species",
         """
         A comma-separated list of the taxonomic ranks you'd like added to the labels
         if adding taxonomic information. E.g., all would be "--lineage-ranks
         domain,phylum,class,order,family,genus,species,strain".

         Note: Strain-level only applicable to NCBI taxonomy (not GTDB).
         """,
         FULL),
    ]),

    ("Filtering settings:", [
        ("-c | --seq-length-cutoff <float>",
         "sequence-length cutoff; default: 0.2",
         """
         A float between 0-1 (inclusive) specifying the range about the median of
         sequences to be retained. For example, if the median length of a target
         gene-set of sequences is 100 AAs, those seqs longer than 120 or shorter than
         80 will be filtered out before alignment of that gene-set with the default
         0.2 setting.
         """,
         FULL),

        ("-r | --gene-rep-cutoff <float>",
         "gene-representation cutoff; def: 0.1",
         """
         A float between 0-1 (inclusive) specifying the minimum proportion of genomes
         that must have hits to a target gene for that gene to be retained and used in
         the final tree. For example, if 100 input genomes are provided, and
         target-gene X only has hits to 9 of them, that target gene would be removed
         from the analysis with the default value of 0.1 for this parameter.

         Note: This is calculated based on the total number of genomes remaining after
         the primary genome retrieval and processing steps, but before genome
         filtering based on the `--genome-hits-cutoff` parameter, and it is not
         revisited after genome-level filtering. This is to avoid iterative pruning
         effects between the `--gene-rep-cutoff` and `--genome-hits-cutoff`
         parameters.
         """,
         FULL),

        ("-G | --genome-hits-cutoff <float>",
         "genome-hits cutoff; default: 0.5",
         """
         A float between 0-1 (inclusive) specifying the minimum proportion of
         target-gene hits a genome must have of the target SCGs that end up in the
         final tree. For example, if 100 target genes are going to comprise the final
         tree, and Genome X only has hits to 49 of them, that genome would be removed
         from the analysis with the default value of 0.5 for this parameter.

         Note: This is calculated based on the total amount of SCGs that will
         contribute to the final tree (after any filtering), not necessarily the total
         starting number of target SCGs.
         """,
         FULL),

        ("-B | --best-hit-mode",
         "best-hit mode",
         """
         Provide this flag if you'd like to run GToTree in "best-hit" mode.
         By default, if a SCG has more than one hit in a given genome, GToTree
         won't include a sequence for that target from that genome in the final
         alignment. With this flag provided, GToTree will allow multiple hits and
         will just use the best hit. See here for more discussion:
             github.com/AstrobioMike/GToTree/wiki/things-to-consider
         """,
         FULL),
    ]),

    ("KO searching:", [
        ("-K | --target-kos <file>",
         "single-column file of KO targets to search for in each genome",
         """
         A table of hit counts, fastas of hit sequences, and files compatible with the
         iToL web-based tree-viewer will be generated for each target KO.
         """,
         FULL),
    ]),

    ("Pfam searching:", [
        ("-p | --target-pfams <file>",
         "single-column file of Pfam targets to search for in each genome",
         """
         A table of hit counts, fastas of hit sequences, and files compatible with the
         iToL web-based tree-viewer will be generated for each target Pfam.
         """,
         FULL),
    ]),

    ("General run settings:", [
        ("-j | --num-jobs <int>",
         "num jobs; default: 4",
         """
         The number of jobs you'd like to run in parallel during steps that are
         parallelizable. This includes things like downloading input genomes,
         running alignments, and portions of the treeing step if using
         FastTreeMP or VeryFastTree.
         """,
         BRIEF),

        ("-T | --tree-program <str>",
         "tree program to use; default: FastTreeMP",
         """
         Which program to use for tree generation. Currently supported are
         "FastTree", "FastTreeMP", "VeryFastTree", and "IQTREE". These run with
         default settings only (and IQTREE includes "-m MFP" and "-B 1000"). To run
         any with more specific options you can use the output alignment file from
         GToTree (and the partitions file if wanted for mixed-model specification) as
         input into a dedicated treeing program (the GToTree `-N` option will generate
         the alignment only and skip internal treeing if wanted).
         """,
         FULL),

        ("-N | --no-tree",
         "do not make a tree",
         """
         Add this flag to stop after producing the concatenated alignment, skipping tree generation.
         """,
         FULL),

        ("-M | --muscle-threads <int>",
         "number of threads passed to muscle; default: 5",
         """
         The number of threads muscle will use during alignment. Keep in mind this
         will be multiplied by the number of jobs running concurrently based on
         the `-j` parameter.
         """,
         FULL),

        ("-X | --no-super5",
         "override super5 alignment",
         """
         If working with greater than 1,000 target genomes, GToTree will by default
         use the 'super5' muscle alignment algorithm to increase the speed of the
         alignments. Provide this flag if you want to prevent that from happening and
         use the standard muscle alignment instead.

         Note: See sections in the link below for more information on working with
         many genomes:
             github.com/AstrobioMike/GToTree/wiki/things-to-consider
         """,
         FULL),

        ("-z | --nucleotide-mode",
         "nucleotide mode",
         """
         Provide this flag to make the alignment and/or tree with nucleotide sequences
         instead of amino-acid sequences. (GToTree still finds target genes based on
         amino-acid HMM searches.)

         Note: This mode cannot accept amino-acid (`-A`) or genbank files (`-g`) as input sources.
         """,
         FULL),

        ("-k | --keep-gene-alignments",
         "keep individual target-gene alignments",
         """
         Add this flag to keep individual-gene alignment files.
         """,
         FULL),

        ("--tmp-dir <path>",
         "temporary directory location; default: <output-dir>/gtt-tmp-*",
         """
         If you want to specify where the temporary working directory will be created,
         you can provide the path to this parameter.
         """,
         FULL),

        ("-d | --debug",
         "debug mode",
         """
         Provide this flag if you'd like to keep the temporary directory.
         """,
         FULL),
    ]),
]


################################################################################
# the two renderings
################################################################################

REQUIRED_INPUTS = [
    ("-w  <str>",
     ("wanted ref tax to include (see full help for details)",
      "wanted ref tax to include (see more details below)")),
    ("-a <file>", "single-column file of NCBI assembly accessions"),
    ("-f <file>", "single-column file with the paths to each fasta file"),
    ("-g <file>", "single-column file with the paths to each genbank file"),
    ("-A <file>", "single-column file with the paths to each amino-acid file"),
]

HMM_INPUT = (
    "-H <file>",
    ("target single-copy gene HMMs to use (can be a path or the name of one of "
    "the pre-packaged sets; run 'gtt hmms' to view available gene-sets; will "
    "be auto-selected if not specified AND the `-w` parameter was used)"),
)

HMM_LEAD = "      2)  "
HMM_CONT = 22


def required_inputs_section(detailed=False):
    out = banner("REQUIRED INPUTS")
    out += "\n      1) Input genomes in one or any combination of the following formats:\n"

    for flags, brief in REQUIRED_INPUTS:
        out += entry(pick(flags, detailed), brief, detailed)

    hmm_flags, hmm_brief = HMM_INPUT
    hmm_flags = pick(hmm_flags, detailed)
    lead_len = len(HMM_LEAD) + 1 + len(hmm_flags) + 2
    wrapped = textwrap.fill(
        " ".join(pick(hmm_brief, detailed).split()),
        width=WIDTH,
        initial_indent=" " * lead_len,
        subsequent_indent=" " * HMM_CONT,
    )
    out += f"\n{HMM_LEAD}[{flag(hmm_flags)}] " + wrapped[lead_len:] + "\n"

    return out


def build_helpmenu(detailed=False):
    out = gtotree_header()

    out += banner("HELP INFO")
    out += """
  This program takes input genomes from various sources and ultimately produces
  a phylogenomic tree. You can find detailed usage information at:
                                  github.com/AstrobioMike/GToTree/wiki
"""
    if not detailed:
        out += condensed_note()

    out += required_inputs_section(detailed)

    out += banner("OPTIONAL SETTINGS")

    first_section = True
    for title, items in ENTRIES:
        shown = [i for i in items if detailed or i[3] is BRIEF]
        if not shown:
            continue

        out += subsection(title, extra_gap=detailed and not first_section)
        first_section = False
        for flags, brief, detail, _ in shown:
            out += ("\n" if detailed else "") + entry(flags, brief, detailed)
            if detailed and detail:
                out += detail_block(detail)

    gap = "\n" if detailed else ""
    out += subsection("Help and version:", extra_gap=detailed)
    out += gap + entry("-h | --help", "show the condensed help menu and exit", detailed)
    out += gap + entry("-S | --show-detailed-help",
                       "show the detailed help menu and exit", detailed)
    out += gap + entry("-v | --version", "show the GToTree version and exit", detailed)

    out += banner("EXAMPLE USAGE")
    out += "\n\tgtotree -w Alteromonas -f fasta-files.txt -H Gammaproteobacteria -D\n"

    if not detailed:
        out += condensed_note()

    out += "\n"

    return out


def print_helpmenu(detailed=False):
    """
    Render and write the menu. Built on demand rather than cached at import time --
    color_text() checks sys.stdout.isatty() when it runs, so a module-level constant
    would freeze whatever the tty state happened to be at import.
    """
    sys.stdout.write(build_helpmenu(detailed=detailed))


if __name__ == "__main__":
    print_helpmenu(detailed=any(a in ("-S", "--show-detailed-help")
                                for a in sys.argv[1:]))
