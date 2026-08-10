"""
This holds the main gtotree workflow stages
"""


class GenomeRemovalStage:
    """Where an input genome left the run."""

    NCBI_LOOKUP     = "ncbi-lookup"        # accession not found at NCBI
    NCBI_DOWNLOAD   = "ncbi-download"      # found, but download/prep failed
    GENBANK_PREP    = "genbank-prep"
    FASTA_PREP      = "fasta-prep"
    AMINO_ACID_PREP = "amino-acid-prep"
    HMM_SEARCH      = "hmm-search"         # search failed, or seq extraction failed
    SCG_HIT_FILTER  = "scg-hit-filter"     # too few SCG hits (`-G`)


# in the order the run can reach them; `alive_through` depends on this being accurate
GENOME_REMOVAL_STAGE_ORDER = (
    GenomeRemovalStage.NCBI_LOOKUP,
    GenomeRemovalStage.NCBI_DOWNLOAD,
    GenomeRemovalStage.GENBANK_PREP,
    GenomeRemovalStage.FASTA_PREP,
    GenomeRemovalStage.AMINO_ACID_PREP,
    GenomeRemovalStage.HMM_SEARCH,
    GenomeRemovalStage.SCG_HIT_FILTER,
)

# the four input types are processed one after another, but they're all "processing"
# as far as the overall input-genome summary is concerned
PREPROCESSING_REMOVAL_STAGES = (
    GenomeRemovalStage.NCBI_LOOKUP,
    GenomeRemovalStage.NCBI_DOWNLOAD,
    GenomeRemovalStage.GENBANK_PREP,
    GenomeRemovalStage.FASTA_PREP,
    GenomeRemovalStage.AMINO_ACID_PREP,
)

# which source list each processing stage can remove from, used to keep the
# per-source reports from counting another source's failures
NCBI_REMOVAL_STAGES = (GenomeRemovalStage.NCBI_LOOKUP,
                       GenomeRemovalStage.NCBI_DOWNLOAD)


class SCGRemovalStage:
    """Where a target SCG-set left the run."""

    NO_HITS       = "no-hits"          # nothing hit it at all
    GENE_FILTER   = "gene-filter"      # `-c` / `-r` filtering
    GENOME_FILTER = "genome-filter"    # nothing left after genomes were dropped
    ALIGNMENT     = "alignment"        # muscle or trimal failed


SCG_REMOVAL_STAGE_ORDER = (
    SCGRemovalStage.NO_HITS,
    SCGRemovalStage.GENE_FILTER,
    SCGRemovalStage.GENOME_FILTER,
    SCGRemovalStage.ALIGNMENT,
)

# what the "N had no hits or were filtered out" report covers
SCG_GENE_FILTERING_STAGES = (SCGRemovalStage.NO_HITS, SCGRemovalStage.GENE_FILTER)


class PipelineStage:
    """
    The stages of a full GToTree run, in order
    """

    PROCESS_GENOMES      = "process-genomes"
    FILTER_GENES         = "filter-genes"
    FILTER_GENOMES       = "filter-genomes"
    ALIGN_SCG_SETS       = "align-SCG-sets"
    CONCATENATE_SCG_SETS = "concatenate-SCG-sets"
    UPDATE_HEADERS       = "update-headers"
    BUILD_TREE           = "build-tree"
    FINALIZE             = "finalize"


PIPELINE_STAGES = (
    PipelineStage.PROCESS_GENOMES,
    PipelineStage.FILTER_GENES,
    PipelineStage.FILTER_GENOMES,
    PipelineStage.ALIGN_SCG_SETS,
    PipelineStage.CONCATENATE_SCG_SETS,
    PipelineStage.UPDATE_HEADERS,
    PipelineStage.BUILD_TREE,
    PipelineStage.FINALIZE,
)


class UnknownStage(ValueError):
    """A stage name that isn't in one of the vocabularies above."""


def validate_stage(stage, order, what):
    """
    Raise unless `stage` is a known member of `order`.

    Library code raises; the CLI layer never sees this one, because a bad stage name is
    a programming error rather than anything a user can cause.
    """
    if stage not in order:
        raise UnknownStage(
            f"unknown {what} {stage!r}; expected one of {list(order)}")
    return stage


def stages_through(stage, order):
    """
    Every stage at or before `stage`, as a set.

    This is what "was this still in the run as of stage X?" reduces to: an input is
    alive through X if it wasn't removed at any stage in here. Removals from *later*
    stages, which a resume has already loaded, don't count yet.
    """
    validate_stage(stage, order, "stage")
    return frozenset(order[:order.index(stage) + 1])
