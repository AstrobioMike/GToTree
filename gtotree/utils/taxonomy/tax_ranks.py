RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]

# fill value for a rank that is absent from a lineage
NA = "NA"

REFERENCE_VALUE = "reference genome"

LINEAGE_NAME_COLUMNS = list(RANKS)
LINEAGE_TAXID_COLUMNS = [f"{r}_taxid" for r in RANKS]


def accession_core(acc):
    """
    Canonical key for an assembly, shared by everything that has to match accessions
    ACROSS sources: GTDB's RS_/GB_ prefixes, the GCA_/GCF_ pairing, and version
    suffixes all normalise away.

        'RS_GCF_000005845.2' -> '000005845'
        'GCA_000005845.3'    -> '000005845'
        'na' / '' / junk     -> ''
    """
    s = str(acc).strip()
    if s.startswith(("RS_", "GB_")):
        s = s[3:]
    if not s.upper().startswith(("GCA_", "GCF_")):
        return ""
    return s.split(".")[0].split("_")[-1]


def rank_index(rank):
    """0 for domain ... 6 for species. Raises ValueError on an unknown rank."""
    r = str(rank).strip().lower()
    if r not in RANKS:
        raise ValueError(
            f"'{rank}' is not one of the supported ranks: {', '.join(RANKS)}")
    return RANKS.index(r)


def validate_derep_rank(wanted_rank, derep_rank):
    """
    The --derep-rank rank must be the SAME as, or FINER than, the rank
    of --wanted-ref-tax. Asking for one genome per family within a genus is
    incoherent (the genus sits below family).

    Returns None if OK, else an explanatory string.
    """
    w = rank_index(wanted_rank)
    d = rank_index(derep_rank)
    if d < w:
        return (f"--derep-rank '{derep_rank}' is a broader rank than "
                f"the target taxon's rank '{wanted_rank}'. The per-rank rank must "
                f"be the same or finer (e.g., target 'domain' + per-rank 'class').")
    return None
