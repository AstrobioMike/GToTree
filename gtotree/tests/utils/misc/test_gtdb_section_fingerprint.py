"""
`--gtdb-section` changes which genomes a run selects, so it has to be fingerprinted.

The subtlety is that it must be fingerprinted RESOLVED. The flag is declared with
`default=None` so that "typed it" stays distinguishable from "got the default", and
build_fingerprint's field loop copies args verbatim -- so the naive version stores
None for `GToTree -w Bacteria` and "reps" for `GToTree -w Bacteria --gtdb-section
reps`. Those two invocations select exactly the same genomes, and a resume across
them would be refused for a difference that doesn't exist.

Fingerprinting the effective value fixes that without collapsing the distinction the
messaging needs, because the messaging reads args, not the fingerprint.
"""

import pytest  # type: ignore

from gtotree.utils.misc.preflight_checks import build_fingerprint


class _Args:
    def __init__(self, **kw):
        self.wanted_ref_tax = ["Bacteria"]
        self.ncbi_accessions = None
        self.genbank_files = None
        self.fasta_files = None
        self.amino_acid_files = None
        self.mapping_file = None
        self.target_pfams_file = None
        self.target_kos_file = None
        self.exclusion_list = None
        self.hmm = "Universal"
        self.target_rank = None
        self.target_domain = None
        self.derep_rank = "auto"
        self.source = "gtdb"
        self.ncbi_section = "both"
        self.gtdb_section = None
        self.add_gtdb_tax = False
        self.add_ncbi_tax = False
        self.lineage = ""
        self.seq_length_cutoff = 0.2
        self.gene_representation_cutoff = 0.1
        self.genome_hits_cutoff = 0.5
        self.best_hit_mode = False
        self.no_super5 = False
        self.no_tree = False
        self.tree_program = "FastTreeMP"
        self.nucleotide_mode = False
        self.keep_gene_alignments = False
        self.__dict__.update(kw)


class TestGtdbSectionFingerprint:

    def test_the_section_is_fingerprinted_at_all(self):
        assert "gtdb_section" in build_fingerprint(_Args())

    def test_changing_the_section_invalidates_a_resume(self):
        # it changes the candidate pool, so it must change what a resume will accept
        reps = build_fingerprint(_Args(gtdb_section="reps"))
        every = build_fingerprint(_Args(gtdb_section="all"))

        assert reps["gtdb_section"] != every["gtdb_section"]

    def test_the_default_and_an_explicit_reps_fingerprint_identically(self):
        # same selection, so it must resume across
        defaulted = build_fingerprint(_Args(gtdb_section=None))
        explicit = build_fingerprint(_Args(gtdb_section="reps"))

        assert defaulted["gtdb_section"] == explicit["gtdb_section"] == "reps"

    def test_an_absent_attribute_does_not_blow_up(self):
        # bit mirrors these fingerprints, and its args objects don't always carry
        # every GToTree flag
        args = _Args()
        del args.gtdb_section

        assert build_fingerprint(args)["gtdb_section"] == "reps"
