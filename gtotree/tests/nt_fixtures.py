"""
Synthetic nucleotide genomes and a matching SCG profile, generated at test time

These exist to exercise the gene-calling (prodigal) and nucleotide-mode (`-z`) paths,
which need real nucleotide input. Generated rather than committed for storage reasons
"""

import random

import pyhmmer.easel as easel  # type: ignore
import pyhmmer.plan7 as plan7  # type: ignore

# One codon per amino acid, deliberately biased: a consistent codon usage on the forward
# strand gives prodigal a clear coding signal to train on, so it calls the genes we
# planted rather than spurious reverse-strand ORFs.
CODON = {"A": "GCG", "R": "CGC", "N": "AAC", "D": "GAC", "C": "TGC", "Q": "CAG",
         "E": "GAG", "G": "GGC", "H": "CAC", "I": "ATC", "L": "CTG", "K": "AAG",
         "M": "ATG", "F": "TTC", "P": "CCG", "S": "AGC", "T": "ACC", "W": "TGG",
         "Y": "TAC", "V": "GTG"}
AMINO_ACIDS = "ARNDCQEGHILKMFPSTWYV"

SCG_LENGTH = 40          # comfortably above prodigal's 30-codon minimum gene length
MIN_GENOME_BP = 20_100   # prodigal will not train below 20,000
N_GENOMES = 4
GATHERING_CUTOFF = 20.0

SCG_NAME = "NT_MOCK"
SCG_ACCESSION = "Mike.2"


def _orf(protein):
    return "ATG" + "".join(CODON[a] for a in protein) + "TGA"


def _intergenic(rng, n):
    # AT-rich, bracketed by stops in all three reverse frames, to suppress spurious
    # reverse-strand ORFs that would otherwise compete during prodigal's training
    return "TTATTAGTCA" + "".join(rng.choice("ATATGC") for _ in range(n)) + "TTAATTAGTA"


def _mutate(protein, rng, n):
    chars = list(protein)
    for _ in range(n):
        chars[rng.randrange(len(chars))] = rng.choice(AMINO_ACIDS)
    return "".join(chars)


def scg_variants():
    """The conserved single-copy gene, one slightly different copy per genome."""
    base = random.Random(99)
    reference = "".join(base.choice(AMINO_ACIDS) for _ in range(SCG_LENGTH))

    variants = {}
    for n in range(1, N_GENOMES + 1):
        rng = random.Random(1000 + n)
        variants[f"nt-mock-{n}"] = reference if n == 1 else _mutate(reference, rng, 3)
    return reference, variants


def _genome_sequence(rng, scg_protein):
    """One planted SCG among filler ORFs, padded past prodigal's training minimum."""
    parts = []
    planted = False
    while sum(len(p) for p in parts) < MIN_GENOME_BP:
        filler = "".join(rng.choice(AMINO_ACIDS) for _ in range(rng.randint(60, 140)))
        parts.append(_orf(filler))
        parts.append(_intergenic(rng, rng.randint(20, 45)))
        if not planted and sum(len(p) for p in parts) > MIN_GENOME_BP // 2:
            parts.append(_orf(scg_protein))
            parts.append(_intergenic(rng, 30))
            planted = True
    return "".join(parts)


def write_genomes(directory):
    """Write the genome FASTAs into `directory`; returns their paths in order."""
    _, variants = scg_variants()
    paths = []
    for n, (name, scg) in enumerate(variants.items(), start=1):
        rng = random.Random(1000 + n)
        _mutate(scg, rng, 0)                      # keep rng in step with scg_variants
        sequence = _genome_sequence(rng, scg)
        path = directory / f"{name}.fna"
        path.write_text(
            f">{name}_contig1\n"
            + "\n".join(sequence[i:i + 70] for i in range(0, len(sequence), 70))
            + "\n")
        paths.append(path)
    return paths


def write_hmm(directory):
    """
    Build the SCG profile in-process and write it out; returns the path.

    Built from the reference protein with pyhmmer rather than shelling out to
    muscle + hmmbuild, so no external tools are needed to produce the fixture.
    """
    reference, _ = scg_variants()

    alphabet = easel.Alphabet.amino()
    sequence = easel.TextSequence(
        name=SCG_NAME.encode(),
        accession=SCG_ACCESSION.encode(),
        sequence=reference).digitize(alphabet)

    builder = plan7.Builder(alphabet, seed=42)
    hmm, _, _ = builder.build(sequence, plan7.Background(alphabet))
    # GToTree searches with gathering thresholds, so the profile needs one
    hmm.cutoffs.gathering = (GATHERING_CUTOFF, GATHERING_CUTOFF)

    path = directory / "nt-mock.hmm"
    with open(path, "wb") as fh:
        hmm.write(fh)
    return path
