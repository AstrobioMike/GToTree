import os
from dataclasses import dataclass
import subprocess
from importlib.metadata import version


@dataclass
class CitationsInfo:
    gtotree = "Lee MD. GToTree: a user-friendly workflow for phylogenomics. Bioinformatics. 2019; (March):1-3. doi.org/10.1093/bioinformatics/btz188\n\n"
    snakemake = "Koster J and Rahmann S. Snakemake - a scalable bioinformatics workflow engine. Bioinformatics. 2012; 28, 2520-2522. doi.org/10.1093/bioinformatics/bts480\n\n"
    hmmer = "Eddy SR. Accelerated profile HMM searches. PLoS Comput. Biol. 2011; (7)10. doi.org/10.1371/journal.pcbi.1002195\n\n"
    muscle = "Edgar RC. MUSCLE v5 enables improved estimates of phylogenetic tree confidence by ensemble bootstrapping. bioRxiv. 2021.06.20.449169. doi.org/10.1101/2021.06.20.449169\n\n"
    trimal = "Gutierrez SC. et al. TrimAl: a Tool for automatic alignment trimming. Bioinformatics. 2009; 25, 1972-1973. doi.org/10.1093/bioinformatics/btp348\n\n"
    prodigal = "Hyatt, D. et al. Gene and translation initiation site prediction in metagenomic sequences. Bioinformatics. 2010; 28, 2223-2230. doi.org/10.1186/1471-2105-11-119\n\n"
    taxonkit = "Shen W and Ren H. TaxonKit: a practical and efficient NCBI Taxonomy toolkit. Journal of Genetics and Genomics. 2021. doi.org/10.1016/j.jgg.2021.03.006\n\n"
    gtdb = "Parks DH et al. A complete domain-to-species taxonomy for Bacteria and Archaea. Nat. Biotech. 2020. doi.org/10.1038/s41587-020-0501-8\n\n"
    fasttree = "Price MN et al. FastTree2 - approximately maximum-likelihood trees for large alignments. PLoS One. 2010; 5. doi.org/10.1371/journal.pone.0009490\n\n"
    veryfasttree = "Pineiro C et al. VeryFastTree: speeding up the estimation of phylogenies for large alignments through parallelization and vectorization strategies. Bioinformatics. 2020. doi.org/10.1093/bioinformatics/btaa582\n\n"
    iqtree = (
        "Nguyen L-T et al. IQ-TREE: a fast and effective stochastic algorithm for estimating maximum likelihood phylogenies. Mol. Biol. Evol. 2015; 32, 268-274. doi.org/10.1093/molbev/msu300\n"
        "Kalyaanamoorthy et al. ModelFinder: fast model selection for accurate phylogenetic estimates. Nat. Methods 2017; 14:587-589. doi.org/10.1038/nmeth.4285\n"
        "Hoang et al. UFBoot2: Improving the ultrafast bootstrap approximation. Mol. Biol. Evol. 2018; 35:518-522. doi.org/10.1093/molbev/msx281\n\n"
    )
    pfam = "Mistry J et al. Pfam: the protein families database in 2021. Nucleic Acids Research. 2021. doi.org/10.1093/nar/gkaa913\n\n"
    kofamscan = "Aramaki, T et al. KofamKOALA: KEGG Ortholog assignment based on profile HMM and adaptive score threshold. Bioinformatics. 2020. doi.org/10.1093/bioinformatics/btz859\n\n"
    universal_SCG_set = "Hug LA et al. A new view of the tree of life. Nat. Microbiol. 2016; 1, 1-6. doi.org/10.1038/NMICROBIOL.2016.48\n\n"


def generate_citations_info(run_data):

    citations_file = os.path.join(run_data.output_dir, "citations.txt")
    citations_info = CitationsInfo

    with open(citations_file, "w") as outfile:

        outfile.write(f"GToTree v{version('GToTree')}\n")
        outfile.write(citations_info.gtotree)

        outfile.write(f"Snakemake v{get_snakemake_version()}\n")
        outfile.write(citations_info.snakemake)

        outfile.write(f"HMMER v{get_hmmer_version()}\n")
        outfile.write(citations_info.hmmer)

        outfile.write(f"MUSCLE v{get_muscle_version()}\n")
        outfile.write(citations_info.muscle)

        outfile.write(f"TrimAl {get_trimal_version()}\n")
        outfile.write(citations_info.trimal)

        if run_data.tools_used.prodigal_used:
            outfile.write(f"Prodigal {get_prodigal_version()}\n")
            outfile.write(citations_info.prodigal)

        if run_data.tools_used.taxonkit_used:
            outfile.write(f"TaxonKit v{get_taxonkit_version()}\n")
            outfile.write(citations_info.taxonkit)

        if run_data.tools_used.gtdb_used:
            outfile.write(f"Genome Taxonomy Database (GTDB)\n")
            outfile.write(citations_info.gtdb)

        if run_data.tools_used.fasttree_used:
            outfile.write(f"FastTree v{get_fasttree_version()}\n")
            outfile.write(citations_info.fasttree)

        if run_data.tools_used.veryfasttree_used:
            outfile.write(f"VeryFastTree v{get_veryfasttree_version()}\n")
            outfile.write(citations_info.veryfasttree)

        if run_data.tools_used.iqtree_used:
            outfile.write(f"IQ-TREE v{get_iqtree_version()}\n")
            outfile.write(citations_info.iqtree)

        if run_data.tools_used.kofamscan_used:
            outfile.write(f"KoFamScan v{get_kofamscan_version()}\n")
            outfile.write(citations_info.kofamscan)

        if run_data.tools_used.pfam_db_used:
            outfile.write(f"Pfam Database\n")
            outfile.write(citations_info.pfam)

        if run_data.tools_used.universal_SCGs_used:
            outfile.write(f"Universal SCG-set\n")
            outfile.write(citations_info.universal_SCG_set)


def get_snakemake_version():
    snakemake_version = subprocess.run('snakemake -v',
                                       shell = True, capture_output = True, text = True)
    return snakemake_version.stdout.strip()

def get_hmmer_version():
    hmm_version = subprocess.run('hmmsearch -h | head -n 2 | tail -n 1 | tr -s " " "\t" | cut -f 3',
                                 shell = True, capture_output = True, text = True)
    return hmm_version.stdout.strip()

def get_muscle_version():
    muscle_version = subprocess.run('muscle -version | tr -s " " "\t" | cut -f 2 | head -n 1',
                                    shell = True, capture_output = True, text = True)
    return muscle_version.stdout.strip()

def get_trimal_version():
    trimal_version = subprocess.run('trimal --version | grep "trim" | tr -s " " "\t" | cut -f 2',
                                    shell = True, capture_output = True, text = True)
    return trimal_version.stdout.strip()

def get_prodigal_version():
    prodigal_version = subprocess.run('prodigal -v 2>&1 | grep Prodigal | tr -s " " "\t" | cut -f 2 | tr -d ":" | sed "s/V/v/"',
                                      shell = True, capture_output = True, text = True)
    return prodigal_version.stdout.strip()

def get_taxonkit_version():
    taxonkit_version = subprocess.run('taxonkit -h | grep Version | tr -s " " "\t" | cut -f 2',
                                      shell = True, capture_output = True, text = True)
    return taxonkit_version.stdout.strip()

def get_fasttree_version():
    fasttree_version = subprocess.run('FastTree -expert 2>&1 | head -n 1 | tr -s " " "\t" | cut -f 5',
                                      shell = True, capture_output = True, text = True)
    return fasttree_version.stdout.strip()

def get_veryfasttree_version():
    veryfasttree_version = subprocess.run('VeryFastTree -h | head -n 1 | cut -f 2 -d " "',
                                          shell = True, capture_output = True, text = True)
    return veryfasttree_version.stdout.strip()

def get_iqtree_version():
    iqtree_version = subprocess.run('iqtree -V | head -n 1 | tr -s " " "\t" | cut -f 4',
                                    shell = True, capture_output = True, text = True)
    return iqtree_version.stdout.strip()

def get_kofamscan_version():
    kofamscan_version = subprocess.run('exec_annotation -v | cut -f 2 -d " "',
                                       shell = True, capture_output = True, text = True)
    return kofamscan_version.stdout.strip()
