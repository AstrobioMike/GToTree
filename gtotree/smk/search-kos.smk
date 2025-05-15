from gtotree.utils.general import (read_run_data,
                                   write_run_data)
from gtotree.main_stages.additional_ko_searching import (run_ko_search)

# snakefile, for each genome, needs to
    # do the KO search
    # parse the output table to get each unique KOs associated gene IDs
    # write those genes to a file under that KOs name (for that genome)
    # keep a count of hits for each KO

# either in the rule all, or back at the main-stage runner, combine all
# seqs from all genomes for each individual KO in one file (seq ids should
# already be assembly names at this point)
