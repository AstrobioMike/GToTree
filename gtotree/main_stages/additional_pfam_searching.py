from gtotree.utils.pfam.pfam_handling import get_additional_pfam_targets
from gtotree.utils.messaging import (report_processing_stage,
                                     report_pfam_searching_update)


def search_pfams(run_data):

    report_processing_stage("additional-pfam-searching", run_data)

    if run_data.additional_pfam_searching_done:
        report_pfam_searching_update(run_data)
        return run_data

    run_data = get_additional_pfam_targets(run_data)

    report_pfam_searching_update(run_data)

    return run_data
