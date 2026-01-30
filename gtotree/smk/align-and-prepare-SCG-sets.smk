import os
import json
from gtotree.utils.general import (read_run_data,
                                   write_run_data)
from gtotree.utils.seqs import (run_muscle,
                                run_trimal,
                                add_needed_gap_seqs,
                                fasta_has_single_record)

run_data = read_run_data(config['run_data_path'])
if run_data is None:
    raise ValueError("Run data not found")

SCG_dict = {SCG.id: SCG for SCG in run_data.get_all_SCG_targets_remaining()}
SCG_ids = list(SCG_dict.keys())

rule all:
    input:
        expand(f"{run_data.found_SCG_seqs_dir}/{{SCG}}-align.json", SCG=SCG_ids)
    run:
        for SCG in SCG_ids:
            SCG_target = SCG_dict[SCG]
            status_path = f"{run_data.found_SCG_seqs_dir}/{SCG}-align.json"

            with open(status_path, 'r') as f:
                obj = json.load(f)

            ID = obj.get('scg', SCG)
            align_failed = bool(obj.get('align_failed', False))
            trimal_failed = bool(obj.get('trimal_failed', False))

            if align_failed:
                SCG_target.mark_removed("alignment failed")
                SCG_target.aligned = False
            elif trimal_failed:
                SCG_target.mark_removed("trimal failed")
                SCG_target.trimmed = False
            else:
                SCG_target.aligned = True
                SCG_target.trimmed = True
                SCG_target.ready_for_cat = True

        run_data.all_SCG_sets_aligned = True

        write_run_data(run_data)


rule align:
    output:
        f"{run_data.found_SCG_seqs_dir}/{{SCG}}-align.json"
    run:

        inpath = run_data.found_SCG_seqs_dir + f"/{wildcards.SCG}-genome-filtered{run_data.general_ext}"

        aligned_path = run_data.found_SCG_seqs_dir + f"/{wildcards.SCG}-aligned{run_data.general_ext}"
        align_log_path = run_data.found_SCG_seqs_dir + f"/{wildcards.SCG}-align.log"

        trimmed_path = run_data.found_SCG_seqs_dir + f"/{wildcards.SCG}-trimmed{run_data.general_ext}"
        trimmal_log_path = run_data.found_SCG_seqs_dir + f"/{wildcards.SCG}-trimmal.log"

        final_path = run_data.found_SCG_seqs_dir + f"/{wildcards.SCG}-final{run_data.general_ext}"

        if fasta_has_single_record(inpath):
            align_failed = False
            trimal_failed = False
            trimmed_source = inpath
            with open(align_log_path, 'w') as f:
                f.write("Single sequence in SCG set; skipping alignment.\n")
        else:
            align_failed = run_muscle(f"{wildcards.SCG}", run_data, inpath, aligned_path, align_log_path)
            trimal_failed = run_trimal(aligned_path, trimmed_path, trimmal_log_path)
            trimmed_source = trimmed_path

        add_needed_gap_seqs(run_data, trimmed_source, final_path)

        status_obj = {
            "scg": wildcards.SCG,
            "align_failed": align_failed,
            "trimal_failed": trimal_failed,
        }

        tmp_path = output[0] + ".tmp"
        with open(tmp_path, 'w') as fh:
            json.dump(status_obj, fh, indent=2, sort_keys=True)
        os.replace(tmp_path, output[0])
