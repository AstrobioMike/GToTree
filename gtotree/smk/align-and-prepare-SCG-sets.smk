import os
from gtotree.utils.general import (read_run_data,
                                   write_run_data)
from gtotree.utils.seqs import (run_muscle,
                                run_trimal,
                                add_needed_gap_seqs)

run_data = read_run_data(config['run_data_path'])
if run_data is None:
    raise ValueError("Run data not found")

SCG_dict = {SCG.id: SCG for SCG in run_data.get_all_SCG_targets_remaining()}
SCG_ids = list(SCG_dict.keys())

rule all:
    input:
        expand(f"{run_data.found_SCG_seqs_dir}/{{SCG}}.aligned", SCG=SCG_ids)
    run:
        for SCG in SCG_ids:
            SCG_target = SCG_dict[SCG]
            status_path = f"{run_data.found_SCG_seqs_dir}/{SCG}.aligned"

            with open(status_path, 'r') as f:
                for line in f:
                    ID, align_failed, trimal_failed = line.strip().split('\t')

                    if int(align_failed):
                        SCG_target.mark_removed()
                        SCG_target.aligned = False
                        SCG_target.reason_removed = "alignment failed"
                    elif int(trimal_failed):
                        SCG_target.mark_removed()
                        SCG_target.trimmed = False
                        SCG_target.reason_removed = "trimal failed"
                    else:
                        SCG_target.aligned = True
                        SCG_target.trimmed = True
                        SCG_target.ready_for_cat = True

        for f in input:
            os.remove(f)

        write_run_data(run_data)

rule align:
    output:
        f"{run_data.found_SCG_seqs_dir}/{{SCG}}.aligned"
    run:
        inpath = run_data.found_SCG_seqs_dir + f"/{wildcards.SCG}-genome-filtered.fasta"

        aligned_path = run_data.found_SCG_seqs_dir + f"/{wildcards.SCG}-aligned.fasta"
        align_log_path = run_data.found_SCG_seqs_dir + f"/{wildcards.SCG}-align.log"

        trimmed_path = run_data.found_SCG_seqs_dir + f"/{wildcards.SCG}-trimmed.fasta"
        trimmal_log_path = run_data.found_SCG_seqs_dir + f"/{wildcards.SCG}-trimmal.log"

        final_path = run_data.found_SCG_seqs_dir + f"/{wildcards.SCG}-final.fasta"

        align_failed = run_muscle(f"{wildcards.SCG}", run_data, inpath, aligned_path, align_log_path)
        trimal_failed = run_trimal(aligned_path, trimmed_path, trimmal_log_path)
        add_needed_gap_seqs(run_data, trimmed_path, final_path)

        with open(output[0], 'w') as f:
            f.write(f'{wildcards.SCG}\t{int(align_failed)}\t{int(trimal_failed)}\n')
