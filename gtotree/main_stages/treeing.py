from gtotree.utils.messaging import report_processing_stage
from gtotree.utils.tree_program_handling import run_tree_building

def make_tree(args, run_data):

    if args.no_tree:
        return run_data

    report_processing_stage("treeing")

    if run_data.final_tree_path is not None:
        print("")
        print(f"The tree was built with {args.tree_program}.".center(82))

    else:
        run_data = run_tree_building(args, run_data)

    return run_data
