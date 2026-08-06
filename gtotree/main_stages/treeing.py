from gtotree.utils.messaging import report_processing_stage
from gtotree.utils.tree_program_handling import run_tree_building
from gtotree.utils.general import file_is_usable_else_clear, write_run_data

def make_tree(args, run_data):

    if args.no_tree:
        return run_data

    report_processing_stage("treeing", run_data)

    if file_is_usable_else_clear(run_data.final_tree_path):
        print("")
        print(f"The tree was built with {args.tree_program}.".center(82))

    else:
        run_data = run_tree_building(args, run_data)
        write_run_data(run_data)


    return run_data
