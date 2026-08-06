"""
Shared implementation behind `gtt search-pfams` and `gtt search-kos`.

Both subcommands are the same program with a different target type: take the same
genome inputs the main GToTree driver accepts, preprocess them, search each one for a
set of user-supplied targets, and write out the same per-genome results, hit-count
matrix, and per-target hit sequences a full GToTree run would have produced with `-p`
or `-K`.

Everything that differs between the two is captured in a `TargetSearchSpec`
(`target_search_spec.py`); everything else is shared. The thin per-subcommand entry
points live next to their respective target code, at
`gtotree/utils/pfam/search_pfams_cli.py` and `gtotree/utils/ko/search_kos_cli.py`,
so the `gtt` dispatcher can find them the same way it finds every other subcommand.
"""
