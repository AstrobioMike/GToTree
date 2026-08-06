"""
`gtt search-kos` -- search a set of input genomes for a list of target KOs.

Everything of substance lives in `gtotree.utils.target_search`, driven by the KO
`TargetSearchSpec`. This module exists so the `gtt` dispatcher can find a
`build_parser`/`main` pair at a predictable path next to the rest of the KO code.
"""

from gtotree.utils.target_search import target_search_cli
from gtotree.utils.target_search.target_search_spec import get_spec


def _spec():
    return get_spec("ko")


def build_parser(parent_subparsers=None):
    return target_search_cli.build_parser(_spec(), parent_subparsers=parent_subparsers)


def main():  # pragma: no cover
    target_search_cli.main(_spec())


if __name__ == "__main__":
    main()
