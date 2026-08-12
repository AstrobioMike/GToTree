"""Unit tests for gtotree/main_stages/filtering_genomes.py."""

import pytest  # type: ignore


def test_uncounted_genome_is_filtered_out_rather_than_crashing():
    """
    A genome whose hits were never counted (num_SCG_hits_after_filtering is None) must
    compare as zero. Comparing None against the threshold directly raises TypeError.
    """
    from gtotree.utils.misc.general import required_count

    threshold = required_count(4, 0.5)
    uncounted = None
    assert (uncounted or 0) < threshold
    with pytest.raises(TypeError):
        _ = uncounted < threshold
