from pathlib import Path

import numpy as np

from rico_dna_sequence_data_handling.main import (
    one_hot_encoding,
    ordinal_encoding,
    print_sequences_data
)


def test_print_sequences_data() -> None:
    fasta_file_path = Path("tests/data/pdb_seqres_small.txt")
    assert print_sequences_data(fasta_file_path)


def test_ordinal_encoding() -> None:
    ordinal_encoded_seq = ordinal_encoding("GAATTCTCGAA")
    np.testing.assert_array_equal(
        ordinal_encoded_seq,
        np.array([0.75, 0.25, 0.25, 1.0, 1.0, 0.5, 1.0, 0.5, 0.75, 0.25, 0.25])
    )


def test_one_hot_encoding() -> None:
    one_hot_encoded_seq = one_hot_encoding("GAATTCTCGAA")
    np.testing.assert_array_equal(
        one_hot_encoded_seq,
        np.array(
            [
                [0, 0, 1],
                [1, 0, 0],
                [1, 0, 0],
                [0, 0, 0],
                [0, 0, 0],
                [0, 1, 0],
                [0, 0, 0],
                [0, 1, 0],
                [0, 0, 1],
                [1, 0, 0],
                [1, 0, 0],
            ]
        ),
    )
