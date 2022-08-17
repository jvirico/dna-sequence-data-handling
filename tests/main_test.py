from pathlib import Path

from rico_dna_sequence_data_handling.main import print_sequences_data


def test_print_sequences_data() -> None:
    fasta_file_path = Path("tests/data/pdb_seqres_small.txt")
    assert print_sequences_data(fasta_file_path)
