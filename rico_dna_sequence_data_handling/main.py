from pathlib import Path

from Bio import SeqIO


def print_sequences_data(fasta_file_path: Path) -> bool:

    try:
        all_sequencies = SeqIO.parse(fasta_file_path, "fasta")
        for sequence in all_sequencies:
            print(sequence.id)
            print(sequence.seq)
            print(f"{len(sequence)}\n")
        return True
    except Exception as e:
        print(f"An error happened: {e}")
        return False


def main(fasta_file_path: Path) -> None:
    print_sequences_data(fasta_file_path)
