import logging
import re
from pathlib import Path

import numpy as np
import numpy.typing as npt
from Bio import SeqIO
from sklearn.preprocessing import LabelEncoder, OneHotEncoder

logger = logging.getLogger(__name__)


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


def ordinal_encoding(sequence_string: str) -> npt.NDArray[np.complex64]:
    """
    We encode each nitrogen base as an ofdinal value.
    E.g. "ATGC" == [0.25, 0.5, 0.75, 1.0], any other
    base such as "N" can be a 0.
    """
    # cast to nump.array
    sequence_str: str = sequence_string.lower()
    sequence_str = re.sub("[^acgt]", "n", sequence_str)
    sequence_npa: npt.NDArray[np.complex64] = np.array(list(sequence_str))

    # create label encoder with 'acgtn' alphabet
    label_encoder = LabelEncoder()
    label_encoder.fit(np.array(["a", "c", "g", "t", "z"]))

    # to ordinal encoding
    integer_encoded = label_encoder.transform(sequence_npa)
    float_encoded = integer_encoded.astype(float)
    float_encoded[float_encoded == 0] = 0.25  # A
    float_encoded[float_encoded == 1] = 0.50  # C
    float_encoded[float_encoded == 2] = 0.75  # G
    float_encoded[float_encoded == 3] = 1.00  # T
    float_encoded[float_encoded == 4] = 0.00  # enything else

    return float_encoded


def one_hot_encoding(sequence_string: str) -> npt.NDArray[np.complex64]:
    # cast to nump.array
    sequence_str: str = sequence_string.lower()
    sequence_str = re.sub("[^acgt]", "n", sequence_str)
    sequence_npa: npt.NDArray[np.complex64] = np.array(list(sequence_str))

    # create label encoder with 'acgtn' alphabet
    label_encoder = LabelEncoder()
    label_encoder.fit(np.array(["a", "c", "g", "t", "z"]))
    integer_encoded = label_encoder.transform(sequence_npa)
    integer_encoded = integer_encoded.reshape(len(integer_encoded), 1)
    onehot_encoded = OneHotEncoder(sparse=False, dtype=int)
    onehot_encoded = onehot_encoded.fit_transform(integer_encoded)
    onehot_encoded = np.delete(onehot_encoded, -1, 1)

    return onehot_encoded


def show_encodings(fasta_file_path: Path) -> None:
    try:
        all_sequences = SeqIO.parse(fasta_file_path, "fasta")
        for sequence in all_sequences:
            logger.info(" Ordinal encoding:\n")
            print(ordinal_encoding(str(sequence.seq)))
            logger.info(" One-hot encoding:\n")
            one_hot_seq = one_hot_encoding(str(sequence.seq))
            print(one_hot_seq)
    except AttributeError as e:
        print(f"An error happened: {e}")


def main(fasta_file_path: Path) -> None:
    logging.basicConfig(level=logging.INFO)
    logger.info("Printing sequences:")
    print_sequences_data(fasta_file_path)
    show_encodings(fasta_file_path)


main(Path("tests/data/norway_rat_HBB.fasta"))
