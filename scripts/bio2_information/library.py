import math


def min_codon_length(n_amino_acids, n_bases=4):
    # n_aa: number of amino acids to encode, n_bases: alphabet size -> int
    L = math.ceil(math.log(n_amino_acids) / math.log(n_bases))
    while n_bases**L < n_amino_acids:
        L += 1
    return L


def library_size(n_positions, alphabet=20):
    # n_positions: randomized residues, alphabet: 20 (natural) or 4 (NT)
    return alphabet ** n_positions


def dna_sequence_space(length):
    # length [nt] -> 4^length
    return 4 ** length


def compare(values):
    # values: list of (label, expression_result) -> sorted desc
    return sorted(values, key=lambda r: r[1], reverse=True)


def bits_in_genome(genome_bp):
    # bp -> bits, bytes
    bits = 2 * genome_bp
    return bits, bits / 8


if __name__ == "__main__":
    print(f"min codon length for 25 aa (4 bases) = {min_codon_length(25)}")
    print(f"min codon length for 25 aa (3 bases) = {min_codon_length(25, n_bases=3)}")
    print(f"library size, 6 randomized aa positions = {library_size(6):.3e}")
    print(f"DNA space, 30 nt = {dna_sequence_space(30):.3e}")

    items = [
        ("4^10",       4**10),
        ("20^4",       20**4),
        ("20^6",       20**6),
        ("4^20",       4**20),
    ]
    for label, v in compare(items):
        print(f"{label:8s} = {v:.3e}")

    bits, bytes_ = bits_in_genome(3e9)
    print(f"3 Gb human genome: {bits:.2e} bits, {bytes_/1e9:.2f} GB")
