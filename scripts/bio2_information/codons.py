import math


DEGENERACY = {
    "A": 4, "R": 6, "N": 2, "D": 2, "C": 2,
    "Q": 2, "E": 2, "G": 4, "H": 2, "I": 3,
    "L": 6, "K": 2, "M": 1, "F": 2, "P": 4,
    "S": 6, "T": 4, "W": 1, "Y": 2, "V": 4,
    "*": 3,
}

CODON_TABLE = {
    "UUU":"F","UUC":"F","UUA":"L","UUG":"L",
    "UCU":"S","UCC":"S","UCA":"S","UCG":"S",
    "UAU":"Y","UAC":"Y","UAA":"*","UAG":"*",
    "UGU":"C","UGC":"C","UGA":"*","UGG":"W",
    "CUU":"L","CUC":"L","CUA":"L","CUG":"L",
    "CCU":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAU":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CGU":"R","CGC":"R","CGA":"R","CGG":"R",
    "AUU":"I","AUC":"I","AUA":"I","AUG":"M",
    "ACU":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAU":"N","AAC":"N","AAA":"K","AAG":"K",
    "AGU":"S","AGC":"S","AGA":"R","AGG":"R",
    "GUU":"V","GUC":"V","GUA":"V","GUG":"V",
    "GCU":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAU":"D","GAC":"D","GAA":"E","GAG":"E",
    "GGU":"G","GGC":"G","GGA":"G","GGG":"G",
}


def num_sequences_encoding(protein_seq, include_stop=True):
    # protein_seq: one-letter codes, e.g. "MAW" -> int
    n = 1
    for aa in protein_seq:
        n *= DEGENERACY[aa]
    if include_stop:
        n *= DEGENERACY["*"]
    return n


def translate(nt_seq):
    # nt_seq: DNA or RNA string -> protein string (stops at *)
    s = nt_seq.upper().replace("T", "U")
    out = []
    for i in range(0, len(s) - 2, 3):
        aa = CODON_TABLE.get(s[i:i+3])
        if aa is None:
            out.append("?")
            continue
        if aa == "*":
            break
        out.append(aa)
    return "".join(out)


def orf_probability(n_codons, p_start=1/64, p_stop=3/64):
    # n_codons: ORF length in codons (incl. start, excl. stop)
    return p_start * (1 - p_stop)**(n_codons - 1) * p_stop


def expected_orfs(genome_bp, n_codons, both_strands=True):
    # genome_bp [bp], n_codons: target length
    p = orf_probability(n_codons)
    return (2 if both_strands else 1) * genome_bp * p


if __name__ == "__main__":
    print(f"MAW: {num_sequences_encoding('MAW')} sequences (incl. stop)")
    print(f"MET-ALA-TRP coded by AUGGCAUGG... -> {translate('AUGGCAUGGUAA')}")
    print(f"P(ORF length 100 codons) = {orf_probability(100):.3e}")
    print(f"E[ORFs >= 100 codons] in 5e6 bp genome = {expected_orfs(5e6, 100):.1f}")
