import math


# residue monoisotopic-ish average masses, water removed
RESIDUE_MW = {
    "G":  57.05, "A":  71.08, "S":  87.08, "P":  97.12, "V":  99.13,
    "T": 101.10, "C": 103.14, "L": 113.16, "I": 113.16, "N": 114.10,
    "D": 115.09, "Q": 128.13, "K": 128.17, "E": 129.12, "M": 131.20,
    "H": 137.14, "F": 147.18, "R": 156.19, "Y": 163.18, "W": 186.21,
}
WATER_MW = 18.02

HYDROPHOBIC = set("AVLIFMWCYP")

# pKa values (Lehninger-style)
PKA_SIDECHAIN = {
    "D": 3.65, "E": 4.25, "H": 6.0,
    "C": 8.18, "Y": 10.07, "K": 10.53, "R": 12.48,
}
PKA_NTERM = 8.0    # alpha-NH3+
PKA_CTERM = 3.1    # alpha-COOH
ACIDIC = set("DECY")   # deprotonated form is negative
BASIC = set("HKR")     # protonated form is positive


def peptide_MW(seq):
    # seq: one-letter codes -> MW [g/mol]
    return sum(RESIDUE_MW[a] for a in seq) + WATER_MW


def peptide_charge(seq, pH):
    # one-letter seq, pH -> net charge
    q = 0.0
    # N-term basic
    q += 1.0 / (1.0 + 10**(pH - PKA_NTERM))
    # C-term acidic
    q -= 1.0 / (1.0 + 10**(PKA_CTERM - pH))
    for a in seq:
        if a in PKA_SIDECHAIN:
            pKa = PKA_SIDECHAIN[a]
            if a in ACIDIC:
                q -= 1.0 / (1.0 + 10**(pKa - pH))
            else:
                q += 1.0 / (1.0 + 10**(pH - pKa))
    return q


def peptide_mz(seq, z):
    # z: charge state (positive int) -> m/z
    return (peptide_MW(seq) + z * 1.00728) / z


def hp_pattern(seq):
    # seq -> "HPHHP..." string
    return "".join("H" if a in HYDROPHOBIC else "P" for a in seq)


def hp_matches(seq, pattern):
    # check if seq matches given HP pattern
    return hp_pattern(seq) == pattern


if __name__ == "__main__":
    s = "ACDEFGHIKLMNPQRSTVWY"
    print(f"MW({s}) = {peptide_MW(s):.2f} g/mol")
    for pH in (2, 7, 12):
        print(f"net charge at pH {pH}: {peptide_charge(s, pH):+.2f}")
    print(f"m/z (z=2) = {peptide_mz(s, 2):.2f}")
    print(f"HP pattern: {hp_pattern(s)}")
