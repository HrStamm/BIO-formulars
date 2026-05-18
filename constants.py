"""Exam reference constants from provided sets.

All values are as stated in the exam instructions. Units in comments.
"""

# Fundamental physical constants
BOLTZMANN_K = 1.38065e-23  # J/K
KBT_ROOM_KJ_PER_MOL = 2.5  # kJ/mol
KBT_ROOM_J_PER_MOL = KBT_ROOM_KJ_PER_MOL * 1e3  # J/mol

# Material and fluid properties
VISCOSITY_WATER_25C_PA_S = 9e-4  # Pa*s
DENSITY_PROTEIN_G_PER_CM3 = 1.4  # g/cm^3
DENSITY_PROTEIN_KG_PER_M3 = DENSITY_PROTEIN_G_PER_CM3 * 1e3  # kg/m^3
DENSITY_GOLD_G_PER_CM3 = 19.3  # g/cm^3
DENSITY_GOLD_KG_PER_M3 = DENSITY_GOLD_G_PER_CM3 * 1e3  # kg/m^3

OSMOTIC_PRESSURE_GLUCOSE_1M_BAR = 25.0  # bar
OSMOTIC_PRESSURE_GLUCOSE_1M_PA = OSMOTIC_PRESSURE_GLUCOSE_1M_BAR * 1e5  # Pa
OSMOTIC_PRESSURE_SUCROSE_100MM_BAR = 2.5  # bar
OSMOTIC_PRESSURE_SUCROSE_100MM_PA = OSMOTIC_PRESSURE_SUCROSE_100MM_BAR * 1e5  # Pa

# Biological and biochemical parameters
ATP_HYDROLYSIS_KJ_PER_MOL = 50.0  # kJ/mol
ATP_HYDROLYSIS_J_PER_MOL = ATP_HYDROLYSIS_KJ_PER_MOL * 1e3  # J/mol

ECOLI_LINEAR_DIM_M = 1e-6  # m (order of magnitude)
ECOLI_DIAMETER_M = 2e-6  # m

D_BACTERIAL_RIBOSOME_M2_S = 1e-11  # m^2/s
D_TYPICAL_PROTEIN_M2_S = 1e-10  # m^2/s
TYPICAL_PROTEIN_RADIUS_M = 2e-9  # m (stated radius)
TYPICAL_PROTEIN_DIAMETER_M = 5e-9  # m (stated diameter)

INVERTASE_MW_KDA = 100.0  # kDa
INVERTASE_MW_G_PER_MOL = INVERTASE_MW_KDA * 1e3  # g/mol
INVERTASE_KCAT_PER_S = 1e4  # 1/s

# Amino acid molecular masses (free amino acids, g/mol).
# Values are residue masses + water (18.02), rounded to 2 decimals.
AMINO_ACID_MW_1LETTER_G_PER_MOL = {
    "G": 75.07,  # Glycine
    "A": 89.10,  # Alanine
    "S": 105.10,  # Serine
    "P": 115.14,  # Proline
    "V": 117.15,  # Valine
    "T": 119.12,  # Threonine
    "C": 121.16,  # Cysteine
    "L": 131.18,  # Leucine
    "I": 131.18,  # Isoleucine
    "N": 132.12,  # Asparagine
    "D": 133.11,  # Aspartate
    "Q": 146.15,  # Glutamine
    "K": 146.19,  # Lysine
    "E": 147.14,  # Glutamate
    "M": 149.22,  # Methionine
    "H": 155.16,  # Histidine
    "F": 165.20,  # Phenylalanine
    "R": 174.21,  # Arginine
    "Y": 181.20,  # Tyrosine
    "W": 204.23,  # Tryptophan
}

AMINO_ACID_MW_3LETTER_G_PER_MOL = {
    "Gly": 75.07,
    "Ala": 89.10,
    "Ser": 105.10,
    "Pro": 115.14,
    "Val": 117.15,
    "Thr": 119.12,
    "Cys": 121.16,
    "Leu": 131.18,
    "Ile": 131.18,
    "Asn": 132.12,
    "Asp": 133.11,
    "Gln": 146.15,
    "Lys": 146.19,
    "Glu": 147.14,
    "Met": 149.22,
    "His": 155.16,
    "Phe": 165.20,
    "Arg": 174.21,
    "Tyr": 181.20,
    "Trp": 204.23,
}
