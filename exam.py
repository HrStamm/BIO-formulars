"""Exam scratchpad. Edit numbers, uncomment the line you need, run:
    python3 exam.py
"""

import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from scripts.constants import *
from scripts.units import *

from scripts.bio1_physical.diffusion import (
    stokes_einstein_D, diffusion_distance, diffusion_time,
    D_from_trajectory, R_from_D, scale_D, boltzmann_entropy,
)
from scripts.bio1_physical.osmosis import (
    osmotic_pressure, effective_osmolarity, rank_osmotic, VANT_HOFF,
)
from scripts.bio1_physical.flow import (
    hagen_poiseuille_Q, hagen_poiseuille_dP, jurin_height, reynolds,
)
from scripts.bio1_physical.membrane import (
    sphere_area, sphere_volume, capacitance, capacitor_charge,
    capacitor_energy, n_charges, n_ions_in_cell,
)
from scripts.bio1_physical.energy import (
    atp_per_day, polymer_force, energy_density,
    fat_energy_density, atp_energy_density,
)

from scripts.bio2_information.codons import (
    DEGENERACY, CODON_TABLE, num_sequences_encoding,
    translate, orf_probability, expected_orfs,
)
from scripts.bio2_information.library import (
    min_codon_length, library_size, dna_sequence_space, compare, bits_in_genome,
)
from scripts.bio2_information.pcr import (
    pcr_copies, cycles_for_copies,
    transcript_concentration, molecules_from_concentration,
)

from scripts.bio3_proteins.kinetics import (
    rate_enhancement_kBT, rate_enhancement_kJmol,
    michaelis_menten, vmax_from_kcat, fractional_occupancy, hill_occupancy,
)
from scripts.bio3_proteins.peptide import (
    peptide_MW, peptide_charge, peptide_mz, hp_pattern, hp_matches,
    RESIDUE_MW, PKA_SIDECHAIN, HYDROPHOBIC,
)
from scripts.bio3_proteins.stoichiometry import (
    enzyme_mass_for_throughput, antibody_mass_for_toxin,
    catalyst_mass_for_production,
)

from scripts.bio4_evolution.virus import (
    virus_volume_fraction, virus_dry_mass_fraction,
    replication_cycles, endogenous_retroviruses,
)
from scripts.bio4_evolution.mutations import (
    mutations_per_division,
    fraction_mutated_intra_species, fraction_diverged_two_species,
)

from scripts.bio5_microbiology.growth import (
    mu_from_doubling, doubling_from_mu, mu_from_endpoints,
    cells_at_time, cells_at_time_doubling,
    inoculum_needed, synchronize_inoculation,
)


# ===================================================================
# BIO1 - PHYSICAL PRINCIPLES
# ===================================================================

# ---- diffusion ----
# print("D  =", stokes_einstein_D(R=2e-9, T=T_room, eta=eta_water_25C))   # R[m] T[K] eta[Pa*s] -> [m^2/s]
# print("x  =", diffusion_distance(D=1e-10, t=1, dim=1))                  # D[m^2/s] t[s] dim=1/2/3 -> [m]
# print("t  =", diffusion_time(x=10e-6, D=1e-10, dim=1))                  # x[m] D[m^2/s] dim=1/2/3 -> [s]
print("D  =", D_from_trajectory(x=6e-6, t=20, dim=2))                    # x[m] t[s] dim=1/2/3 -> D[m^2/s]
print("R  =", R_from_D(D=9.000000000000001e-13, T=T_room, eta=eta_water_25C))           # D[m^2/s] -> R[m]
# print("D2 =", scale_D(D1=1e-10, R1=2e-9, R2=4e-9, eta1=eta_water_25C, eta2=eta_water_25C, T1=T_room, T2=T_room))
# print("S  =", boltzmann_entropy(Omega=928))                             # microstates -> [J/K]

# ---- osmosis ----
# print("Pi =", osmotic_pressure(c_M=0.15, i=2, T=T_room))                # c[mol/L] i[vant Hoff] -> [Pa]
# print("osm=", effective_osmolarity(c_M=0.1, solute="NaCl"))             # -> [Osm/L]
# print(rank_osmotic([("A", 0.1, "NaCl"), ("B", 0.1, "glucose"), ("C", 0.05, "CaCl2")]))
# print(VANT_HOFF)                                                         # lookup table

# ---- flow ----
# print("Q  =", hagen_poiseuille_Q(dP=1000, r=10e-6, L=1e-3, eta=eta_water_25C))   # [m^3/s]
# print("dP =", hagen_poiseuille_dP(Q=1e-12, r=10e-6, L=1e-3, eta=eta_water_25C))  # [Pa]
# print("h  =", jurin_height(r=0.5e-3, theta_deg=0, gamma=gamma_water, rho=1000))  # [m]
# print("Re =", reynolds(rho=1000, u=0.01, L=1e-3, eta=eta_water_25C))             # dimensionless

# ---- membrane / cell ----
# print("A  =", sphere_area(R=10e-6))                                     # R[m] -> [m^2]
# print("V  =", sphere_volume(R=10e-6))                                   # R[m] -> [m^3]
# print("C  =", capacitance(A=1.26e-9, d=membrane_thickness, eps_r=eps_membrane))  # [F]
# print("Q  =", capacitor_charge(C=1e-11, V=0.07))                        # [C]
# print("E  =", capacitor_energy(C=1e-11, V=0.07))                        # [J]
# print("N  =", n_charges(Q=1e-12))                                       # count of e
# print("N  =", n_ions_in_cell(c_M=0.14, V_m3=4.19e-15))                  # ions in cell

# ---- energy ----
# print(atp_per_day(kcal_per_day=2000, fraction_to_atp=0.40))             # -> (mol, g)
# print("F  =", polymer_force(R=20e-9, L_contour=100e-9, l_persist=50e-9, T=T_room))
# print("e  =", fat_energy_density(kcal_per_g_fat=9.0, rho_fat=900))      # [J/m^3]
# print("e  =", atp_energy_density(c_ATP_M=5e-3))                         # [J/m^3]


# ===================================================================
# BIO2 - INFORMATION FLOW
# ===================================================================

# print("N  =", num_sequences_encoding("MAW", include_stop=True))         # protein -> # DNA sequences
# print("aa =", translate("AUGGCAUGGUAA"))                                # DNA/RNA -> protein
# print("L  =", min_codon_length(n_amino_acids=25, n_bases=4))            # alt genetic code
# print("S  =", library_size(n_positions=6, alphabet=20))                 # randomized library
# print("S  =", dna_sequence_space(length=30))                            # 4^L
# print(compare([("4^10", 4**10), ("20^6", 20**6), ("4^20", 4**20)]))
# print(bits_in_genome(genome_bp=3e9))                                    # -> (bits, bytes)
# print("P  =", orf_probability(n_codons=100, p_start=1/64, p_stop=3/64))
# print("N  =", expected_orfs(genome_bp=5e6, n_codons=100, both_strands=True))
# print("N  =", pcr_copies(n_cycles=30, efficiency=1.0, n_initial=1))
# print("n  =", cycles_for_copies(n_target=1e9, n_initial=1, efficiency=1.0))
# print("c  =", transcript_concentration(n_transcripts=10, V_cell_m3=1e-15))   # [mol/L]
# print("N  =", molecules_from_concentration(c_M=100e-9, V_m3=1e-15))


# ===================================================================
# BIO3 - PROTEINS
# ===================================================================

# print("k  =", rate_enhancement_kBT(dG_uncat_kBT=20, dG_cat_kBT=10))
# print("k  =", rate_enhancement_kJmol(dG_uncat_kJmol=120, dG_cat_kJmol=70, T=T_room))
# print("v  =", michaelis_menten(S_M=1e-3, Vmax=1.0, Km_M=1e-3))
# print("V  =", vmax_from_kcat(kcat_per_s=100, E_total_M=1e-6))
# print("th =", fractional_occupancy(L_M=1e-6, Kd_M=1e-6))
# print("th =", hill_occupancy(L_M=1e-6, Kd_M=1e-6, n=4))

# print("MW =", peptide_MW("ACDEFGHIKLMNPQRSTVWY"))                       # one-letter codes
# print("q  =", peptide_charge("ACDEFGHIKLMNPQRSTVWY", pH=7))
# print("mz =", peptide_mz("ACDEFGHIKLMNPQRSTVWY", z=2))
# print("HP =", hp_pattern("ACDEFGHIKLMNPQRSTVWY"))
# print("ok =", hp_matches("AVLI", "HHHH"))

# print(enzyme_mass_for_throughput(c_substrate_M=10e-3, V_m3=1e-3, kcat_per_s=100, t_s=3600, MW_enzyme=40_000))
# print(antibody_mass_for_toxin(mass_toxin_g=1e-3, MW_toxin=7000, MW_antibody=150_000, stoich_ab_per_toxin=2))
# print("m  =", catalyst_mass_for_production(moles_product=10, rate_mol_per_g_per_h=0.05, time_h=24))


# ===================================================================
# BIO4 - EVOLUTION AND ECOLOGY
# ===================================================================

# print("f  =", virus_volume_fraction(n_virus=1000, R_virus_m=50e-9, V_cell_m3=1e-15))
# print("f  =", virus_dry_mass_fraction(n_virus=1000, R_virus_m=50e-9, V_cell_m3=1e-15))
# print("n  =", replication_cycles(N_initial=1, N_final=1e9, fold_per_cycle=100))
# print("N  =", endogenous_retroviruses(host_genome_bp=3e9, retro_fraction=0.08, retrovirus_genome_bp=9000))
# print("m  =", mutations_per_division(genome_bp=3e9, error_rate_per_bp=1e-9))
# print("f  =", fraction_mutated_intra_species(mut_per_gen=3, T_years=1e6, gen_time_years=20, genome_bp=3e9))
# print("f  =", fraction_diverged_two_species(mut_per_gen=3, T_years=6e6, gen_time_years=20, genome_bp=3e9))


# ===================================================================
# BIO5 - MICROBIOLOGY
# ===================================================================

# print("mu =", mu_from_doubling(t_double=20))                             # any time unit
# print("td =", doubling_from_mu(mu=0.0347))
# print("mu =", mu_from_endpoints(N0=1e5, N1=8e5, t=90))
# print("N  =", cells_at_time(N0=1e6, mu=0.0347, t=60))
# print("N  =", cells_at_time_doubling(N0=1e6, t_double=20, t=60))
# print("N0 =", inoculum_needed(N_target=5e9, t_double=80, t_grow=24*60))
# print("t  =", synchronize_inoculation(N_target=1e9, N0=1e3, t_double=20, t_final=10*60))
