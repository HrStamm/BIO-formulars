import math


def sphere_volume(R):
    # R [m] -> V [m^3]
    return (4.0 / 3.0) * math.pi * R**3


def virus_volume_fraction(n_virus, R_virus_m, V_cell_m3):
    # n_virus per cell, R_virus [m], V_cell [m^3] -> fraction
    return n_virus * sphere_volume(R_virus_m) / V_cell_m3


def virus_dry_mass_fraction(n_virus, R_virus_m, V_cell_m3,
                            rho_virus=1300.0, rho_cell=1100.0,
                            dry_frac_cell=0.30, dry_frac_virus=1.0):
    # dry mass viruses / dry mass cell
    V_virus = n_virus * sphere_volume(R_virus_m)
    dry_virus = V_virus * rho_virus * dry_frac_virus
    dry_cell = V_cell_m3 * rho_cell * dry_frac_cell
    return dry_virus / dry_cell


def replication_cycles(N_initial, N_final, fold_per_cycle=100):
    # N_initial -> N_final at fold_per_cycle growth -> n cycles
    return math.log(N_final / N_initial) / math.log(fold_per_cycle)


def endogenous_retroviruses(host_genome_bp, retro_fraction,
                            retrovirus_genome_bp=9000):
    # fraction of host genome retroviral, ERV genome size -> count
    return host_genome_bp * retro_fraction / retrovirus_genome_bp


if __name__ == "__main__":
    V_cell = 1e-15        # 1 fL
    R_v = 50e-9
    print(f"vol fraction (1000 virus, 50 nm) = {virus_volume_fraction(1000, R_v, V_cell):.3e}")
    print(f"dry mass fraction               = {virus_dry_mass_fraction(1000, R_v, V_cell):.3e}")
    print(f"cycles 1->1e9 at 100x/cycle     = {replication_cycles(1, 1e9):.2f}")
    print(f"ERVs in 3e9 bp, 8% retroviral   = {endogenous_retroviruses(3e9, 0.08):.2e}")
