def enzyme_mass_for_throughput(c_substrate_M, V_m3, kcat_per_s, t_s, MW_enzyme):
    # c [mol/L], V [m^3], kcat [1/s], time [s], MW [g/mol] -> mass [g]
    V_L = V_m3 * 1000.0
    moles_substrate = c_substrate_M * V_L
    moles_enzyme = moles_substrate / (kcat_per_s * t_s)
    return moles_enzyme * MW_enzyme, moles_enzyme


def antibody_mass_for_toxin(mass_toxin_g, MW_toxin, MW_antibody, stoich_ab_per_toxin=1):
    # mass [g], MW [g/mol], stoich antibodies per toxin -> mass antibody [g]
    moles_toxin = mass_toxin_g / MW_toxin
    moles_ab = moles_toxin * stoich_ab_per_toxin
    return moles_ab * MW_antibody, moles_ab


def catalyst_mass_for_production(moles_product, rate_mol_per_g_per_h,
                                 time_h, MW_catalyst=None):
    # moles_product target, rate [mol/(g*h)], time [h] -> mass [g]
    mass = moles_product / (rate_mol_per_g_per_h * time_h)
    return mass


if __name__ == "__main__":
    mass, moles = enzyme_mass_for_throughput(
        c_substrate_M=10e-3, V_m3=1e-3, kcat_per_s=100, t_s=3600, MW_enzyme=40_000,
    )
    print(f"enzyme: {mass:.3f} g  ({moles*1e6:.2f} umol) for 10 mM in 1 L, 1 h, kcat=100")

    mass_ab, n_ab = antibody_mass_for_toxin(
        mass_toxin_g=0.001, MW_toxin=7000, MW_antibody=150_000, stoich_ab_per_toxin=2,
    )
    print(f"antibody: {mass_ab*1000:.2f} mg ({n_ab*1e9:.1f} nmol) for 1 mg toxin")

    m = catalyst_mass_for_production(moles_product=10, rate_mol_per_g_per_h=0.05, time_h=24)
    print(f"catalyst: {m:.2f} g for 10 mol in 24 h at 0.05 mol/(g*h)")
