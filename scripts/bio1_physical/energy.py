import math
from scripts.constants import (
    kB, kBT_room, T_room, N_A,
    MW_ATP, kJ_per_mol_ATP_hydrolysis, kJ_per_kcal,
)


def atp_per_day(kcal_per_day, fraction_to_atp=0.40,
                kJ_per_mol_atp=kJ_per_mol_ATP_hydrolysis):
    # kcal_per_day, fraction usable, kJ/mol per ATP -> (moles, mass [g])
    kJ_total = kcal_per_day * kJ_per_kcal
    kJ_to_atp = kJ_total * fraction_to_atp
    moles = kJ_to_atp / kJ_per_mol_atp
    mass_g = moles * MW_ATP
    return moles, mass_g


def polymer_force(R, L_contour, l_persist, T=T_room):
    # R [m], L_contour [m], l_persist [m] -> F [N]
    l_kuhn = 2 * l_persist
    N_kuhn = L_contour / l_kuhn
    return 3 * kB * T * R / (N_kuhn * l_kuhn * l_kuhn)


def energy_density(E_total_J, V_m3):
    # E [J], V [m^3] -> J/m^3
    return E_total_J / V_m3


def fat_energy_density(kcal_per_g_fat=9.0, rho_fat=900.0):
    # kcal/g, rho [kg/m^3] -> J/m^3
    J_per_g = kcal_per_g_fat * kJ_per_kcal * 1000.0
    J_per_kg = J_per_g * 1000.0
    return J_per_kg * rho_fat


def atp_energy_density(c_ATP_M, kJ_per_mol_atp=kJ_per_mol_ATP_hydrolysis):
    # c_ATP [mol/L] -> J/m^3
    c_per_m3 = c_ATP_M * 1000.0
    return c_per_m3 * kJ_per_mol_atp * 1000.0


if __name__ == "__main__":
    moles, mass = atp_per_day(kcal_per_day=2000, fraction_to_atp=0.40)
    print(f"ATP/day: {moles:.1f} mol  ({mass/1000:.2f} kg)")

    F = polymer_force(R=20e-9, L_contour=100e-9, l_persist=50e-9)
    print(f"F (dsDNA, R=20 nm) = {F:.3e} N")

    e_fat = fat_energy_density()
    e_atp = atp_energy_density(c_ATP_M=5e-3)
    print(f"E_fat = {e_fat:.3e} J/m^3, E_ATP = {e_atp:.3e} J/m^3  -> ratio {e_fat/e_atp:.1f}")
