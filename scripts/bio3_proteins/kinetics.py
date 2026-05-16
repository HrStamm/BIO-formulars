import math
from scripts.constants import kBT_room, RT_room_kJ, T_room, kB, R_gas


def rate_enhancement_kBT(dG_uncat_kBT, dG_cat_kBT):
    # both in units of kBT -> dimensionless ratio
    return math.exp(dG_uncat_kBT - dG_cat_kBT)


def rate_enhancement_kJmol(dG_uncat_kJmol, dG_cat_kJmol, T=T_room):
    # both in kJ/mol -> dimensionless ratio
    RT = R_gas * T / 1000.0
    return math.exp((dG_uncat_kJmol - dG_cat_kJmol) / RT)


def michaelis_menten(S_M, Vmax, Km_M):
    # S [mol/L], Vmax [whatever / s], Km [mol/L] -> v
    return Vmax * S_M / (Km_M + S_M)


def vmax_from_kcat(kcat_per_s, E_total_M):
    # kcat [1/s], [E]_total [mol/L] -> Vmax [mol/(L*s)]
    return kcat_per_s * E_total_M


def fractional_occupancy(L_M, Kd_M):
    # L [mol/L], Kd [mol/L] -> theta (0..1)
    return L_M / (Kd_M + L_M)


def hill_occupancy(L_M, Kd_M, n):
    # cooperativity coefficient n
    return L_M**n / (Kd_M**n + L_M**n)


if __name__ == "__main__":
    print(f"enhancement (10 kBT lower) = {rate_enhancement_kBT(20, 10):.3e}")
    print(f"enhancement (50 kJ/mol lower) = {rate_enhancement_kJmol(120, 70):.3e}")
    print(f"v at S=Km: {michaelis_menten(1e-3, Vmax=1.0, Km_M=1e-3):.2f} (= Vmax/2)")
    print(f"theta at L=Kd: {fractional_occupancy(1e-6, 1e-6):.2f}")
    print(f"hill theta at L=Kd, n=4: {hill_occupancy(1e-6, 1e-6, 4):.2f}")
