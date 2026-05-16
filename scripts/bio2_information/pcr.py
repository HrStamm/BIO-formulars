import math
from scripts.constants import N_A


def pcr_copies(n_cycles, efficiency=1.0, n_initial=1):
    # n_cycles, efficiency in (0,1], n_initial: starting molecules
    return n_initial * (1 + efficiency) ** n_cycles


def cycles_for_copies(n_target, n_initial=1, efficiency=1.0):
    # solve copies = n0 * (1+e)^n -> n
    return math.log(n_target / n_initial) / math.log(1 + efficiency)


def transcript_concentration(n_transcripts, V_cell_m3):
    # N transcripts in volume [m^3] -> c [mol/L]
    V_L = V_cell_m3 * 1000.0
    return n_transcripts / (V_L * N_A)


def molecules_from_concentration(c_M, V_m3):
    # c [mol/L], V [m^3] -> count
    V_L = V_m3 * 1000.0
    return c_M * V_L * N_A


if __name__ == "__main__":
    print(f"PCR 30 cycles from 1 molecule = {pcr_copies(30):.3e}")
    print(f"Cycles needed to reach 1e9 = {cycles_for_copies(1e9):.1f}")

    V = 1e-15   # 1 fL bacterium
    c = transcript_concentration(10, V)
    print(f"10 transcripts in 1 fL = {c*1e9:.2f} nM")
    print(f"100 nM in 1 fL = {molecules_from_concentration(100e-9, V):.1f} molecules")
