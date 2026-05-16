import math
from scripts.constants import (
    eps0, e_charge, N_A,
    membrane_thickness, eps_membrane,
)


def sphere_area(R):
    # R [m] -> A [m^2]
    return 4 * math.pi * R * R


def sphere_volume(R):
    # R [m] -> V [m^3]
    return (4.0 / 3.0) * math.pi * R**3


def capacitance(A, d=membrane_thickness, eps_r=eps_membrane):
    # A [m^2], d [m] -> C [F]
    return A * eps0 * eps_r / d


def capacitor_charge(C, V):
    # C [F], V [V] -> Q [C]
    return C * V


def capacitor_energy(C, V):
    # C [F], V [V] -> E [J]
    return 0.5 * C * V * V


def n_charges(Q):
    # Q [C] -> count of elementary charges
    return Q / e_charge


def n_ions_in_cell(c_M, V_m3):
    # c [mol/L], V [m^3] -> N ions
    V_L = V_m3 * 1000.0
    return c_M * V_L * N_A


if __name__ == "__main__":
    R = 10e-6                              # 10 um cell
    A = sphere_area(R)
    V = sphere_volume(R)
    C = capacitance(A)
    Vmem = -70e-3
    Q = capacitor_charge(C, abs(Vmem))
    E = capacitor_energy(C, abs(Vmem))
    nQ = n_charges(Q)
    nK = n_ions_in_cell(c_M=0.140, V_m3=V)
    print(f"A = {A:.3e} m^2, V = {V:.3e} m^3")
    print(f"C = {C:.3e} F, Q = {Q:.3e} C, E = {E:.3e} J")
    print(f"charges across membrane = {nQ:.3e}")
    print(f"K+ ions in cell        = {nK:.3e}")
    print(f"fraction on membrane   = {nQ/nK:.3e}")
