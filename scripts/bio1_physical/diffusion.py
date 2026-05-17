import math
from scripts.constants import kB, T_room, eta_water_25C


def stokes_einstein_D(R, T=T_room, eta=eta_water_25C):
    # R: radius [m], T: [K], eta: viscosity [Pa*s] -> D [m^2/s]
    return kB * T / (6 * math.pi * eta * R)


def diffusion_distance(D, t, dim=1):
    # rms: D [m^2/s], t [s], dim = 1/2/3 -> x [m]
    # MSD = 2*dim*D*t
    return math.sqrt(2 * dim * D * t)


def diffusion_time(x, D, dim=1):
    # x [m], D [m^2/s], dim = 1/2/3 -> t [s]
    return x * x / (2 * dim * D)


def D_from_trajectory(x, t, dim=1):
    # x [m] = sqrt(MSD), t [s], dim = 1/2/3 -> D [m^2/s]
    # 1D: D = x^2 / (2t),  2D: D = x^2 / (4t),  3D: D = x^2 / (6t)
    return x * x / (2 * dim * t)


def R_from_D(D, T=T_room, eta=eta_water_25C):
    # D [m^2/s] -> R [m]
    return kB * T / (6 * math.pi * eta * D)


def scale_D(D1, R1, R2, eta1=eta_water_25C, eta2=eta_water_25C, T1=T_room, T2=T_room):
    # known D1 at (R1, eta1, T1) -> D2 at (R2, eta2, T2)
    return D1 * (R1 / R2) * (eta1 / eta2) * (T2 / T1)


def boltzmann_entropy(Omega):
    # Omega: number of microstates -> S [J/K]
    return kB * math.log(Omega)


if __name__ == "__main__":
    # protein R = 2 nm in water at 25 C
    R = 2e-9
    D = stokes_einstein_D(R)
    print(f"D = {D:.3e} m^2/s")
    print(f"x in 1 s = {diffusion_distance(D, 1):.3e} m")
    print(f"t to travel 10 um = {diffusion_time(10e-6, D):.3f} s")

    # back-calc radius from trajectory (2D tracking)
    D_meas = D_from_trajectory(x=5e-6, t=1.0, dim=2)
    print(f"R from trajectory (2D) = {R_from_D(D_meas)*1e9:.2f} nm")

    # scale: same particle, viscosity doubles
    D2 = scale_D(D, R1=R, R2=R, eta2=2*eta_water_25C)
    print(f"D in 2x viscosity = {D2:.3e} m^2/s")
