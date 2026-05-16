import math
from scripts.constants import kB, T_room, eta_water_25C


def stokes_einstein_D(R, T=T_room, eta=eta_water_25C):
    # R: radius [m], T: [K], eta: viscosity [Pa*s] -> D [m^2/s]
    return kB * T / (6 * math.pi * eta * R)


def diffusion_distance(D, t):
    # 1D rms: D [m^2/s], t [s] -> x [m]
    return math.sqrt(2 * D * t)


def diffusion_time(x, D):
    # x [m], D [m^2/s] -> t [s]
    return x * x / (2 * D)


def D_from_trajectory(x, t):
    # x [m], t [s] -> D [m^2/s]
    return x * x / (2 * t)


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

    # back-calc radius from trajectory
    D_meas = D_from_trajectory(x=5e-6, t=1.0)
    print(f"R from trajectory = {R_from_D(D_meas)*1e9:.2f} nm")

    # scale: same particle, viscosity doubles
    D2 = scale_D(D, R1=R, R2=R, eta2=2*eta_water_25C)
    print(f"D in 2x viscosity = {D2:.3e} m^2/s")
