import math
from scripts.constants import g, gamma_water, eta_water_25C


def hagen_poiseuille_Q(dP, r, L, eta=eta_water_25C):
    # dP [Pa], r [m], L [m], eta [Pa*s] -> Q [m^3/s]
    return math.pi * dP * r**4 / (8 * eta * L)


def hagen_poiseuille_dP(Q, r, L, eta=eta_water_25C):
    # Q [m^3/s] -> dP [Pa]
    return 8 * eta * L * Q / (math.pi * r**4)


def jurin_height(r, theta_deg=0, gamma=gamma_water, rho=1000.0):
    # r [m], contact angle [deg], gamma [N/m], rho [kg/m^3] -> h [m]
    theta = math.radians(theta_deg)
    return 2 * gamma * math.cos(theta) / (rho * g * r)


def reynolds(rho, u, L, eta=eta_water_25C):
    # rho [kg/m^3], u [m/s], L [m] -> Re
    return rho * u * L / eta


if __name__ == "__main__":
    Q = hagen_poiseuille_Q(dP=1000, r=10e-6, L=1e-3)
    print(f"Q (10 um capillary, 1 kPa, 1 mm) = {Q:.3e} m^3/s")
    Q_half = hagen_poiseuille_Q(dP=1000, r=5e-6, L=1e-3)
    print(f"Q with half radius = {Q_half:.3e} m^3/s (ratio {Q/Q_half:.0f}x)")

    h = jurin_height(r=0.5e-3, theta_deg=0)
    print(f"Jurin h (0.5 mm radius) = {h*1000:.2f} mm")

    Re = reynolds(rho=1000, u=0.01, L=1e-3)
    print(f"Re = {Re:.2f}")
