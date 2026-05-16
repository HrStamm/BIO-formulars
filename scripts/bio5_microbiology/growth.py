import math


def mu_from_doubling(t_double):
    # t_double in any time unit -> mu in 1/(same unit)
    return math.log(2) / t_double


def doubling_from_mu(mu):
    return math.log(2) / mu


def mu_from_endpoints(N0, N1, t):
    # exponential growth between two timepoints
    return math.log(N1 / N0) / t


def cells_at_time(N0, mu, t):
    # N0 * exp(mu*t)
    return N0 * math.exp(mu * t)


def cells_at_time_doubling(N0, t_double, t):
    # N0 * 2^(t/t_double)
    return N0 * 2 ** (t / t_double)


def inoculum_needed(N_target, t_double, t_grow):
    # back-calculate N0
    return N_target / (2 ** (t_grow / t_double))


def synchronize_inoculation(N_target, N0, t_double, t_final):
    # time after t=0 to inoculate so culture reaches N_target at t_final
    mu = mu_from_doubling(t_double)
    t_growth = math.log(N_target / N0) / mu
    return t_final - t_growth


if __name__ == "__main__":
    mu = mu_from_doubling(20)           # 20 min doubling -> 1/min
    print(f"mu (20 min doubling) = {mu:.4f} /min")
    print(f"t_double (mu=0.5/h) = {doubling_from_mu(0.5):.2f} h")
    print(f"mu from N0=1e5 -> N1=8e5 in 90 min: {mu_from_endpoints(1e5, 8e5, 90):.4f} /min")
    print(f"N0=1e6, 60 min growth at 20 min doubling: {cells_at_time_doubling(1e6, 20, 60):.2e}")
    print(f"inoculum for 5e9 in 24 h, 80 min doubling: {inoculum_needed(5e9, 80, 24*60):.2e}")

    # synchronize: 3 strains, target 1e9 at t=10 h, N0=1e3 each
    for tdouble in (20, 30, 60):
        t_ino = synchronize_inoculation(N_target=1e9, N0=1e3,
                                        t_double=tdouble, t_final=10*60)
        print(f"strain t_d={tdouble} min -> inoculate at t = {t_ino:.1f} min")
