from scripts.constants import R_gas, T_room


VANT_HOFF = {
    "glucose": 1, "sucrose": 1, "fructose": 1, "urea": 1, "glycerol": 1,
    "NaCl": 2, "KCl": 2, "LiCl": 2, "NH4Cl": 2,
    "CaCl2": 3, "MgCl2": 3, "Na2SO4": 3, "K2SO4": 3,
    "AlCl3": 4, "FeCl3": 4,
    "Al2(SO4)3": 5,
}


def osmotic_pressure(c_M, i=1, T=T_room):
    # c_M [mol/L], i: van't Hoff factor, T [K] -> Pi [Pa]
    c_per_m3 = c_M * 1000.0
    return i * c_per_m3 * R_gas * T


def effective_osmolarity(c_M, solute):
    # solute: key in VANT_HOFF -> osmolarity [Osm/L]
    return VANT_HOFF[solute] * c_M


def rank_osmotic(solutions):
    # solutions: list of (label, c_M, solute) -> sorted (label, osmolarity, Pi)
    out = []
    for label, c, sol in solutions:
        osm = effective_osmolarity(c, sol)
        out.append((label, osm, osmotic_pressure(c, VANT_HOFF[sol])))
    return sorted(out, key=lambda r: r[1], reverse=True)


if __name__ == "__main__":
    print(f"0.15 M NaCl Pi = {osmotic_pressure(0.15, i=2)/1000:.2f} kPa")
    print(f"0.3 M glucose Pi = {osmotic_pressure(0.3, i=1)/1000:.2f} kPa")

    sols = [
        ("A: 0.1 M NaCl",       0.10, "NaCl"),
        ("B: 0.1 M glucose",    0.10, "glucose"),
        ("C: 0.05 M CaCl2",     0.05, "CaCl2"),
        ("D: 0.05 M Al2(SO4)3", 0.05, "Al2(SO4)3"),
    ]
    for label, osm, Pi in rank_osmotic(sols):
        print(f"{label:25s}  osm={osm:.3f}  Pi={Pi/1000:.2f} kPa")
