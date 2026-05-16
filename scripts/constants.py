import math

kB = 1.380649e-23          # J/K
N_A = 6.02214076e23        # 1/mol
R_gas = 8.314462618        # J/(mol*K)
e_charge = 1.602176634e-19 # C
eps0 = 8.8541878128e-12    # F/m
g = 9.80665                # m/s^2

T_room = 298.15            # K (25 C)
T_body = 310.15            # K (37 C)
kBT_room = kB * T_room     # J  ~ 4.11e-21
kBT_body = kB * T_body     # J
RT_room_kJ = R_gas * T_room / 1000   # kJ/mol  ~ 2.479
RT_body_kJ = R_gas * T_body / 1000   # kJ/mol  ~ 2.577

eta_water_20C = 1.002e-3   # Pa*s
eta_water_25C = 0.890e-3   # Pa*s
eta_water_37C = 0.692e-3   # Pa*s

membrane_thickness = 5e-9  # m
eps_membrane = 5           # dimensionless
gamma_water = 0.0728       # N/m surface tension

MW_ATP = 507.18            # g/mol
kJ_per_mol_ATP_hydrolysis = 30.5  # kJ/mol (standard ~ -30.5, can be ~50 in cell)
kJ_per_kcal = 4.184
