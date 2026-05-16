def nm_to_m(x):      return x * 1e-9
def um_to_m(x):      return x * 1e-6
def mm_to_m(x):      return x * 1e-3
def cm_to_m(x):      return x * 1e-2
def angstrom_to_m(x):return x * 1e-10

def m_to_nm(x):      return x * 1e9
def m_to_um(x):      return x * 1e6

def pL_to_m3(x):     return x * 1e-15    # 1 pL = 1e-15 m^3
def fL_to_m3(x):     return x * 1e-18
def uL_to_m3(x):     return x * 1e-9
def mL_to_m3(x):     return x * 1e-6
def L_to_m3(x):      return x * 1e-3

def m3_to_L(x):      return x * 1e3

def mM_to_M(x):      return x * 1e-3
def uM_to_M(x):      return x * 1e-6
def nM_to_M(x):      return x * 1e-9
def pM_to_M(x):      return x * 1e-12

def M_to_per_m3(c_M):       return c_M * 1000.0   # mol/L -> mol/m^3
def per_m3_to_M(c_per_m3):  return c_per_m3 / 1000.0

def kcal_to_kJ(x):   return x * 4.184
def kJ_to_kcal(x):   return x / 4.184

def kDa_to_g_per_mol(x): return x * 1000.0
def g_per_mol_to_kDa(x): return x / 1000.0

def fN_to_N(x):      return x * 1e-15
def pN_to_N(x):      return x * 1e-12

def C_to_K(t_c):     return t_c + 273.15
def K_to_C(t_k):     return t_k - 273.15
