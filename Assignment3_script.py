# --- ASSIGNMENT 3 ---
import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
import math
from scipy.special import j0
from scipy.optimize import fsolve
import numpy as np

# Complete vectors of the inner and outer cladding temperatures (from Assignment 2)
T_ci_vec = np.array([
    280.79350629, 283.91485079, 287.06656817, 290.24524538, 293.44745721,
    296.66977038, 299.90874758, 303.16095129, 306.4229476 , 309.69130975,
    312.96262162, 316.23348103, 319.50050296, 322.76032258, 326.00959815,
    329.24501385, 332.46328237, 335.66114755, 338.83538674, 341.98281316,
    345.10027818, 348.18467339, 351.23293278, 354.24203468, 357.20900377,
    360.130913  , 363.00488543, 365.82809618, 368.5977742 , 371.31120418,
    373.9657284 , 376.55874863, 379.08772801, 381.55019305, 383.9437356 ,
    386.26601493, 388.51475982, 390.68777076, 392.78292221, 394.7981649 ,
    396.73152828, 398.58112294, 400.34514322, 402.02186977, 403.60967231,
    405.1070123 , 405.77787877, 405.97489394, 406.11573303, 406.20028208,
    406.2284727 , 406.20028208, 406.11573303, 405.97489394, 405.77787877,
    405.52484697, 405.21600339, 404.85159817, 404.43192658, 403.9573289 ,
    403.42819015, 402.84493992, 402.20805208, 401.5180445 , 400.77547872,
    399.98095961, 399.13513492, 398.23869491, 397.29237182, 396.29693933,
    395.25321203, 394.16204468, 393.02433162, 391.84100587, 390.61303833,
    389.34143679, 388.02724485, 386.67154068, 385.27543567, 383.8400728 ,
    382.36662487, 380.85629232, 379.31030073, 377.7298978 , 376.11634964,
    374.47093621, 372.79494552, 371.08966621, 369.35637755, 367.59633606,
    365.8107566 , 364.00078506, 362.16745715, 360.31163345, 358.43389132,
    356.5343321 , 354.61220638, 352.66509335, 350.68676447, 348.65976526,
    346.09862569
])

T_co_vec = np.array([
    279.88988559, 281.01565898, 282.18255964, 283.3892164 , 284.63420518,
    285.91605106, 287.23323052, 288.58417364, 289.96726641, 291.38085297,
    292.82323798, 294.29268889, 295.78743832, 297.30568642, 298.84560317,
    300.40533077, 301.98298596, 303.57666243, 305.1844331 , 306.80435253,
    308.43445923, 310.07277808, 311.71732264, 313.3660976 , 315.01710113,
    316.66832733, 318.31776869, 319.96341858, 321.60327377, 323.23533701,
    324.85761964, 326.46814428, 328.06494758, 329.64608299, 331.20962367,
    332.75366544, 334.27632976, 335.77576685, 337.25015889, 338.69772325,
    340.1167158 , 341.50543433, 342.86222202, 344.18547093, 345.47362557,
    346.72518656, 347.16675553, 347.1687235 , 347.17012749, 347.17096919,
    347.17124965, 347.17096919, 347.17012749, 347.1687235 , 347.16675553,
    347.16422116, 347.16111724, 347.15743989, 347.15318446, 347.14834548,
    347.14291663, 347.13689068, 347.13025944, 347.12301368, 347.11514306,
    347.10663603, 347.09747973, 347.08765984, 347.07716047, 347.06596397,
    347.05405072, 347.04139896, 347.02798446, 347.01378022, 346.99875614,
    346.98287855, 346.96610966, 346.94840699, 346.92972259, 346.91002208,
    346.88918359, 346.86719628, 346.84395863, 346.81937621, 346.79333885,
    346.76571693, 346.73635657, 346.70507312, 346.67164229, 346.6357878 ,
    346.59716369, 346.55532826, 346.5097042 , 346.45951522, 346.40367962,
    346.3406196 , 346.26788916, 346.18135737, 346.07307868, 345.9238899 ,
    345.24637571
])

# --- AP1000 DATA ---
# === POWER AND HEAT FLUX ===
reactor_type = "PWR"
coolant = "Light water"
moderator = "Light water"

thermal_power = 3400 # Reactor core heat output [MWt]
heat_in_fuel = 0.974  # [-] Fraction of heat generated in fuel (97.4%)
q_flux_avg     = 628710 * 10**-6  # Average heat flux [MW/m^2] (convertito da 199,300 Btu/hr-ft^2)
q_flux_max     = 1634709  # Maximum heat flux for normal op. [W/m^2] (da 518,200 Btu/hr-ft^2)
q_lin_avg      = 18766  # Average linear power [W/m] (convertito da 5.72 kW/ft)
q_lin_max      = 48884  # Peak linear power for normal op. [W/m] (da 14.9 kW/ft)
FQ = 2.60
power_density  = 109.7e6  # Power density [W/m^3] (convertito da 109.7 kW/liter)
g = 9.81

# === PRESSURE ===
P = 155.13 # System pressure, nominal [bar] (convertito da 2250 psia)
P_Pa = P * 1e5  # System pressure, nominal [Pa]
min_pressure = 150.99 # System pressure, minimum steady-state [bar] (da 2190 psia)
dp_core_bar    = 2.75  # Pressure drop across core [bar] (convertito da 39.9 psi)
dp_vessel_bar  = 4.30  # Pressure drop across vessel [bar] (convertito da 62.3 psi)

# === CORE and ASSEMBLY GEOMETRY ===
num_fuel_assemblies = 157  # Number of fuel assemblies [-]
core_diameter= 3.040 # m
core_height = 4.2672  # Core height, cold, active fuel [m] (convertito da 168.0 in)

# === FUEL ASSEMBLY ===
fuel_assembly_type = "17x17 XL Robust"
fuel_rods_per_assembly = 264  # Uranium dioxide rods per assembly [-]
rod_pitch = 0.01260  # Rod pitch [m] (convertito da 0.496 in)
assembly_width = 0.21402  # Overall dimensions (square) [m] (convertito da 8.426 in)
fuel_weight = 95974.70238 # kg
clad_weight = 19552 # kg

# Grids
num_grids_top_bottom = 2
num_grids_intermediate = 8
num_grids_IFM = 4
num_grids_protective = 1

# === FUEL RODS ===
total_rods = 41448
rod_diameter = 0.00950 # External cladding diameter (D_co) [m]
gap = 0.000165 # m
clad_thickness = 0.0005715  # Cladding thickness (t) [m]

# === FUEL PELLETS ===
pellet_material = "UO2"
pellet_density_percent = 0.955
pellet_diameter = 0.00819 # m
pellet_length = 0.00983 # m

# === COOLANT FLOW ===
total_flow_kg_s = 14300 # kg/s
effective_flow_kg_s = 13460 # kg/s
flow_area_m2 = 3.883 # m^2
velocity_m_s = 4.82 # m/s
mass_velocity = 3460 # kg/(s*m^2)

# === TEMPERATURES ===
inlet_temp_C = 279.4 # °C
rise_vessel_C = 42.9 # °C
rise_core_C = 45.2 # °C
avg_core_temp_C = 303.4 # °C
avg_vessel_temp_C = 300.89 # °C

# === HEAT TRANSFER ===
heat_area_m2 = 5268 # m^2
max_heat_flux_W_m2 = 1635000 # W/m^2
peak_linear_power_kW_m = 48.88 # kW/m
max_peak_linear_power_kW_m = 73.66 # kW/m

# === CONTROL RODS ===
RCCA_clusters = 53
GRCA_clusters = 16
absorber_material = "Ag-In-Cd"
rodlets_per_cluster = 24

# === SAFETY LIMITS ===
max_clad_temp_C = 1204 # °C
max_fuel_temp_C = 2593 # °C
max_control_speed_m_min = 1.143 # m/min

T_cl_max       = 2593.33  # Peak fuel center line temp. prevention of melt [°C] (da 4700 °F)
DNBR_nom_typ   = 2.80  # Minimum DNBR at nominal conditions (Typical flow channel) [-]
DNBR_nom_thim  = 2.74  # Minimum DNBR at nominal conditions (Thimble cold wall) [-]
void_avg       = 0.0  # Core average void fraction [%]
void_hot_avg   = 0.1  # Hot subchannel average void fraction [%]
void_hot_max   = 0.9  # Hot subchannel maximum void fraction [%]

# === STRENGHT ===
sigma_y = 241  # Yield strength [MPa]
sigma_u = 413  # Ultimate tensile strength [MPa]

sigma_y_Pa = sigma_y * 1e6  # Conversion in [Pa]
sigma_u_Pa = sigma_u * 1e6  # Conversion in [Pa]

# =============================================================================
# TASK 1: PRELIMINARY VERIFICATION AGAINST BUCKLING
# =============================================================================
print("\n --- TASK 1: BUCKLING VERIFICATION ---")
# Geometric parameters of the cladding
r_o = rod_diameter / 2  # Outer radius [m]
r_i = r_o - clad_thickness  # Inner radius [m]
r_avg = (r_o + r_i) / 2  # Average radius of the cladding [m]
t = clad_thickness  # cladding thickness [m]

z_axis = np.linspace(-2.1336, 2.1336, len(T_ci_vec))

# Functions for mechanical properties of Zircaloy-4 (T in Kelvin)
Th_exp_coeff_zir = lambda T: (6.72e-6 * T - 2.07e-3) / (T - 308)
Young_modulus_zir = lambda T: 9.81 * (9.9e3 - 5.669 * (T - 273.15))  # Young's Modulus [MPa]
Poisson_ratio_zir = lambda T: 0.3303 + 8.376e-5 * (T - 273.15)  # Poisson's Modulus [-]

T_cl_in_K = T_ci_vec + 273.15
T_cl_out_K = T_co_vec + 273.15
T_cl_ave_K = (T_cl_in_K + T_cl_out_K) / 2.0

E_cl_in = Young_modulus_zir(T_cl_in_K)  # Young's modulus of the inner wall [MPa]
E_cl_out = Young_modulus_zir(T_cl_out_K)  # Young's modulus of the outer wall [MPa]
E_cl_ave = Young_modulus_zir(T_cl_ave_K)

poisson_cl_in = Poisson_ratio_zir(T_cl_in_K)  # Poisson's ratio of the inner wall [-]
poisson_cl_out = Poisson_ratio_zir(T_cl_out_K)  # Poisson's ratio of the outer wall [-]
poisson_cl_ave = Poisson_ratio_zir(T_cl_ave_K)

# critical pressure (p_cr) for buckling
# We multiply by 10 because 1 MPa = 10 bar
p_cr_in_bar = (E_cl_in / (4 * (1 - poisson_cl_in**2))) * (t / r_avg)**3 * 10
p_cr_out_bar = (E_cl_out / (4 * (1 - poisson_cl_out**2))) * (t / r_avg)**3 * 10
p_cr_bar_ave = (E_cl_ave / (4 * (1 - poisson_cl_ave**2))) * (t / r_avg)**3 * 10

#min_p_cr_bar = min(np.min(p_cr_in_bar), np.min(p_cr_out_bar))
min_p_cr_bar = np.min(p_cr_bar_ave)

idx_min_buckling = np.argmin(p_cr_bar_ave)

print(f"Average radius of the cladding (r_avg): {r_avg * 1000:.4f} mm")
print(f"Critical pressure for Buckling (p_cr): {min_p_cr_bar:.2f} bar  (localized at z = {z_axis[idx_min_buckling]:.2f} m)")
print(f"External system pressure (P): {P:.2f} bar")

# Check against external system pressure (P)
if min_p_cr_bar > P:
    print(
        f"CHECK PASSED: The sheath is SAFE against buckling (p_cr_min > P). "
        f"The difference between the two pressures is: {min_p_cr_bar - P:.2f} bar. ")
else:
    print("FAILED CHECK: Risk of elastic instability / buckling (p_cr_min <= P)."
          f"Deficit: {P - min_p_cr_bar:.2f} bar. ")


# =============================================================================
# TASK 2: MAXIMUM INTERNAL PRESSURE (depressurized reactor)
# =============================================================================
print("\n --- TASK 2: MAXIMUM INTERNAL PRESSURE (Thick-wall Lamé formulation) ---")

ratio = t / r_avg
print(f"t / r_avg : {ratio:.2f}")
if ratio <= 0.1:
    print("The cladding is classified as THIN WALL (t/r_avg <= 0.1)")
else:
    print("The cladding is classified as a THICK WALL (t/r_avg > 0.1)")
print("\n")

# p_i_max = p_i
p_o = 0.0      # [Pa]  Depressurized reactor (outer pressure ~ 0)

# Rigorous evaluation of maximum internal pressure
p_i_max_MPa = sigma_y * (r_o**2 - r_i**2) / (2 * r_o**2)
p_i_max_bar = p_i_max_MPa * 10
p_i_max_Pa = p_i_max_MPa * 1e6

sigma_h = p_i_max_MPa * (r_o**2 + r_i**2) / (r_o**2 - r_i**2)  # Tangential (hoop) stress [MPa]
sigma_l = p_i_max_MPa * r_i**2 / (r_o**2 - r_i**2)  # Longitudinal stress [MPa]
sigma_r = -p_i_max_MPa    # Radial stress [MPa]

print(f"Yield Strength (Limit for inner hoop stress): {sigma_y} MPa")
print(f"Longitudinal stress: {sigma_l:.2f} MPa")
print(f"radial stress: {sigma_r:.2f} MPa")
print(f"Rigorous maximum internal pressure (Lamé): {p_i_max_bar:.2f} bar ({p_i_max_MPa:.2f} MPa)")


# ==========================================================================================================================================================
# TASK 3: MECHANICAL STRESSES   (THICK WALLS) -->  for both the inner and outer cladding wall, and for each stress component
# ==========================================================================================================================================================
print("\n --- TASK 3: MECHANICAL STRESSES ---")

# Let's define the external refrigerant pressure under operating conditions (not hp depressurized reactor anymore)
p_o_operating = P / 10  # P = 155.13 bar -> 15.513 MPa

# --- CALCULATION OF THE SCALAR VALUES (LAME'S EQUATIONS for thick walls) ---
# Common denominator term of the Lamé equations
den = r_o**2 - r_i**2

# Inner cladding wall (r = r_i)
sigma_h_mec_in = (p_i_max_MPa * (r_o**2 + r_i**2) - 2.0 * p_o_operating * r_o**2) / den
sigma_r_mec_in = -p_i_max_MPa
sigma_l_mec_in = (p_i_max_MPa * r_i**2 - p_o_operating * r_o**2) / den

# Outer cladding wall (r = r_o)
sigma_h_mec_out = (2 * p_i_max_MPa * r_i**2 - p_o_operating * (r_o**2 + r_i**2)) / den
sigma_r_mec_out = - p_o_operating
sigma_l_mec_out = (p_i_max_MPa * r_i**2 - p_o_operating * r_o**2) / den

print(f"Operating External Pressure (p_o): {P:.2f} bar ({p_o_operating:.3f} MPa)")
print(f"Internal Gas Pressure (p_i): {p_i_max_bar:.2f} bar ({p_i_max_MPa:.3f} MPa) \n")

print("Summary table with mechanical stresses considering THICK WALLS (LAME'S EQUATIONS)")
print(f"{'Component':<20} | {'Inner cladding wall [MPa]':<26} | {'Outer cladding wall [MPa]':<26}")
print("-" * 78)
print(f"{'Hoop (tangential)':<20} | {sigma_h_mec_in:<26.2f} | {sigma_h_mec_out:<26.2f}")
print(f"{'Radial':<20} | {sigma_r_mec_in:<26.2f} | {sigma_r_mec_out:<26.2f}")
print(f"{'Longitudinal':<20} | {sigma_l_mec_in:<26.2f} | {sigma_l_mec_out:<26.2f}\n")

# --- CALCULATION OF LOCAL VECTOR PROFILES ALONG THE ENTIRE Z AXIS ---
sigma_h_mec_in_vec = sigma_h_mec_in * np.ones_like(z_axis)
sigma_r_mec_in_vec = sigma_r_mec_in * np.ones_like(z_axis)
sigma_l_mec_in_vec = sigma_l_mec_in * np.ones_like(z_axis)

sigma_h_mec_out_vec = sigma_h_mec_out * np.ones_like(z_axis)
sigma_r_mec_out_vec = sigma_r_mec_out * np.ones_like(z_axis)
sigma_l_mec_out_vec = sigma_l_mec_out * np.ones_like(z_axis)


# ==========================================================================================================================================================
# TASK 3.1: PRIMARY STRESSES
# ==========================================================================================================================================================
print("\n --- PRIMARY STRESSES ---")
# Evaluation of the average mechanical stresses (primary stresses)
# Primary stresses have 3 components: membrane, bending and local, but we neglet the last two components

# --- CALCULATION OF THE AVERAGE SCALAR VALUES ---
sigma_h_mec_avg = (sigma_h_mec_in + sigma_h_mec_out) / 2
sigma_r_mec_avg = (sigma_r_mec_in + sigma_r_mec_out) / 2
sigma_l_mec_avg = (sigma_l_mec_in + sigma_l_mec_out) / 2

# --- CALCULATION OF LOCAL VECTOR PROFILES ALONG THE ENTIRE Z AXIS ---
sigma_h_mec_avg_vec = sigma_h_mec_avg * np.ones_like(z_axis)
sigma_r_mec_avg_vec = sigma_r_mec_avg * np.ones_like(z_axis)
sigma_l_mec_avg_vec = sigma_l_mec_avg * np.ones_like(z_axis)

print('Summary table with average mechanical stresses (PRIMARY STRESSES)')
print(f"Average tangential (hoop) stress: {sigma_h_mec_avg:.2f} MPa")
print(f"Average radial stress: {sigma_r_mec_avg:.2f} MPa")
print(f"Average longitudinal stress: {sigma_l_mec_avg:.2f} MPa ")


# ==========================================================================================================================================================
# TASK 4: THERMAL STRESSES  (THICK WALLS) -->  for both the inner and outer cladding wall, and for each stress component
# ==========================================================================================================================================================
print("\n --- TASK 4: THERMAL STRESSES ---")

# --- CALCULATION OF LOCAL VECTOR PROFILES ALONG THE ENTIRE Z AXIS ---

# Material properties calculated locally point by point
E_in_vec = Young_modulus_zir(T_ci_vec + 273.15)
E_out_vec = Young_modulus_zir(T_co_vec + 273.15)
nu_in_vec = Poisson_ratio_zir(T_ci_vec + 273.15)
nu_out_vec = Poisson_ratio_zir(T_co_vec + 273.15)
alpha_in_vec = Th_exp_coeff_zir(T_ci_vec + 273.15)
alpha_out_vec = Th_exp_coeff_zir(T_co_vec + 273.15)

delta_T_cl = np.array(T_ci_vec) - np.array(T_co_vec)
ln_ratio = np.log(r_o / r_i)
den_raggi = r_o**2 - r_i**2

# Analytical thermal factors (Thick Walls) alog the z-axis
factor_th_in_vec = (alpha_in_vec * E_in_vec * delta_T_cl) / (2.0 * (1.0 - nu_in_vec) * ln_ratio)
factor_th_out_vec = (alpha_out_vec * E_out_vec * delta_T_cl) / (2.0 * (1.0 - nu_out_vec) * ln_ratio)

# Lamé geometric terms for thick walls
term_in = (2.0 * r_o**2 / den_raggi) * ln_ratio - 1.0
term_out = 1.0 - (2.0 * r_i**2 / den_raggi) * ln_ratio

sigma_h_th_in_vec = -factor_th_in_vec * term_in
sigma_l_th_in_vec = -factor_th_in_vec * term_in
sigma_r_th_in_vec = np.zeros_like(z_axis)

sigma_h_th_out_vec = factor_th_out_vec * term_out  # for the external wall we use the positive sign since the fibres are taut
sigma_l_th_out_vec = factor_th_out_vec * term_out
sigma_r_th_out_vec = np.zeros_like(z_axis)


# --- CALCULATION OF THE THERMAL STRESSES IN THE PEAK POINTS ---

# WE EXTRACT THE PEAK POINTS
max_delta_T_cl = np.max(delta_T_cl)
indice_dt_max = np.argmax(delta_T_cl)

T_peak_in_K = T_ci_vec[indice_dt_max] + 273.15
T_peak_out_K = T_co_vec[indice_dt_max] + 273.15

# Evaluation of the mechanical and thermal properties of Zircaloy-4
E_in_peak = Young_modulus_zir(T_peak_in_K)          # [MPa]
E_out_peak = Young_modulus_zir(T_peak_out_K)        # [MPa]
nu_in_peak = Poisson_ratio_zir(T_peak_in_K)          # [-]
nu_out_peak = Poisson_ratio_zir(T_peak_out_K)        # [-]
alpha_in_peak = Th_exp_coeff_zir(T_peak_in_K)        # [1/K] Average linear thermal expansion coefficient of Zircaloy-4 [1/K]
alpha_out_peak = Th_exp_coeff_zir(T_peak_out_K)      # [1/K]


sigma_h_th_in = sigma_h_th_in_vec[indice_dt_max]
sigma_l_th_in = sigma_l_th_in_vec[indice_dt_max]
sigma_r_th_in = 0.0

sigma_h_th_out = sigma_h_th_out_vec[indice_dt_max]
sigma_l_th_out = sigma_l_th_out_vec[indice_dt_max]
sigma_r_th_out = 0.0

# Summary table with thermal stresses based on peak thermal conditions
print('Summary table with thermal stresses')
print(f"{'Component':<20} | {'Inner cladding wall [MPa]':^26} | {'Outer cladding wall [MPa]':^26}")
print("-" * 78)
print(f"{'Hoop (tangential)':<20} | {sigma_h_th_in:^26.2f} | {sigma_h_th_out:^26.2f}")
print(f"{'Radial':<20} | {sigma_r_th_in:^26.2f} | {sigma_r_th_out:^26.2f}")
print(f"{'Longitudinal':<20} | {sigma_l_th_in:^26.2f} | {sigma_l_th_out:^26.2f} ")


# =============================================================================
# TASK 4.1: SECONDARY STRESSES (mechanical and thermal components)
# =============================================================================
print("\n --- SECONDARY STRESSES ---")

# --- CALCULATION OF LOCAL VECTOR PROFILES ALONG THE ENTIRE Z AXIS ---

# SECONDARY STRESSES composed by mechanical and thermal components
sigma_sec_h_in_vec = (sigma_h_mec_in_vec - sigma_h_mec_avg_vec) + sigma_h_th_in_vec
sigma_sec_r_in_vec = (sigma_r_mec_in_vec - sigma_r_mec_avg_vec) + sigma_r_th_in_vec
sigma_sec_l_in_vec = (sigma_l_mec_in_vec - sigma_l_mec_avg_vec) + sigma_l_th_in_vec

sigma_sec_h_out_vec = (sigma_h_mec_out_vec - sigma_h_mec_avg_vec) + sigma_h_th_out_vec
sigma_sec_r_out_vec = (sigma_r_mec_out_vec - sigma_r_mec_avg_vec) + sigma_r_th_out_vec
sigma_sec_l_out_vec = (sigma_l_mec_out_vec - sigma_l_mec_avg_vec) + sigma_l_th_out_vec

# --- CALCULATION OF THE SECONDARY STRESSES IN THE PEAK POINTS ---
sigma_sec_h_in = sigma_sec_h_in_vec[indice_dt_max]
sigma_sec_r_in = sigma_sec_r_in_vec[indice_dt_max]
sigma_sec_l_in = sigma_sec_l_in_vec[indice_dt_max]

sigma_sec_h_out = sigma_sec_h_out_vec[indice_dt_max]
sigma_sec_r_out = sigma_sec_r_out_vec[indice_dt_max]
sigma_sec_l_out = sigma_sec_l_out_vec[indice_dt_max]

# Summary table with secondary stresses
print('Summary table with secondary stresses')
print(f"{'Component':<20} | {'Inner cladding wall [MPa]':^26} | {'Outer cladding wall [MPa]':^26}")
print("-" * 78)
print(f"{'Hoop (tangential)':<20} | {sigma_sec_h_in:^26.2f} | {sigma_sec_h_out:^26.2f}")
print(f"{'Radial':<20} | {sigma_sec_r_in:^26.2f} | {sigma_sec_r_out:^26.2f}")
print(f"{'Longitudinal':<20} | {sigma_sec_l_in:^26.2f} | {sigma_sec_l_out:^26.2f} ")


# =============================================================================
# TASK 4.2: PRIMARY + SECONDARY STRESSES AND COMPLETE SUMMARY TABLE
# =============================================================================
print("\n --- PRIMARY + SECONDARY STRESSES ---")

# TOTAL (PRIMARY + SECONDARY) STRESSES CALCULATION (Algebraic sum of the components)
# Average mechanical stresses cancel out

# --- CALCULATION OF LOCAL VECTOR PROFILES ALONG THE ENTIRE Z AXIS ---
# Inner cladding wall (r = r_i)
sigma_h_tot_in_vec = sigma_h_mec_in_vec + sigma_h_th_in_vec
sigma_r_tot_in_vec = sigma_r_mec_in_vec + sigma_r_th_in_vec
sigma_l_tot_in_vec = sigma_l_mec_in_vec + sigma_l_th_in_vec

# Outer cladding wall (r = r_o)
sigma_h_tot_out_vec = sigma_h_mec_out_vec + sigma_h_th_out_vec
sigma_r_tot_out_vec = sigma_r_mec_out_vec + sigma_r_th_out_vec
sigma_l_tot_out_vec = sigma_l_mec_out_vec + sigma_l_th_out_vec

# --- CALCULATION OF THE PRIMARY + SECONDARY STRESSES IN THE PEAK POINTS ---
sigma_h_tot_in = sigma_h_tot_in_vec[indice_dt_max]
sigma_r_tot_in = sigma_r_tot_in_vec[indice_dt_max]
sigma_l_tot_in = sigma_l_tot_in_vec[indice_dt_max]

sigma_h_tot_out = sigma_h_tot_out_vec[indice_dt_max]
sigma_r_tot_out = sigma_r_tot_out_vec[indice_dt_max]
sigma_l_tot_out = sigma_l_tot_out_vec[indice_dt_max]

# Summary table with primary + secondary stresses
print('Summary table with primary + secondary stresses')
print(f"{'Component':<25} | {'Inner Cladding Wall [MPa]':^25} | {'Outer Cladding Wall [MPa]':^25}")
print("-" * 81)
print(f"{'Hoop (tangential)':<25} | {sigma_h_tot_in:^25.2f} | {sigma_h_tot_out:^25.2f}")
print(f"{'Radial':<25} | {sigma_r_tot_in:^25.2f} | {sigma_r_tot_out:^25.2f}")
print(f"{'Longitudinal':<25} | {sigma_l_tot_in:^25.2f} | {sigma_l_tot_out:^25.2f}")


# =============================================================================
# TASK 5: ASME BPVC SECTION III VERIFICATION
# =============================================================================
print("\n --- TASK 5: ASME SECTION III VERIFICATION ---")

# Allowable stress according to the ASME code (sigma_a)
S_m = min((2/3) * sigma_y_Pa, (1/3) * sigma_u_Pa)  # [Pa] (for primary membrane stresses sigma_a = S_m)
S_m_MPa = S_m / 1e6  # [MPa]

S_m_2 = 3 * S_m  # [Pa] (for primary + secondary stresses sigma_a = 3 * S_m)
S_m_2_MPa = 3 * S_m_MPa  # [MPa]

print(f"Allowable stresses defined by ASME BPVC Section III (sigma_a):")
print(f"For primary membrane stresses (S_m): {S_m_MPa:.2f} MPa")
print(f"For primary + secondary stresses (3*S_m): {S_m_2_MPa:.2f} MPa\n")

# Equivalent stress with the Tresca’s criterion (PRIMARY STRESSES = Average mechanical stresses)
tresca_pri = max(abs(sigma_h_mec_avg - sigma_r_mec_avg), abs(sigma_r_mec_avg - sigma_l_mec_avg), abs(sigma_h_mec_avg - sigma_l_mec_avg))

# Primary + Thermal stresses
# Equivalent stress with the Tresca’s criterion, vector profiles (PRIMARY + SECONDARY STRESSES)
tresca_tot_in_vec = np.max([
    abs(sigma_h_tot_in_vec - sigma_r_tot_in_vec),
    abs(sigma_r_tot_in_vec - sigma_l_tot_in_vec),
    abs(sigma_h_tot_in_vec - sigma_l_tot_in_vec)
], axis=0)

tresca_tot_out_vec = np.max([
    abs(sigma_h_tot_out_vec - sigma_r_tot_out_vec),
    abs(sigma_r_tot_out_vec - sigma_l_tot_out_vec),
    abs(sigma_h_tot_out_vec - sigma_l_tot_out_vec)
], axis=0)

tresca_tot_in = tresca_tot_in_vec[indice_dt_max]
tresca_tot_out = tresca_tot_out_vec[indice_dt_max]

status_pri = 'FAILED' if (tresca_pri > S_m_MPa) else 'PASSED'
status_tot = 'FAILED' if (tresca_tot_in > S_m_2_MPa or tresca_tot_out > S_m_2_MPa) else 'PASSED'

print(f"{'Component (TRESCA)':<35} | {'Inner wall [MPa]':^20} | {'Outer wall [MPa]':^20} | {'Limit [MPa]':^15} | {'Status':^10}")
print("-" * 108)
print(f"{'Primary (S_m)':<35} | {tresca_pri:^43.2f} | {S_m_MPa:^15.2f} | {status_pri:^10}")
print(f"{'Primary + Secondary (3*S_m)':<35} | {tresca_tot_in:^20.2f} | {tresca_tot_out:^20.2f} | {S_m_2_MPa:^15.2f} | {status_tot:^10}")

# If the check fails, we recalculate the maximum internal pressure allowed
if tresca_pri > S_m_MPa:
    print("\n[!] NOTICE: Primary stress limit exceeded using the bounding pressure from Task 2.")

    # Analytical inversion of Lamé-Tresca under pressure: we obtain the maximum p_i so that Tresca = S_m
    p_i_max_allowed_MPa = (S_m_MPa * (r_o**2 - r_i**2) + 2.0 * p_o_operating * r_o**2) / (2.0 * r_o**2)
    p_x = p_i_max_allowed_MPa  # [MPa] DEFINITION OF p_x FOR STEP 6 (Converted to MPa consistently with p_plenum)
    print(f"--> Maximum allowable internal pressure to satisfy ASME limits: {p_x:.2f} MPa")

else:
    print("\n--> CHECK PASSED: Primary stress limits are satisfied.")
    # If the check passes, p_x coincides with the maximum pressure calculated in Task 2
    p_x = p_i_max_MPa
    print(f"--> p_x set to Step 2 maximum pressure: {p_x * 10:.2f} bar ({p_x:.2f} MPa)")


# =============================================================================
# TASK 6: SIZE THE GAS PLENUM OF THE FUEL PIN
# =============================================================================
print("\n --- TASK 6: GAS PLENUM SIZING ---")

# DATA
Burnup = 60000  # Reactor burn-up [MWD/t_U]
Burnup_Jt = Burnup * 1e6 * 24 * 60 * 60  # From MWD/t_U to J/t_U
E_fiss_MeV = 200   # Energy released per fission [MeV]
E_fiss_J = E_fiss_MeV * 1.60218e-13  # [J/fission]

# Impurity gases calculations (converted from ppm to absolute weight fraction --> 1 ppm = mg/Kg)
ppm_N2 = 25 * 1e-6  # [-]
ppm_H2O = 75 * 1e-6   # [-]

fission_yield = 0.28  # Fission yield for Xe and Kr [-]
f_R = 0.40  # Release rate of Xe and Kr from the fuel pellet [-]

# Molar masses [kg/mol]
# Let's assume that the enrichment of U-235 is 3%
enrichment = 0.03
M_U = (enrichment * 235.0 + (1 - enrichment) * 238.0) * 1e-3
M_O = 16.00e-3
M_UO2 = M_U + 2 * M_O
M_N2 = 28.02e-3
M_H2O = 18.016e-3

# CALCULATIONS OF THE TOTAL MOLES OF GAS IN PLENUM

pellet_density_percent = 0.955  # (porosity 4.5%)
dens_th_fuel = 10960.0  # UO2 theoretical density [kg/m^3]
dens_real_fuel = pellet_density_percent *  dens_th_fuel  # Actual density considering porosity

#tot_fuel_rods = num_fuel_assemblies * fuel_rods_per_assembly
fuel_area = np.pi * (pellet_diameter / 2) ** 2  # [m^2]
vol_fuel = fuel_area * core_height  # Fuel volume (UO2) [m^3]
m_fuel_uo2 = vol_fuel * dens_real_fuel  # Mass of UO2 per single element

#m_fuel_uo2 = fuel_weight / tot_fuel_rods # # Mass of UO2 per rod [kg]

m_uranio = m_fuel_uo2 * (M_U / M_UO2)  # Mass of Uranium per rod [kg]
tonn_u_rod = m_uranio / 1000.0  # Mass of Uranium per rod [tons]

total_energy_rod_J = Burnup_Jt * tonn_u_rod  # Total thermal energy per rod [J]

n_fissions_rod = total_energy_rod_J / E_fiss_J  # Total fissions per rod
n_rilasciate = n_fissions_rod *  fission_yield * f_R
N_A = 6.02214e23  # Avogadro's number [atoms/mol]

moles_fiss_gas= n_rilasciate / N_A  # Moles of Xe+Kr [mol]

m_N2_rod = m_fuel_uo2 * ppm_N2
m_H2O_rod = m_fuel_uo2 * ppm_H2O

moles_N2 = m_N2_rod / M_N2
moles_H2O = m_H2O_rod / M_H2O

total_moles_gas = moles_fiss_gas + moles_N2 + moles_H2O

print(f"UO2 Mass per rod: {m_fuel_uo2:.4f} kg")
print(f"U Mass per rod: {m_uranio:.4f} kg ({tonn_u_rod*1e6:.2f} g)")
print(f"Total fissions per rod: {n_fissions_rod:.4f} ")

print(f"Moles of fission gas (Xe+Kr) released: {moles_fiss_gas:.6f} mol")
print(f"Moles of impurity N2: {moles_N2:.6f} mol")
print(f"Moles of impurity H2O: {moles_H2O:.6f} mol")
print(f"Total moles of gas in plenum: {total_moles_gas:.6f} mol")

# --- GEOMETRIC SIZING ---
R = 8.3144  # Ideal gas constant [J/(mol*K)]
A_internal = np.pi * r_i**2  # Internal cross-sectional area of cladding [m^2]
p_plenum = p_x

#Extracting the temperature value at the first node (Inlet)
T_ci_plenum = T_ci_vec[-1]
T_plenum_K = T_ci_plenum + 273.15  # [K]

V_plenum = (total_moles_gas * R * T_plenum_K) / (p_plenum * 1e6)  # Minimum volume of the plenum [m^3]  (p_pl converted to Pa via 1e6 because p_pl is in MPa)
H_plenum = V_plenum / A_internal  # Minimum height of the plenum [m]
H = H_plenum + 15e-2  # we add 15cm to account for the spring [cm]

print(f"At local inner cladding temperature (from Assignment 2): {T_ci_plenum:.2f} °C / {T_plenum_K:.2f} K")
print(f"  Gas plenum volume: {V_plenum * 1e6:.2f} cm³")
print(f"  Minimum gas plenum height: {H_plenum * 100:.2f} cm")
print(f"  Gas plenum Height accounting for the spring: {H* 100:.2f} cm")


# ========================================================================================================
# TASK 6.1: CHECK THE IDEAL GAS ASSUMPTION FOR THE H2O VAPOR IN THE GAS PLENUM
# ========================================================================================================
print("\n --- IDEAL GAS ASSUMPTION CHECK FOR H2O ---")

# To understand whether water behaves like an ideal gas or is at risk of condensing (becoming a liquid), we calculate its ideal partial pressure in the plenum. Next, we use the CoolProp library to calculate the saturation pressure of the water at the plenum temperature.
# If the calculated partial pressure is less than the saturation pressure, the ideal gas hypothesis is acceptable

p_H2O_id_Pa = (moles_H2O * R * T_plenum_K) / V_plenum  # Ideal partial pressure OF H2O in the plenum [Pa]
p_H2O_id_bar = p_H2O_id_Pa / 1e5  # [bar]

print(f"Ideal H2O partial pressure in the plenum: {p_H2O_id_bar:.4f} bar / {p_H2O_id_Pa:.4f} Pa")

P_critica_H2O = 22.06e6  # [Pa]
T_critica_H2O = 647.1    # [K]  373.95 °C + 273.25 K

# Calculation of reduced properties
P_ridotta = p_H2O_id_Pa / P_critica_H2O
T_ridotta = T_plenum_K / T_critica_H2O

print(f"H2O reduced properties -> P_r = {P_ridotta:.3f}, T_r = {T_ridotta:.3f}")

try:
    p_sat_H2O_Pa = CP.PropsSI('P', 'T', T_plenum_K, 'Q', 0, 'Water')  # [Pa]
    p_sat_H2O_bar = p_sat_H2O_Pa / 1e5  # [bar]
    Z_water = CP.PropsSI('Z', 'T', T_plenum_K, 'P', p_H2O_id_Pa, 'Water')
    coolprop_success = True
except Exception as e:
    print(f"CoolProp failed to calculate p_sat (Temperature might be above critical point: 373.95 °C).")
    coolprop_success = False

# Hypothesis testing

tol = 0.10  # Tollerance of 10%

if coolprop_success:
    print(f"Water saturation pressure at {T_ci_plenum:.2f} °C: {p_sat_H2O_bar:.2f} bar")
    print(f"Compressibility factor from CoolProp (Z): {Z_water:.4f}")

    if p_H2O_id_Pa < p_sat_H2O_Pa and abs(1.0 - Z_water) <= tol:
        print("--> CHECK PASSED: p_H2O < p_sat and Z is within 10% error (water remains in vapor phase). The ideal gas assumption is VALID.")
    else:
        print("--> CHECK FAILED: p_H2O >= p_sat. Water would condense. The ideal gas assumption is NOT valid.")
else:
    # If the plenum temperature is higher than 373.95 °C (647.1 KV), the water is in the supercritical state.
    # Above the critical temperature there is no liquid-vapor phase transition, so it cannot condense.
    if T_plenum_K > T_critica_H2O:
        print(f"--> TEMPERATURE IS SUPERCRITICAL ({T_ci_plenum:.2f} °C > 373.95 °C).")
        print(
            "    Water is a supercritical fluid, it cannot condense because it is above its critical point. The ideal gas assumption is acceptable.")

print("\n" + "="*70)




# =============================================================================
# 1. FIGURE: HOOP STRESS ANALYSIS (SUBPLOTS: 1 ROW, 3 COLUMNS)
# =============================================================================
fig1, axs1 = plt.subplots(1, 3, figsize=(18, 8), sharey=True)
(ax1, ax2, ax3) = axs1

# 1. SUBPLOT: Mechanical Stresses (Lamé)
ax1.plot(sigma_h_mec_in_vec, z_axis, color='blue', linestyle='-', label='Mechanical Hoop (Inner Wall)', linewidth=2)
ax1.plot(sigma_h_mec_out_vec, z_axis, color='red', linestyle='-', label='Mechanical Hoop (Outer Wall)', linewidth=2)
ax1.set_title("1. Mechanical Hoop Stress (Lamé)", fontsize=11, fontweight='bold')
ax1.set_xlabel("Hoop Stress [MPa]", fontsize=10)
ax1.set_ylabel("Core Height [m]", fontsize=10)
ax1.grid(True, linestyle=':', alpha=0.6)
ax1.legend(loc='best', fontsize=9)

# 2. SUBPLOT: Thermal Stresses
ax2.plot(sigma_h_th_in_vec, z_axis, color='blue', linestyle='-', label='Thermal Hoop (Inner Wall)', linewidth=2)
ax2.plot(sigma_h_th_out_vec, z_axis, color='red', linestyle='-', label='Thermal Hoop (Outer Wall)', linewidth=2)
ax2.axvline(x=0, color='black', linestyle=':', linewidth=1) # Reference zero-line
ax2.set_title("2. Thermal Hoop Stress Profile", fontsize=11, fontweight='bold')
ax2.set_xlabel("Hoop Stress [MPa]", fontsize=10)
ax1.set_ylabel("Core Height [m]", fontsize=10)
ax2.grid(True, linestyle=':', alpha=0.6)
ax2.legend(loc='best', fontsize=9)

# 3. SUBPLOT: Total Combined Stresses (Primary + Secondary)
ax3.plot(sigma_h_tot_in_vec, z_axis, color='blue', linestyle='-', label='Total Stress Hoop (Inner Wall)', linewidth=2.5)
ax3.plot(sigma_h_tot_out_vec, z_axis, color='red', linestyle='-', label='Total Stress Hoop (Outer Wall)', linewidth=2.5)
ax3.set_title("3. Total Hoop Stresses (Prim + Sec)", fontsize=11, fontweight='bold')
ax3.set_xlabel("Hoop Stress [MPa]", fontsize=10)
ax3.set_ylabel("Core Height [m]", fontsize=10)
ax3.grid(True, linestyle=':', alpha=0.6)
ax3.legend(loc='best', fontsize=9)

plt.suptitle("Cladding hoop stress distribution profiles along core height", fontsize=14, fontweight='bold', y=0.98)
plt.tight_layout()


# =============================================================================
# 2. FIGURE: LONGITUDINAL STRESS  ANALYSIS (SUBPLOTS: 1 ROW, 2 COLUMNS)
# =============================================================================
fig3, lxs1 = plt.subplots(1, 3, figsize=(18, 8), sharey=True)
(lx1, lx2, lx3) = lxs1

# SUBPLOT 1: Pure Mechanical Longitudinal Stress (Lamé)
lx1.plot(sigma_l_mec_in_vec, z_axis, color='purple', linestyle='-', label='Mechanical Longitudinal (Inner wall)', linewidth=2)
lx1.plot(sigma_l_mec_out_vec, z_axis, color='orange', linestyle='-', label='Mechanical Longitudinal (Outer Wall)', linewidth=2)
lx1.set_title(r"1. Mechanical Longitudinal Stress ($\sigma_l$)", fontsize=11, fontweight='bold')
lx1.set_xlabel("Longitudinal Stress [MPa]", fontsize=10)
lx1.set_ylabel("Core Height [m]", fontsize=10)
lx1.grid(True, linestyle=':', alpha=0.6)
lx1.legend(loc='best', fontsize=9)

# 2. SUBPLOT: Thermal Stresses
lx2.plot(sigma_l_th_in_vec, z_axis, color='purple', linestyle='--', label='Thermal Longitudinal (Inner Wall)', linewidth=2)
lx2.plot(sigma_l_th_out_vec, z_axis, color='orange', linestyle='--', label='Thermal Longitudinal (Outer Wall)', linewidth=2)
lx2.axvline(x=0, color='black', linestyle=':', linewidth=1) # Reference zero-line
lx2.set_title("2. Thermal Longitudinal Stress Profile", fontsize=11, fontweight='bold')
lx2.set_xlabel("Longitudinal Stress [MPa]", fontsize=10)
lx1.set_ylabel("Core Height [m]", fontsize=10)
lx2.grid(True, linestyle=':', alpha=0.6)
lx2.legend(loc='best', fontsize=9)

# 3. SUBPLOT: Total Combined Stresses (Primary + Secondary)
lx3.plot(sigma_l_tot_in_vec, z_axis, color='purple', linestyle='-', label='Total Stress Longitudinal (Inner)', linewidth=2.5)
lx3.plot(sigma_l_tot_out_vec, z_axis, color='orange', linestyle='-', label='Total Stress Longitudinal (Outer)', linewidth=2.5)
lx3.axvline(x=0, color='black', linestyle=':', linewidth=1)
lx3.set_title(r"3. Total Longitudinal Stresses (Prim + Sec)", fontsize=11, fontweight='bold')
lx3.set_xlabel("Stress [MPa]", fontsize=10)
lx3.grid(True, linestyle=':', alpha=0.6)
lx3.legend(loc='best', fontsize=9)

plt.suptitle("Cladding longitudinal stress distribution profiles along core height", fontsize=14, fontweight='bold', y=0.98)
plt.tight_layout()
plt.show()

# =============================================================================
# 3. FIGURE: RADIAL STRESS ANALYSIS (SUBPLOTS: 1 ROW, 3 COLUMNS)
# =============================================================================
fig4, rxs1 = plt.subplots(1, 3, figsize=(18, 8), sharey=True)
(rx1, rx2, rx3) = rxs1  # Spacchettamento corretto dall'asse corrente rxs1

# SUBPLOT 1: Pure Mechanical Radial Stress (Lamé)
rx1.plot(sigma_r_mec_in_vec, z_axis, color='pink', linestyle='-', label='Mechanical Radial (Inner Wall)', linewidth=2)
rx1.plot(sigma_r_mec_out_vec, z_axis, color='limegreen', linestyle='-', label='Mechanical Radial (Outer Wall)', linewidth=2)
rx1.set_title(r"1. Mechanical Radial Stress ($\sigma_r$)", fontsize=11, fontweight='bold')
rx1.set_xlabel("Radial Stress [MPa]", fontsize=10)
rx1.set_ylabel("Core Height [m]", fontsize=10)
rx1.grid(True, linestyle=':', alpha=0.6)
rx1.legend(loc='best', fontsize=9)

# SUBPLOT 2: Thermal Stresses
rx2.plot(sigma_r_th_in_vec, z_axis, color='pink', linestyle='--', label='Thermal Radial (Inner Wall)', linewidth=2)
rx2.plot(sigma_r_th_out_vec, z_axis, color='limegreen', linestyle='--', label='Thermal Radial (Outer Wall)', linewidth=2)
rx2.axvline(x=0, color='black', linestyle=':', linewidth=1) # Reference zero-line
rx2.set_title("2. Thermal Radial Stress Profile", fontsize=11, fontweight='bold')
rx2.set_xlabel("Radial Stress [MPa]", fontsize=10)
rx2.grid(True, linestyle=':', alpha=0.6)
rx2.legend(loc='best', fontsize=9)

# SUBPLOT 3: Total Combined Radial Stresses (Primary + Secondary)
rx3.plot(sigma_r_tot_in_vec, z_axis, color='pink', linestyle='-', label='Total Stress Radial (Inner Wall)', linewidth=2.5)
rx3.plot(sigma_r_tot_out_vec, z_axis, color='limegreen', linestyle='-', label='Total Stress Radial (Outer Wall)', linewidth=2.5)
rx3.axvline(x=0, color='black', linestyle=':', linewidth=1) # Reference zero-line
rx3.set_title(r"3. Total Radial Stresses (Prim + Sec)", fontsize=11, fontweight='bold')
rx3.set_xlabel("Radial Stress [MPa]", fontsize=10)
rx3.grid(True, linestyle=':', alpha=0.6)
rx3.legend(loc='best', fontsize=9)

plt.suptitle("Cladding radial stress distribution profiles along core height", fontsize=14, fontweight='bold', y=0.98)
plt.tight_layout()

plt.show()


# =============================================================================
# 2. FIGURE: COMPREHENSIVE ASME SECTION III VERIFICATION (Sm & 3*Sm CONTROLS)
# =============================================================================
fig2, ax4 = plt.subplots(figsize=(9, 8.5))

# --- VERIFICA 1: PRIMARY STRESSES ---
ax4.axvline(x=tresca_pri, color='navy', linestyle='-', label=r'Primary Stress (Tresca)', linewidth=2)
ax4.axvline(x=S_m_MPa, color='teal', linestyle='--', label=r'ASME Limit for Primary Stress ($S_m$)', linewidth=2)

# --- VERIFICA 2: PRIMARY + SECONDARY STRESSES ---
ax4.plot(tresca_tot_in_vec, z_axis, color='crimson', linestyle='-', label=r'Total Tresca Stress (Inner Wall)', linewidth=2.5)
ax4.plot(tresca_tot_out_vec, z_axis, color='darkorange', linestyle='-', label=r'Total Tresca Stress (Outer Wall)', linewidth=2.5)
ax4.axvline(x=S_m_2_MPa, color='darkred', linestyle='-.', label=r'ASME Limit for Primary+Secondary ($3 \cdot S_m$)', linewidth=2.5)

ax4.set_title("ASME Section III Structural Verification", fontsize=12, fontweight='bold', pad=15)
ax4.set_xlabel("Equivalent Tresca Stress [MPa]", fontsize=11)
ax4.set_ylabel("Core Height [m]", fontsize=11)

max_stress_presente = max(S_m_2_MPa, np.max(tresca_tot_in_vec), np.max(tresca_tot_out_vec))
ax4.set_xlim(0, max_stress_presente * 1.15)

ax4.grid(True, linestyle=':', alpha=0.6)
ax4.legend(loc='lower right', fontsize=9, framealpha=0.95, edgecolor='gray')

plt.tight_layout()
plt.show()