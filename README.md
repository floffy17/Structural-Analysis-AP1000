# AP1000 Fuel Pin Thermo-Mechanical Verification

![Politecnico di Torino](https://img.shields.io/badge/University-Politecnico_di_Torino-blue.svg)
![Field](https://img.shields.io/badge/Field-Nuclear_Engineering-green.svg)
![Standard](https://img.shields.io/badge/Standard-ASME_BPVC_Sec_III-red.svg)

[English](#english) | [Italiano](#italiano)

---

<a name="english"></a>
## 🇬🇧 English

### Overview
This repository contains the computational model, analytical formulations, and Python routines developed for **Assignment 3** of the **Nuclear Fission Plants** course at **Politecnico di Torino** (A.Y. 2025/2026).

The project performs a comprehensive thermo-mechanical structural verification of a fuel pin cladding and pellet for the **Westinghouse AP1000** Pressurized Water Reactor (PWR) operating under nominal steady-state conditions. Operating in extreme in-core environments requires conservative verification of elastic instability (buckling), internal pressure limits, and combined primary/secondary stress states in strict compliance with the **ASME Boiler and Pressure Vessel Code (BPVC) Section III**. In addition, a detailed gas plenum geometric sizing and thermodynamic validation of internal gas species are carried out.

### Key Objectives & Analysis Features
- **Preliminary Buckling Verification:** Evaluates the critical elastic instability pressure ($p_{cr}$) using Zircaloy-4 temperature-dependent properties evaluated at peak cladding temperatures to ensure resistance against primary system collapse ($155.13\text{ bar}$).
- **Thick-Walled Pressure Limits:** Accounts for thick-walled cylindrical geometry ($t/r_{avg} = 0.13 > 0.10$) via Lamé's equations and Tresca's yield criterion under a conservative bounding depressurization scenario ($p_o = 0\text{ bar}$).
- **Axial Mapping of Thermo-Mechanical Stresses:** Computes continuous axial vectors for hoop ($\sigma_h$), radial ($\sigma_r$), and longitudinal ($\sigma_l$) stress components along the active core height ($4.2672\text{ m}$), accounting for temperature-dependent material properties ($E(z)$, $\nu(z)$, $\alpha(z)$).
- **ASME BPVC Section III Compliance:** Verifies primary membrane stress intensities against allowable design limits ($S_m$) and combined Primary + Secondary stresses against cyclic plastic thresholds ($3S_m$).
- **Fission Gas & Plenum Sizing:** Quantifies noble fission gas generation ($\text{Xe}+\text{Kr}$ at $60000\text{ MWD/t}_U$ burnup) and manufacturing impurities ($\text{N}_2$, $\text{H}_2\text{O}$ vapor) to size the active gas plenum height ($H_{plenum}$) and total pin length including the hold-down spring ($H_{spring} = 15\text{ cm}$).
- **Thermodynamic Gas Phase Validation:** Dynamically checks the ideal gas assumption for $\text{H}_2\text{O}$ vapor at plenum temperature ($346.10^\circ\text{C}$) against CoolProp saturation pressure ($p_{sat,H_2O}$) and compressibility factor limits ($Z_{water}$).

### Methodology & Standards
1. **Geometry & Operating Parameters:** AP1000 $17 \times 17$ XL fuel rod design ($D_{ro} = 9.50\text{ mm}$, $t_c = 0.5715\text{ mm}$, Zircaloy-4 cladding, $3\%$ enriched $\text{UO}_2$ pellet).
2. **Stress Categorization (ASME Section III):**
   - *Primary Membrane Stress ($P_m$)*: Evaluated as wall-averaged mechanical stresses driven by internal pressure ($p_{i,max} = 27.25\text{ MPa}$) and external coolant pressure ($p_o = 15.513\text{ MPa}$).
   - *Secondary Stress ($Q$)*: Evaluated at inner ($r_i$) and outer ($r_o$) cladding boundaries, driven by nonlinear radial thermal gradients $\Delta T_{cl}(z)$ and localized mechanical deviations.
3. **Failure & Design Limits:** 
   - Primary stress limit: $S_m = \min(\frac{2}{3}\sigma_y, \frac{1}{3}\sigma_u) = 137.67\text{ MPa}$.
   - Primary + Secondary stress limit: $3S_m = 413.00\text{ MPa}$.
   - Buckling pressure limit: $p_{cr} > P_{sys}$.

### Repository Structure
```text
.
├── src/                # Python calculation routines and CoolProp integration
├── data/               # AP1000 core parameters, material constants, and ASME limits
├── docs/               # Technical documentation and Assignment formulation
├── results/            # Axial stress profile plots and ASME verification figures
└── README.md           # Project documentation
