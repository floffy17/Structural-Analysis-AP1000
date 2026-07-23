# AP1000 Fuel Pin Thermo-Mechanical Verification

![Politecnico di Torino](https://img.shields.io/badge/University-Politecnico_di_Torino-blue.svg)
![Field](https://img.shields.io/badge/Field-Nuclear_Engineering-green.svg)
![Standard](https://img.shields.io/badge/Standard-ASME_BPVC_Sec_III-red.svg)

---

## Overview
This repository contains the computational model, analytical formulations, and Python routines developed for **Assignment 3** of the **Nuclear Fission Plants** course at **Politecnico di Torino** (A.Y. 2025/2026).

The project performs a comprehensive thermo-mechanical structural verification of a fuel pin cladding and pellet for the **Westinghouse AP1000** Pressurized Water Reactor (PWR) operating under nominal steady-state conditions. Operating in extreme in-core environments requires conservative verification of elastic instability (buckling), internal pressure limits, and combined primary/secondary stress states in strict compliance with the **ASME Boiler and Pressure Vessel Code (BPVC) Section III**. In addition, a detailed gas plenum geometric sizing and thermodynamic validation of internal gas species are carried out.

## Key Objectives & Analysis Features
- **Preliminary Buckling Verification:** Evaluates the critical elastic instability pressure (p_cr) using Zircaloy-4 temperature-dependent properties evaluated at peak cladding temperatures to ensure resistance against primary system collapse (155.13 bar).
- **Thick-Walled Pressure Limits:** Accounts for thick-walled cylindrical geometry (t/r_avg = 0.13 > 0.10) via Lamé's equations and Tresca's yield criterion under a conservative bounding depressurization scenario (p_o = 0 bar).
- **Axial Mapping of Thermo-Mechanical Stresses:** Computes continuous axial vectors for hoop (σ_h), radial (σ_r), and longitudinal (σ_l) stress components along the active core height (4.2672 m), accounting for temperature-dependent material properties (E(z), ν(z), α(z)).
- **ASME BPVC Section III Compliance:** Verifies primary membrane stress intensities against allowable design limits (S_m) and combined Primary + Secondary stresses against cyclic plastic thresholds (3S_m).
- **Fission Gas & Plenum Sizing:** Quantifies noble fission gas generation (Xe + Kr at 60000 MWD/t_U burnup) and manufacturing impurities (N2, H2O vapor) to size the active gas plenum height (H_plenum) and total pin length including the hold-down spring (H_spring = 15 cm).
- **Thermodynamic Gas Phase Validation:** Dynamically checks the ideal gas assumption for H2O vapor at plenum temperature (346.10 °C) against CoolProp saturation pressure (p_sat,H2O) and compressibility factor limits (Z_water).

## Methodology & Standards
1. **Geometry & Operating Parameters:** AP1000 17x17 XL fuel rod design (D_ro = 9.50 mm, t_c = 0.5715 mm, Zircaloy-4 cladding, 3% enriched UO2 pellet).
2. **Stress Categorization (ASME Section III):**
   - *Primary Membrane Stress (P_m)*: Evaluated as wall-averaged mechanical stresses driven by internal pressure (p_i,max = 27.25 MPa) and external coolant pressure (p_o = 15.513 MPa).
   - *Secondary Stress (Q)*: Evaluated at inner (r_i) and outer (r_o) cladding boundaries, driven by nonlinear radial thermal gradients ΔT_cl(z) and localized mechanical deviations.
3. **Failure & Design Limits:** 
   - Primary stress limit: S_m = min(2/3 σ_y, 1/3 σ_u) = 137.67 MPa.
   - Primary + Secondary stress limit: 3S_m = 413.00 MPa.
   - Buckling pressure limit: p_cr > P_sys.

## Repository Structure
```text
.
├── src/                # Python calculation routines and CoolProp integration
├── data/               # AP1000 core parameters, material constants, and ASME limits
├── docs/               # Technical documentation and Assignment formulation
├── results/            # Axial stress profile plots and ASME verification figures
└── README.md           # Project documentation
