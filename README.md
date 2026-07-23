# AP1000 Fuel Pin Thermo-Mechanical Verification

![Politecnico di Torino](https://img.shields.io/badge/University-Politecnico_di_Torino-blue.svg)
![Field](https://img.shields.io/badge/Field-Nuclear_Engineering-green.svg)
![Standard](https://img.shields.io/badge/Standard-ASME_BPVC_Sec_III-red.svg)

---

## Overview
This repository contains the computational model, technical documentation, and calculation scripts developed for **Assignment 3** of the **Nuclear Fission Plants** course at **Politecnico di Torino** (A.Y. 2025/2026).

The project focuses on the thermo-mechanical verification and structural integrity assessment of a fuel pin cladding and fuel pellet for the **Westinghouse AP1000** Pressurized Water Reactor under steady-state nominal operating conditions. Given the extreme in-core operating environment, the study guarantees structural safety against mechanical loads, severe temperature gradients, and internal pressure buildup, fully adhering to the nuclear design criteria established by the **ASME Boiler and Pressure Vessel Code (BPVC) Section III**.

## Key Objectives & Analysis Features
- **Elastic Instability Verification:** Evaluates the cladding resistance against external coolant pressure during early-life operation and depressurization events to prevent structural collapse.
- **Thick-Walled Cylinder Assessment:** Treats the cladding geometry with thick-walled structural criteria to model stress distributions accurately across the wall thickness.
- **Axial Stress Mapping:** Computes and tracks the continuous axial distribution of tangential, radial, and longitudinal stress components along the active core height, incorporating temperature-dependent material behavior.
- **ASME Section III Structural Compliance:** Evaluates primary membrane stresses against design thresholds and verifies that combined primary and secondary stresses remain within safe limits against cyclic plastic failure.
- **Fission Gas Production & Plenum Sizing:** Estimates total noble gas generation at high fuel burnup along with initial fabrication impurities, providing geometric sizing for the upper gas plenum and hold-down spring[cite: 1].
- **Thermodynamic Moisture Validation:** Performs a thermodynamic phase check on trapped moisture inside the plenum volume to verify that internal gases remain entirely in the vapor phase without risk of condensation[cite: 1].

## Methodology & Standards
1. **Operating Context & Core Geometry:** Based on standard AP1000 fuel assembly specifications, Zircaloy-4 cladding properties, and enriched uranium dioxide pellets[cite: 1].
2. **Stress Categorization:** Separates primary mechanical loads caused by system pressures from self-limiting secondary stresses driven by radial thermal gradients across the cladding[cite: 1].
3. **Regulatory Verification:** Compares the calculated stress states against allowable nuclear design intensity limits to ensure wide safety margins during steady-state operations[cite: 1].
