# CO₂ Methanation Kinetic Modeling

**Course:** ChE 401 – Senior Chemical Engineering Laboratory  
**Institution:** University of Massachusetts Amherst  
**Date:** Fall 2025  
**Note:** Experimental data collected as part of a team lab. 
Background and theoretical analysis section of the lab report written 
independently. Reaction order determination, rate constant calculation, 
and Arrhenius analysis — including all MATLAB code — conducted 
independently.

---

## Overview

Kinetic analysis of CO₂ methanation via the Sabatier reaction in a 
300 mL CSTR with 10 wt.% Ni/Al₂O₃ catalyst. Reaction order, rate 
constants, and activation energy were determined by fitting experimental 
concentration-time data to nonlinear kinetic models in MATLAB.

---

## Computational Approach

Seven reaction order models (-1, 0, 0.5, 1, 1.5, 2, 3) were developed 
and fit to experimental CA/CA₀ vs CA₀ data using least-squares 
optimization. Explicit models were fit using `lsqcurvefit`; implicit 
models (1.5 and 3rd order) were solved iteratively using `fsolve` and 
`fminsearch`. Model selection was based on SSE and RMSE comparison across 
all reaction orders at three temperatures.

---

## Key Results

- First order reaction confirmed across 250–300°C temperature range
- Rate constants: k = 0.301, 0.474, 0.629 s⁻¹ at 250, 270, 300°C
- Activation energy: 25.81 kJ/mol (250–270°C) and 37.57 kJ/mol 
(270–300°C), indicating transition from kinetic to mass transfer control
- Reactor volume: 305.38 ± 24.45 mL, consistent with manufacturer 
specification of 300 mL

---

## Files

- `rxn_order.m` — model fitting at 250°C
- `rxn_order_270_fixed.m` — model fitting at 270°C
- `rxn_order_300.m` — model fitting at 300°C
- `rxn_order_eqns.jpeg` — kinetic model equations reference

---

## Dependencies

MATLAB with Optimization Toolbox (`lsqcurvefit`, `fminsearch`, `fsolve`)
