# ABOUT😀

this repository contains code of paper "Longitudinal vibration modeling of hybrid cable-pulley systems with time-varying cable length", the main structure is listed below:

- data: contains the data of this research
- doc: contains some notes about this repository
- src: contains the main code of this research
  - maplesim: maplesim files
  - others: julia code
- README.md
  - this file contains the main instructions about this repository

## REPOSITORY STRUCTURE🤖

### Data File Name Rules👀

- full/partial: full model or partial model data
- N3/5/10: cable segement N=3/5/10
- T50/100/200: cable tenstion is 50/100/200N, represent the torque T=0.75/1.5/3Nm
- ideal: indicates the data in numerical simulation(with ideal) or physical experiment(no symbol)
- 1p/3p: indicates the pulley number, 1 pulley or 3 pulleys

for example, the `full_N3_T50_ideal_1p.mat` represents the data is in experiment with conditions:

- full model
- cable segement N=3
- winch torque T=0.75Nm
- is in numerical simulation
- with 1 pulley

the `partial_N5_3p.mat` represents the data is in experiment with conditions:

- partial model
- cable segement N=5
- winch torque T=0.726Nm(excited by mass gravity)
- is in physical simulation
- with 3 pulley

### CODE FILE STRUCTURES👻

**julia**

- `bc`: contains the initial conditions of the ODE
- `main`: contains the main solve codes, also provide data export and figure illustration generation
- `ode_full`: the odes of the full model
- `ode_partial`: the odes of the partial model
- `param_1p`: the parameters of 1 pulley case of physical experiment
- `param_3p`: the parameters of 3 pulley case of physical experiment
- `param_ideal_1p`: the parameters of 1 pulley case of numberical simulation
- `param_ideal_3p`: the parameters of 3 pulley case of numberical simulation
- `post`: generate data from solver results with given sample frequency(20kHz)
- `util_full`: util functions of full model
- `util_partial`: util functions of partial model

**maplesim**

- `WinchPulleyModel_1p`: 1 pulley file of physical experiment
- `WinchPulleyModel_3p`: 3 pulley file of physical experiment
- `WinchPulleyModel_ideal_1p`: 1 pulley file of numerical simulation
- `WinchPulleyModel_ideal_3p`: 3 pulley file of numerical simulation

## HOW TO USE🤠

### FULL MODEL AND PARTIAL MODEL(Julia code)

execute the `main.jl` to simulate the model.

> Because julia has JIT, the first run could be extremely slow due to the compile of packages and codes

> The suggest way to run julia code is to execute it in VSCode, the manuals can be found in:
>
> - [julia official site](https://discourse.julialang.org/)
> - [configuration of julia in VSCode](https://code.visualstudio.com/docs/languages/julia)

- change the `model_type` in `main.jl` from `full` to `partial` to switch the model as full model to partial model
- change the `param_type` to switch pulley types
- change the `N` to change cable segements N
- change the `tspan` to change the simulation time duration
- change the `tol` to change the solver tolerance
- comment on lines 91/92 to switch between single simulation or benchmark simulation

### RIGID MODEL(Maplesim code)

The configuration can be changed in maplesim the property panel. Cick the triangular `run` button to simulate the mode.

![maplesim configuration](asset/maplesim_panel.png)
