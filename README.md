# Thyrosim Architecture

## thyrosim.jl
---
Collects the ODE's which comprise the thyrosim model with the initial conditions  

## thyrosim_odes.jl
---
System of ordinary differential equations and initial conditions comprosing the model. To be modified with D1 and D2 (Graves model) and Immune subystem (Graves and Hahshimotos models). Details of the model can be found [here](https://www.frontiersin.org/articles/10.3389/fendo.2022.888429). 

## utilities.jl
---
* reads data from blakesley, Jonklås (2 sets), Schneider**
* provides functions for plotting the results of Thyrosim.jl, optionally overlaying the Blakesley, Jonklås, or Schneider data**
* contains the simulation function
* initializes the IC variables and parameter values
* sets adjustable parameters (i.e. T3/T4 infusion doising, T3/T4 secretion and excretion rates)
* scales plasma (Feldchush), slow, fast compartment parameters based on sex, height, weight

## Key functions
---
### ---Simulate---
    Note: some conditions seem to have identical default values in both the simulate funciton (utilities.jl) and in the initialize funciton (thyrosim_odes.jl)... not necessary?

*Inputs:* 
- h: height in meters
- w: weight in kg
- sex: true --> male, false --> female
- days: default to 30 (must be integer)
- dial: [T4 secretion; T4 excretion; T3 secretion; T3 excretion]
    - *defaults to dial=[1.0; 0.88; 1.0; 0.88], unclear units....*
- T4dose: micrograms
- T3dose: micrograms
- dosing_interval: time between doses, in hours
- warmup (*default=true*): runs model assuming healthy thyroid function for 30 days to approximate ICs
- fitting_index: ???
- parameters: ???
- fixed_parameters: ???

### ---Thyrosim---
dq[1]: $\frac{dT4}{dt}$

dq[4]: $\frac{dT3}{dt}$

dq[7]: $\frac{dTSHp}{dt}$

Corresponding T4, T3, and TSHp are in the same number column in the solution vector. Other columns correspond to T3/T4 in fast, slow compartemnts, the gut, exogenous input of T3/T4, and delay fucntions. Thyrosim uses the built in `solve` to compute the solution vector. We can grab a lot of useful information from the solution (sol), for example:

`sol.t` - array of time points solution was saved at\
`sol.u` - array of solution values\
`sol.u[i]` - value of the solution at `sol.t[i]`\
`sol.prob` - returns the ODE which was solved. Thyrosim has defaults `timespan: (0.0, 720.0)` and 19 initial values contained in `u0`. 

More documentation can be found [here](https://docs.sciml.ai/DiffEqDocs/stable/).

## Changelog
1. (not yet done) Added placeholder equations for D1 and D2
    1.1 D1 ... dq[14] and D2 ... dq[15] in thyrosim_ODES.jl <br>
    1.2 Renumbered dq[14...19] to dq[16...21] across <br>
    1.3 ...

## TODO
---
* Collect Fitting data from Ben: Blaksley, Jonklås, Jonklås new, Schneider? Need the data to use some of the plotting functions given in utilities.jl... talk with Ben
* Update with D1, D2. Should go in thyrosim_odes.jl block starting at 456.
* Separate old thyrosim Model from current
* Implement items in changelog


## Notes
1. To use data, added folder "data" at C:\Users\Aidan\.julia\packages\Thyrosim\4xPMT\data and inserted the 4 data folders from ben.
2. To display interactive plots in browser, run `plotly()` in REPL after initalizing Thyrosim.

## Questions
1. Method Error with `plot_jonklaas`. Also not sure if it is necessary to provide T4data, T3data, TSHdata speratley?
**Note: vectorizing leads to bounds error but allows function call**
