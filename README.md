# Karim and Hudson

Once Joe gets you guys a good dataset, you will have to do the following. 

1. Successfuly solve and plot ODE with some set of parameters and ICs
2. Modify the solver to accept dose input at specified timepoints during solution to emulate patient taking LT4
3. Come up with loss function and parameter estimation method to fit the model to dataset
4. Sensitivity analysis on parameters/parameter ranking

Can try in whatever language works best, you can see my attemps under the `thyrosimIM/src/notebooks` folder but the code is a mess and I don't really remember what everything does so would reccomend starting from scratch with the ODEs and parameters/ICs provided below.

### ODEs


Here are the up to date ODEs for the model formatted for Julia. The `q` vector holds state variables, `p` holds paramters. `dq` holds the differentials. I would refrain from changing anything in `dq[1:19]` as this represents the original thyrosim model which works well.

```
@muladd function thyrosimIM!(dq, q, p, t)
    kdelay = 5/8

    # scaling the mass/concentration of compartments
    plasma_volume_ratio = p[69]^p[71]
    slow_volume_ratio = p[74]^p[71]
    fast_volume_ratio = p[75]^p[71]

    # scale comparment sizes
    q1 = q[1] * 1 / p[69]
    q2 = q[2] * 1 / p[75]
    q3 = q[3] * 1 / p[74]
    q4 = q[4] * 1 / p[69]
    q5 = q[5] * 1 / p[75]
    q6 = q[6] * 1 / p[74]
    q7 = q[7] * 1 / p[69]

    # Auxillary equations
    q4F = (p[24]+ p[25] * q1 + p[26] * q1^2 + p[27] * q1^3) * q4 #FT3p
    q1F = (p[7] + p[8] * q1 + p[9] * q1^2 + p[10] * q1^3) * q1  #FT4p
    SR3 = (q[24]/p[100])*(p[19] * p[59] * q[19]) # Scaled (q[24]/p[99]) Brain delay (dial 3)
    SR4 = (q[24]/p[100])*(p[1] * p[57] * q[19])  # Scaled (q[24]/p[99]) Brain delay (dial 1)
    fCIRC = q[9]^p[51] / (q[9]^p[51] + p[49]^p[51])
    SRTSH = (p[30]+p[31]*fCIRC*sin(pi/12*t-p[33]))*(p[50]^p[52]/(p[50]^p[52] + q[9]^p[52]))
    fdegTSH = p[34] + p[35] / (p[36] + q7)
    fLAG = p[41] + 2*q[8]^11 / (p[42]^11 + q[8]^11)
    f4 = p[37]*(1 + 5*(p[53]^p[54]) / (p[53]^p[54]+q[8]^p[54]))
    NL = p[13] / (p[14] + q2)

    # ODEs
    dq[1]  = p[81] + (SR4 + p[3] * q2 + p[4] * q3 - (p[5] + p[6]) * q1F) * plasma_volume_ratio + p[11] * q[11] #T4dot (need to remove u1)
    dq[2]  = (p[6] * q1F - (p[3] + p[12] + NL) * q2) * fast_volume_ratio                                    #T4fast
    dq[3]  = (p[5] * q1F -(p[4] + p[15] / (p[16] + q3) + p[17] /(p[18] + q3)) * q3) * slow_volume_ratio  #T4slow
    dq[4]  = p[82] + (SR3 + p[20] * q5 + p[21] * q6 - (p[22] + p[23]) * q4F) * plasma_volume_ratio + p[28] * q[13] #T3pdot
    dq[5]  = (p[23] * q4F + NL * q2 - (p[20] + p[29]) * q5) * fast_volume_ratio                         #T3fast
    dq[6]  = (p[22] * q4F + p[15] * q3 / (p[16] + q3) + p[17] * q3 / (p[18] + q3) -(p[21])*q6) * slow_volume_ratio #T3slow
    dq[7]  = (SRTSH - fdegTSH * q7) * plasma_volume_ratio                                           #TSHp
    dq[8]  = f4 / p[38] * q1 + p[37] / p[39] * q4 - p[40] * q[8]          #T3B
    dq[9]  = fLAG * (q[8] - q[9])                                             #T3B LAG
    dq[10] = -p[43] * q[10]                                                   #T4PILLdot
    dq[11] =  p[43] * q[10] - (p[44] * p[58]+ p[11]) * q[11]                  #T4GUTdot: note p[44] * p[58] = p[44] * dial[2] = k4excrete
    dq[12] = -p[45] * q[12]                                                   #T3PILLdot
    dq[13] =  p[45] * q[12] - (p[46] * p[60] + p[28]) * q[13]                 #T3GUTdot: note p[46] * p[60] = p[46] * dial[4] = k3excrete

    # Delay ODEs -- why do we have so many chained together like this??
    dq[14] = kdelay * (q7 - q[14]) 
    dq[15] = kdelay * (q[14] - q[15])                                         #delay2: TSH delay
    dq[16] = kdelay * (q[15] - q[16])                                         #delay3
    dq[17] = kdelay * (q[16] - q[17])                                         #delay4
    dq[18] = kdelay * (q[17] - q[18])                                         #delay5
    dq[19] = kdelay * (q[18] - q[19])                                         #delay6

    # ---------- IMMUNE ODEs ---------- q[20:Bcell, 21:Pcell, 22:Tcell, 23:Cytokine, 24:FTS, 25:TPOAb]

    dq[20] = p[83]*(q[23]/(q[23]+p[95]))*q[22]+(p[97]*q1F)-(p[89]+p[84])*q[20] # Bdot -- changed to T4 stimulation
    dq[21] = p[84]*q[20]-p[90]*q[21] # Pdot
    dq[22] = p[85]*q[24]+p[98]*q1F-p[99]*(q[23]/(q[23]+p[96]))*q[22]-p[91]*q[22] # Tdot
    dq[23] = p[86]*q[22]-p[92]*q[23] # Cdot
    dq[24] = p[87]*((q7/q[24])*p[100])-p[93]*(q[24])*q[25] #p[87]*((q7/q[24])*p[100])-p[93]*(q[24]/p[100])*q[25] # FTSdot MODIFIED
    dq[25] = p[88]*q[21]-q[25]*(p[94]+p[93]*q[24]) # Abdot

    return nothing
end

```

### Parameters and ICs

Here are some initial conditions and parameter values to start with. 

```
# OG parameters
p[1] = 0.0027785399344 #S4
p[2] = 8               #tau
p[3] = 0.868           #k12
p[4] = 0.108           #k13
p[5] = 584             #k31free
p[6] = 1503            #k21free
p[7] = 0.000289        #A
p[8] = 0.000214        #B
p[9] = 0.000128        #C
p[10] = -8.83*10^-6    #D
p[11] = 0.88           #k4absorb
p[12] = 0.0189         #k02
p[13] = 0.012101809339 #VmaxD1fast
p[14] = 2.85           #KmD1fast
p[15] = 6.63*10^-4     #VmaxD1slow
p[16] = 95             #KmD1slow
p[17] = 0.00074619     #VmaxD2slow
p[18] = 0.075          #KmD2slow
p[19] = 3.3572*10^-4   #S3
p[20] = 5.37           #k45
p[21] = 0.0689         #k46
p[22] = 127            #k64free
p[23] = 2043           #k54free
p[24] = 0.00395        #a
p[25] = 0.00185        #b
p[26] = 0.00061        #c
p[27] = -0.000505      #d
p[28] = 0.88           #k3absorb
p[29] = 0.184972339613 #k05
p[30] = 450            #Bzero
p[31] = 219.7085301388 #Azero
p[32] = 0              #Amax
p[33] = -3.71          #phi
p[34] = 0.53           #kdegTSH-HYPO
p[35] = 0.226          #VmaxTSH
p[36] = 23             #K50TSH
p[37] = 0.058786935033 #k3
p[38] = 0.29           #T4P-EU
p[39] = 0.006          #T3P-EU
p[40] = 0.037          #KdegT3B
p[41] = 0.0034         #KLAG-HYPO
p[42] = 5              #KLAG
p[43] = 1.3            #k4dissolve
p[44] = 0.12           #k4excrete
p[45] = 1.78           #k3dissolve
p[46] = 0.12           #k3excrete
p[47] = 3.2            #Vp
p[48] = 5.2            #VTSH
p[49] = 3.001011022378 #K_circ
p[50] = 3.094711690204 #K_SR_tsh
p[51] = 5.674773816316 #n_hillcirc
p[52] = 6.290803221796 #m_hillTSH
p[53] = 8.498343729591 #K_f4 for f4
p[54] = 14.36664496926 #l_hillf3
p[57] = dial[1] # controls T4 secretion rate
p[58] = dial[2] # controls T4 excretion rate
p[59] = dial[3] # controls T3 secretion rate
p[60] = dial[4] # controls T3 excretion rate
p[61] = 5.003761571969437   # σT4
p[62] = 0.11122955089297369 # σT3
p[63] = 0.4                 # σTSH
p[64] = 0.1                 # σFT4
p[65] = 21.82854404275587 # maleBMI_ref
p[66] = 22.99050845201536 # femaleBMI_ref
p[67] = 1.0 #Vtsh_scale
p[69] = 1.0 # PV_ratio
p[70] = -1.0 # PV
p[71] = 1.0 # PV_allometric_exp
p[72] = 1.0 # fat_free
p[73] = 0.0 # fat
p[74] = 1.0 # slow_scale
p[75] = 1.0 # fast_scale
p[76] = 0.75 # male_allometric
p[77] = 0.75 # female_allometric
p[78] = 1.7608716659237555 # male_ref_height
p[79] = 1.6696106891941103 # female_ref_height
p[80] = 1.0499391485135692 # male_clearace
p[81] = 0.0 # T4 infusion
p[82] = 0.0 # T3 infusion
# Immune paramters
p[83] = 3e-4 # 83 B-cell activation rate, will probably be lower due to T3 term p[15]
p[84] = 1e-2 # 84 Plasma cell transformation rate
p[85] = 8.05e-1 # 85 CD4+ activation rate
p[86] = 51.84e5 # 86 Cytokine production rate
p[87] = 1e5 # 87 relative growth rate of FTS
p[88] = 2e2 # 88 combined antibody production rate
p[89] = 5e-1 # 89 B-cell death rate
p[90] = 1.0e-2 # 90 Plasma cell death rate
p[91] = 2.41e-1 # 91 CD4+ cell death rate
p[92] = 1.189e3  # 92 Cytokine degredation rate
p[93] = 1e-2 # 93 Functional thyroid destruction rate
p[94] = 5.74e1 # 94 Blood Ab degredation rate
p[95] = 18e5 # 95 B-cell cytokine binding activation threshold
p[96] = 2e6 # 96 CD4+ T-cell cytokine binding activation threshold
p[97] = 1e5 # 97 FT4 B-cell stimulation
p[98] = 1e3 # 98 FT4 T-cell stimulation
p[99] = 9.1e-4 # 99 CD4+ T-cell inhibition rate
p[100] = 13.5 # 100 Euthyroid FTS

# Fitting Variance
p[101] = 1 # 101 T4 Variance
p[102] = 1 # 102 T3 Variance
p[103] = 1 # 103 TSH Variance
p[104] = 1 # 104 Lymphocyte Variance
p[105] = 1 # 105 Antibody Variance

# Initial Conditions
ic[1] = 0.322114215761171 #T4dot
ic[2] = 0.201296960359917 #T4fast
ic[3] = 0.638967411907560 #T4slow
ic[4] = 0.00663104034826483 #T3pdot
ic[5] = 0.0112595761822961 #T3fast
ic[6] = 0.0652960640300348 #T3slow
ic[7] = 1.78829584764370 #TSHp
ic[8] = 7.05727560072869 #T3B
ic[9] = 7.05714474742141 #T3B_lag
ic[10] = 0 #T4PILLdot
ic[11] = 0 #T4GUTdot
ic[12] = 0 #T3PILLdot
ic[13] = 0 #T3GUTdot
ic[14] = 3.34289716182018 #delay1
ic[15] = 3.69277248068433 #delay2
ic[16] = 3.87942133769244 #delay3
ic[17] = 3.90061903207543 #delay4
ic[18] = 3.77875734283571 #delay5
ic[19] = 3.55364471589659 #delay6
#immune
ic[20] = 320 # B-cells
ic[21] = 80 # Plasma cells 
ic[22] = 680 # CD4+ cells
ic[23] = 5e9 # Cytokines
ic[24] = 5 # FTS
ic[25] = 240*10e-13 # Antibodies

```
You will almost definitley need to change some or all of `q[20:25]` and `p[83:end]` as you fit the model. Depending on the loss function you use for fitting and how you do parameter estimation, you can get rid of p[101:105] as these were for maximum liklihood criterion (see [past paper](https://doi.org/10.3389/fendo.2022.888429) from Joe's lab that uses this technique). You might also notice some parameters like `p[81]` for exogenous T3/T4 input. Can try to use these or come up with your own method for doses.

Feel free to reach out if you have any questions.


# _Outdated_: ThyrosimIM Workflow:

_All of new code is in `thyrosimIM/src/notebooks`. All of Ben's code is in `Thyrosim.jl/src`_. If possible, I would try starting from scratch and working through the steps below on your own, using my code and Ben's code as reference as needed. Trying to modify the existing code will probably lead to errors that are difficult to fix.

All documents (i.e. paper drafts, diagrams, etc.) are in `thyrosimIM/src/docs`. All data currently pulled from the python notebook are in `thyrosimIM/src/data`. Again, would reccomend grabbing data using the python notebook yourself so you know what is going on.

1.	**Import patient data to Julia notebook.** I would recommend using the DataFrames.jl library
    - Patient data should all be in terms of days. However the model is in terms of hours so make sure you adjust the timescale of the data either in Julia or in the original python notebook that pulls the data.
    - Plotting the data will help you know if you have adjusted correctly. Make sure to check both lab and medications
2.	**Parameters and initial conditions.** Ben initialized them in a very complicated way to allow parameter variation between different sex/BMI. For now, would recommend initializing them without patient specific variation and adding them back in later - makes fitting easier
3.	**ODEs**. System to be solved. Takes in derivatives, state variables, parameters, and time.
4.	**Loss function.** In order to fit with an optimizer, we try to minimize loss (model predicted – patient data). To do so, we have to come up with a loss function. I use log least-squares loss because we have state variables on different magnitude scales. If you normalize the model, you can use a non-log loss function.
a.	You must convert model output to FT3, FT4 when calculating loss (the data is for free hormones) Here are the correct functions for conversion. Do not try to modify them!

    ```
    """
    Converts TT4 (μmol) to FT4 (μg/dL) using Thyrosim's internal 4th order polynomial
    """
    function FT4(TT4::Float64)
        return (0.000289 + (0.000214 * TT4) + (0.000128 * TT4^2) - (8.83e-6 * TT4^3)) * TT4 * 24281.25 * .45  # Convert from TT4 to FT4 then from μmol FT4 to μg/dL FT4 via scaling factor 
    end

    """
    Converts TT3 (μmol) to FT3 (μmol) using Thyrosim's internal 4th order polynomial
    """
    function FT3(TT3::Float64)
        return (0.00395 + (0.00185 * TT3) + (0.000610 * TT3^2) - (-0.000505 * TT3^3)) * TT3 * 28304000.35 * 0.5
    end
    ```

5.	**Fitting Routine.** This, along with proper import of patient data, is most difficult. You have to pass free parameters (along with an initial guess for these free parameters), fixed parameters, initial conditions, and data to the objective function (loss function) within the scope of the optimizer. To do so, you must pass the free parameters each time the optimizer runs (as they will be updated by the fitting function) and the rest only once.  
    - You also will probably need something to bound the parameters. You can do this manually via making loss infinity when parameters go out of bound or use Optim.jl’s bounded models.
    - One of the tricky aspects is that you have to deal with patient medications using something called a callback function. These meds must be passed through the gut compartment by adding them to the 10th state variable u[10]. To keep it simple, I’d reccomend using a periodic callback with a period of 24 hours. You will have to adjust the meds data to fit this, but due to the half life of T4 and long time scale of the patient data it is ok to combine 2x daily doses into one larger dose or split up less frequent doses. 
        - Most patients change dose over the course of the data so we cannot pass dose to the function like Ben did via parameters. Instead, make sure the medication data is a global variable and then do a time-match on the date of the dose to find out how much needs to be added to u[10]
    - As far as fitting algorithm, Nelder-Mead is the way to go initially. For parameter ID we will use Newton’s (or Newton family i.e. BFGS) to get an approximate inverse hessian which itself approximates covariance matrix. However to do this we need a different loss function (see 7 below)
6.	**Plot loss/results.** To determine if things are working, easiest thing is to graph the fit over the data. Would recommend doing this throughout the fitting process to see how many iterations are actually needed for fitting. Plotting loss is good sanity check as well, but you will have to initialize loss as a global variable and save it at each step if you want to keep track of it during fitting as it is thrown out at each step normally.

![Fit example](/thyrosimIM/src/notebooks/animation/optimization_08.27.1.png)

7.	**Parameter identifiability.** Once fitting looks ok (does not have to be perfect by any means, just tracks with general trends in T3, T4, TSH and fairly close to datapoints) we need to do another round of fitting using these good parameters. With this round, we need to use maximum likelihood estimation and add variance parameters with a new loss function. Ben has a good example in his code, and there may be an old version in mine. This allows us to track variance and thus get covariance estimate by fitting the new loss function with a Newton or Newton family method. 



### Summary of issues to be solved

|     Component    |     Issues / Todo    |
|---|---|
|     Patient Data    |     Convert   everything to hours, either in OG python notebook before you import to julia   or in julia itself. Make sure med dosing frequency is standardized at 24   hours and data corresponds accordingly.    |
|     Parameters   and IC    |     I think my   solution for free-fixed parameters is much easier to work with than Ben’s   initialization process, but feel free to experiment.     |
|     ODEs    |     None as far   as I know. Be careful if you modify as it can be hard to check ODEs are   correct without re-writing everything and comparing parameter numbers.   Current model is not perfect so consider alternatives (i.e. Meha’s) though   they also have issues.    |
|     Loss   functions    |     No problem   with loss function itself, but you have to make sure everything is on the   same timescale and units (i.e. everything in hours and FT4/FT3 not total   hormones). The maximum likelihood function needs work.    |
|     Fitting   routine    |     Callback   needs to be validated. Feed everything into u[10] and do not modify the   medication data frame when you do so.           Boundary   conditions should be implemented in one way or another. Ben and I both bypass   loss function and return INF loss if any parameters go out of bounds, but   Optim.jl has functions to handle this more gracefully.           Be sure to   write the fitting function modularly, do not try to make everything   one function or it is impossible to see where the problem is. Also be sure to   make it easy to change free and fixed parameters, do not have hard-coded   parameters in any of the fitting functions.    |
|     Plotting    |     The plots   look pretty good. You might have to modify them based on your other   functions, though.    |
|     Parameter ID    |     Maximum   likelihood loss function needs to be revamped. Note that you need to run   several iteration of BFGS or Newton’s method even if your parameters are good   from Nelder-Mead already – you aren’t trying to refit the model, but it takes   time to adjust the hessian and get valid variability estimates. Some way to   visualize this (i.e. heatmap) is pretty nice for sanity check.     |