# Last updated from dev: May 08, 11:58

using DifferentialEquations, Optim, Plots, DataFrames, LinearAlgebra, CSV

# Placeholder until some way to pass fixed, free parameters separatley is implemented
function fixed_parameters()
    p = zeros(Float64, 100)
    dial=[1.0; 0.88; 1.0; 0.88]
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

    return p
end

function ics()
    ic    = zeros(Float64, 25)
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
    ic[20] = 100 # B-cells
    ic[21] = 20 # Plasma cells 
    ic[22] = 100 # CD4+ cells
    ic[23] = 5e9 # Cytokines
    ic[24] = 5 # FTS
    ic[25] = 2e9 # Antibodies

    return ic
end

# Placeholder until some way to pass fixed, free parameters separatley is implemented
function initialize_free()
    p = zeros(17)
    p[1] = 3e-3 # B-cell activation rate, will probably be lower due to T3 term p[15]
    p[2] = 1e-2 # Plasma cell transformation rate
    p[3] = 8.05e-1 # CD4+ activation rate
    p[4] = 51.84e5 # Cytokine production rate
    p[5] = 1e6 # relative growth rate of FTS
    p[6] = 1e6 # combined antibody production rate
    p[7] = 2e-6 # B-cell death rate
    p[8] = 4.0e-2 # Plasma cell death rate
    p[9] = 8.91e-3 # CD4+ cell death rate
    p[10] = .189  # Cytokine degredation rate
    p[11] = 1e-2 # Functional thyroid destruction rate
    p[12] = 1.74e-3 # Blood Ab degredation rate
    p[13] = 18e5 # B-cell cytokine binding activation threshold
    p[14] = 2e6 # CD4+ T-cell cytokine binding activation threshold
    p[15] = 1e3 # NOTE: NEED TO FIT and CHANGE
    p[16] = 9.1e-4 # CD4+ T-cell stimulation rate
    p[17] = 13.5 # Euthyroid FTS
    return p
end

# Preturbed free variables to see if fitting works
function initialize_free_varied()
    p = zeros(17)
    p[1] = 3.9e-3 # B-cell activation rate, will probably be lower due to T3 term p[15]
    p[2] = 1.2e-2 # Plasma cell transformation rate
    p[3] = 3.05e-1 # CD4+ activation rate
    p[4] = 56.84e5 # Cytokine production rate
    p[5] = 1.1e6 # relative growth rate of FTS
    p[6] = 1e6 # combined antibody production rate
    p[7] = 2.2e-6 # B-cell death rate
    p[8] = 3.0e-2 # Plasma cell death rate
    p[9] = 1.91e-4 # CD4+ cell death rate
    p[10] = .389  # Cytokine degredation rate
    p[11] = 1.4e-2 # Functional thyroid destruction rate
    p[12] = 1.34e-3 # Blood Ab degredation rate
    p[13] = 78e4 # B-cell cytokine binding activation threshold
    p[14] = 23e5 # CD4+ T-cell cytokine binding activation threshold
    p[15] = 3.2e3 # NOTE: NEED TO FIT and CHANGE
    p[16] = 9.1e-4 # CD4+ T-cell stimulation rate
    p[17] = 12.5 # Euthyroid FTS
    return p
end

function thyrosimIM_estimate(dq, q, p, t)
    kdelay = 5/8
    fixed_p = fixed_parameters()
    dial=[1.0; 0.88; 1.0; 0.88] # can probably remove this?

    # scaling the mass/concentration of compartments
    fixed_plasma_volume_ratio = fixed_p[69]^fixed_p[71]
    slow_volume_ratio = fixed_p[74]^fixed_p[71]
    fast_volume_ratio = fixed_p[75]^fixed_p[71]

    # scale comparment sizes
    q1 = q[1] * 1 / fixed_p[69]
    q2 = q[2] * 1 / fixed_p[75]
    q3 = q[3] * 1 / fixed_p[74]
    q4 = q[4] * 1 / fixed_p[69]
    q5 = q[5] * 1 / fixed_p[75]
    q6 = q[6] * 1 / fixed_p[74]
    q7 = q[7] * 1 / fixed_p[69]

    # Auxillary equations
    q4F = (fixed_p[24]+ fixed_p[25] * q1 + fixed_p[26] * q1^2 + fixed_p[27] * q1^3) * q4 #FT3p
    q1F = (fixed_p[7] + fixed_p[8] * q1 + fixed_p[9] * q1^2 + fixed_p[10] * q1^3) * q1  #FT4p
    SR3 = (q[24]/p[17])*(fixed_p[19] * fixed_p[59] * q[19]) # Scaled (q[24]/fixed_p[100]) Brain delay (dial 3)
    SR4 = (q[24]/p[17])*(fixed_p[1] * fixed_p[57] * q[19])  # Scaled (q[24]/fixed_p[100]) Brain delay (dial 1)
    fCIRC = q[9]^fixed_p[51] / (q[9]^fixed_p[51] + fixed_p[49]^fixed_p[51])
    SRTSH = (fixed_p[30]+fixed_p[31]*fCIRC*sin(pi/12*t-fixed_p[33]))*(fixed_p[50]^fixed_p[52]/(fixed_p[50]^fixed_p[52] + q[9]^fixed_p[52]))
    fdegTSH = fixed_p[34] + fixed_p[35] / (fixed_p[36] + q7)
    fLAG = fixed_p[41] + 2*q[8]^11 / (fixed_p[42]^11 + q[8]^11)
    f4 = fixed_p[37]*(1 + 5*(fixed_p[53]^fixed_p[54]) / (fixed_p[53]^fixed_p[54]+q[8]^fixed_p[54]))
    NL = fixed_p[13] / (fixed_p[14] + q2)

    # ODEs
    dq[1]  = fixed_p[81] + (SR4 + fixed_p[3] * q2 + fixed_p[4] * q3 - (fixed_p[5] + fixed_p[6]) * q1F) * fixed_plasma_volume_ratio + fixed_p[11] * q[11] #T4dot (need to remove u1)
    dq[2]  = (fixed_p[6] * q1F - (fixed_p[3] + fixed_p[12] + NL) * q2) * fast_volume_ratio                                    #T4fast
    dq[3]  = (fixed_p[5] * q1F -(fixed_p[4] + fixed_p[15] / (fixed_p[16] + q3) + fixed_p[17] /(fixed_p[18] + q3)) * q3) * slow_volume_ratio  #T4slow
    dq[4]  = fixed_p[82] + (SR3 + fixed_p[20] * q5 + fixed_p[21] * q6 - (fixed_p[22] + fixed_p[23]) * q4F) * fixed_plasma_volume_ratio + fixed_p[28] * q[13] #T3pdot
    dq[5]  = (fixed_p[23] * q4F + NL * q2 - (fixed_p[20] + fixed_p[29]) * q5) * fast_volume_ratio                         #T3fast
    dq[6]  = (fixed_p[22] * q4F + fixed_p[15] * q3 / (fixed_p[16] + q3) + fixed_p[17] * q3 / (fixed_p[18] + q3) -(fixed_p[21])*q6) * slow_volume_ratio #T3slow
    dq[7]  = (SRTSH - fdegTSH * q7) * fixed_plasma_volume_ratio                                           #TSHfixed_p
    dq[8]  = f4 / fixed_p[38] * q1 + fixed_p[37] / fixed_p[39] * q4 - fixed_p[40] * q[8]          #T3B
    dq[9]  = fLAG * (q[8] - q[9])                                             #T3B LAG
    dq[10] = -fixed_p[43] * q[10]                                                   #T4PILLdot
    dq[11] =  fixed_p[43] * q[10] - (fixed_p[44] * fixed_p[58]+ fixed_p[11]) * q[11]                  #T4GUTdot: note fixed_p[44] * fixed_p[58] = fixed_p[44] * dial[2] = k4excrete
    dq[12] = -fixed_p[45] * q[12]                                                   #T3PILLdot
    dq[13] =  fixed_p[45] * q[12] - (fixed_p[46] * fixed_p[60] + fixed_p[28]) * q[13]                 #T3GUTdot: note fixed_p[46] * fixed_p[60] = fixed_p[46] * dial[4] = k3excrete

    # Delay ODEs
    dq[14] = kdelay * (q7 - q[14]) 
    dq[15] = kdelay * (q[14] - q[15])                                         #delay2: TSH delay
    dq[16] = kdelay * (q[15] - q[16])                                         #delay3
    dq[17] = kdelay * (q[16] - q[17])                                         #delay4
    dq[18] = kdelay * (q[17] - q[18])                                         #delay5
    dq[19] = kdelay * (q[18] - q[19])
    
    # Immune ODEs
    dq[20] = p[1]*(q[23]/(q[23]+p[13]))*q[22]+p[15]*q4F-(p[7]+p[2])*q[20] # Bdot
    dq[21] = p[2]*q[20]-p[8]*q[21] # Pdot
    dq[22] = p[3]*q[24]+p[16]*(q[23]/(q[23]+p[14]))*q[22]-p[9]*q[22] # Tdot
    dq[23] = p[4]*q[22]-p[10]*q[23] # Cdot
    dq[24] = p[5]*((q7/q[24])*p[17])-p[11]*(q[24])*q[25] #p[9]*((q7/q[24])*p[17])-p[10]*(q[24]/p[17])*q[25] # FTSdot MODIFIED
    dq[25] = p[6]*q[21]-q[25]*(p[12]+p[11]*q[24]) # Abdot

    return dq
end

"""
Run parameter estimation on immune parameters. 
"""
function fit_params(data::DataFrame, p0::Vector, lb::Vector, ub::Vector)

    t = data.t
    ic = ics()    
    sol_index = [1,4,7,20,21,22,23,24,25]
    cols = names(data)[2:10]

    function error(p,t,data)
        tspan = (t[1],t[end])
        prob = ODEProblem(thyrosimIM_estimate, ic, tspan, p)
        sol = solve(prob, Rosenbrock23())

        loss = 0
        for i in 1:size(cols)[1] # can probably speed this up with eachcol notation
            dp = 1
            for time in t
                loss += ((sol(time)[sol_index[i]]-data[!,cols[i]][dp])/mean(data[!,cols[i]]))^2
                dp += 1
            end
        end
        return loss
    end

    function set_lower_bounds(p::Vector, lb::Vector)
        for i in eachindex(p)
            if p[i] < lb[i]
                p[i] = lb[i]
            end
        end
        return p
    end

    function set_upper_bounds(p::Vector, lb::Vector)
        for i in eachindex(p)
            if p[i] >= ub[i]
                p[i] = ub[i]
            end
        end
        return p
    end

    function objective(p_free, lb, ub)
        p_free = set_lower_bounds(p_free, lb)
        p_free = set_upper_bounds(p_free, ub)
        return error(p_free, t, data)
    end

    # Adjust so fitting indicies and parameters can be specified
    result = optimize(p -> objective(p, lb, ub), p0, NelderMead(),
     Optim.Options(time_limit = 300.0, iterations = 10000, g_tol=1e-5,
     show_trace = false, allow_f_increases=true))
end

"""
Calculate error prior to parameter fitting. 
"""
function prefit_error(data::DataFrame, p0::Vector)

    t = data.t
    ic = ics()    
    sol_index = [1,4,7,20,21,22,23,24,25]
    cols = names(data)[2:10]

    function error(p,t,data)
        tspan = (t[1],t[end])
        prob = ODEProblem(thyrosimIM_estimate, ic, tspan, p)
        sol = solve(prob, Rosenbrock23())

        loss = 0
        for i in 1:size(cols)[1] # can probably speed this up with eachcol notation
            dp = 1
            for time in t
                loss += ((sol(time)[sol_index[i]]-data[!,cols[i]][dp])/mean(data[!,cols[i]]))^2
                dp += 1
            end
            print("Column: ", cols[i], "\n")
            print("Current Error: ", loss, "\n\n")
        end
        return loss
    end

    function objective(p_free)
        return error(p_free, t, data)
    end

    objective(p0)

end

" Outputs a scatterplot of the data overlayed with model predictions pre and post-fitting"
function output_plotIM(sol_init, sol_fit, data; title::AbstractString = "ThyrosimIM Parameter Fitting", automargins::Bool=true)

    # can definitley clean up this code quite a bit later (looped plotting, get rid of unnecessary plots...)
    # parameters to adjust figure limits
    time = sol_init.t
    time_fit = sol_fit.t
    tvals = data.t./24

    t4lim, t3lim, tshlim = 140, 4, 10
    T4 = 777.0 * sol_init[1, :] / 3.2; T4_fit = 777.0 * sol_fit[1, :] / 3.2 
    T3 = 651.0 * sol_init[4, :] / 3.2; T3_fit = 651.0 * sol_fit[4, :] / 3.2
    TSH = 5.6 * sol_init[7, :] / 5.2; TSH_fit = 5.6 * sol_fit[7, :] / 5.2
    Bcell = sol_init[20, :]; Bcell_fit = sol_fit[20, :]
    PlasmaCell = sol_init[21, :]; PlasmaCell_fit = sol_fit[21, :]
    CD4Cell = sol_init[22, :]; CD4Cell_fit = sol_fit[22, :]
    Cytokines = sol_init[23, :]/6.022e11; Cytokines_fit = sol_fit[23, :]/6.022e11 # convert to nM from molecules/mL
    FTS = sol_init[24, :]; FTS_fit = sol_fit[24, :]
    TPOAb = sol_init[25, :]/6.022e8; TPOAb_fit = sol_fit[25, :]/6.022e8 # convert to pM from molecules/mL

    xlim=(0,time[end]/24) 
    if automargins
        t4lim = max(1.2maximum(T4), 130.0)
        t3lim = max(1.2maximum(T3), 2.5)
        tshlim = max(1.2maximum(TSH), 5.5)
        Blim = 1.2maximum(Bcell)
        Plim = 1.2maximum(PlasmaCell)
        Tlim = 1.2maximum(CD4Cell)
        Clim = 1.2maximum(Cytokines)
        Flim = 1.2maximum(FTS)
        Alim = 1.2maximum(TPOAb)
    end

    p1 = plot(time / 24.0, T4, ylim=(0, t4lim), xlim=xlim, linewidth=2, label="T4 Initial", ylabel="T4 (mcg/L)")
    p1 = scatter!(tvals, (777.0/3.2).*data.T4, ylabel="T4 (mcg/L)", markersize=1, label="")    
    p1 = plot!(time_fit / 24.0, T4_fit, ylim=(0, t4lim), xlim=xlim, linewidth=2, label="T4 Fit", ylabel="T4 (mcg/L)", color=:darkorchid)
    p2 = plot(time / 24.0, T3, ylim=(0, t3lim), xlim=xlim, linewidth=2, label="T3 Initial", ylabel="T3 (mcg/L)", title=title)
    p2 = scatter!(tvals, (651.0/3.2).*data.T3, ylabel="T3 (mcg/L)", markersize=1, label="")    
    p2 = plot!(time_fit / 24.0, T3_fit, ylim=(0, t3lim), xlim=xlim, linewidth=2, label="T3 Fit", ylabel="T3 (mcg/L)", title=title, color=:darkorchid)
    p3 = plot(time / 24.0, TSH, ylim=(0, tshlim), xlim=xlim, linewidth=2, label="TSH Initial", ylabel="TSH (mU/L)")
    p3 = scatter!(tvals, (5.6/5.2).*data.TSH, ylabel="TSH (mU/L)", markersize=1, label="")
    p3 = plot!(time_fit / 24.0, TSH_fit, ylim=(0, tshlim), xlim=xlim, linewidth=2, label="TSH Fit", ylabel="TSH (mU/L)", color=:darkorchid)
    p4 = plot(time / 24.0, Bcell, label="Bcell Initial", ylabel="Bcell (cell/mL)", ylim=(0,Blim), xlim=xlim, linewidth=2)
    p4 = scatter!(tvals, data.Bcells, ylabel="Bcell (cell/mL)", markersize=1, label="")
    p4 = plot!(time_fit / 24.0, Bcell_fit, label="Bcell Fit", ylabel="Bcell (cell/mL)", ylim=(0,Blim), xlim=xlim, linewidth=2, color=:darkorchid)
    p5 = plot(time / 24.0, PlasmaCell, label="PlasmaCell Initial", ylabel="Plasma Cell (cell/mL)", ylim=(0,Plim), xlim=xlim, linewidth=2)
    p5 = scatter!(tvals, data.Plasma, ylabel="Plasma Cell (cell/mL)", markersize=1, label="")
    p5 = plot!(time_fit / 24.0, PlasmaCell_fit, label="PlasmaCell Fit", ylabel="Plasma Cell (cell/mL)", ylim=(0,Plim), xlim=xlim, linewidth=2, color=:darkorchid)
    p6 = plot(time / 24.0, CD4Cell, label="CD4 Initial", ylabel="CD4+ Cell (cell/mL)", ylim=(0,Tlim), xlim=xlim, linewidth=2)
    p6 = scatter!(tvals, data.CD4, ylabel="CD4+ Cell (cell/mL)", markersize=1, label="")
    p6 = plot!(time_fit / 24.0, CD4Cell_fit, label="CD4 Fit", ylabel="CD4+ Cell (cell/mL)", ylim=(0,Tlim), xlim=xlim, linewidth=2, color=:darkorchid)
    p7 = plot(time / 24.0, Cytokines, label="Cytokine Initial", ylabel="Cytokines (nM)", xlabel="time [days]", ylim=(0,Clim), xlim=xlim, linewidth=2)
    p7 = scatter!(tvals, data.Cytokine./6.022e11, ylabel="Cytokines (nM)", markersize=1, label="")
    p7 = plot!(time_fit / 24.0, Cytokines_fit, label="Cytokine Fit", ylabel="Cytokines (nM)", xlabel="time [days]", ylim=(0,Clim), xlim=xlim, linewidth=2, color=:darkorchid)
    p8 = plot(time / 24.0, FTS, label="FTS Initial", ylabel="Functional Thyroid Size (mL)", xlabel="time [days]", ylim=(0,Flim), xlim=xlim, linewidth=2)
    p8 = scatter!(tvals, data.FTS, ylabel="Functional Thyroid Size (mL)", markersize=1, label="")
    p8 = plot!(time_fit / 24.0, FTS_fit, label="FTS Fit", ylabel="Functional Thyroid Size (mL)", xlabel="time [days]", ylim=(0,Flim), xlim=xlim, linewidth=2, color=:darkorchid)
    p9 = plot(time / 24.0, TPOAb, label="TPOAb Initial", ylabel="TPOAb (pM)", xlabel="time [days]", ylim=(0,Alim), xlim=xlim, linewidth=2)
    p9 = scatter!(tvals, data.Ab./6.022e8, ylabel="TPOAb (pM)", markersize=1, label="")
    p9 = plot!(time_fit / 24.0, TPOAb_fit, label="TPOAb Fit", ylabel="TPOAb (pM)", xlabel="time [days]", ylim=(0,Alim), xlim=xlim, linewidth=2, color=:darkorchid)

    plot!(size=(900,900))
    plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, layout=(3, 3))
end

"""
Run parameter estimation on immune parameters. 
"""
function CV_estim(data::DataFrame, p0::Vector, lb::Vector, ub::Vector)

    t = data.t
    ic = ics()    
    sol_index = [1,4,7,20,21,22,23,24,25]
    cols = names(data)[2:10]

    function error(p,t,data)
        tspan = (t[1],t[end])
        prob = ODEProblem(thyrosimIM_estimate, ic, tspan, p)
        sol = solve(prob, Rosenbrock23())

        loss = 0
        for i in 1:size(cols)[1]
            dp = 1
            for time in t
                loss += ((sol(time)[sol_index[i]]-data[!,cols[i]][dp])/mean(data[!,cols[i]]))^2
                dp += 1
            end
        end
        return loss
    end

    function set_lower_bounds(p::Vector, lb::Vector)
        for i in eachindex(p)
            if p[i] < lb[i]
                p[i] = lb[i]
            end
        end
        return p
    end

    function set_upper_bounds(p::Vector, ub::Vector)
        for i in eachindex(p)
            if p[i] >= ub[i]
                p[i] = ub[i]
            end
        end
        return p
    end

    function objective(p_free, lb, ub)
        p_free = set_lower_bounds(p_free, lb)
        p_free = set_upper_bounds(p_free, ub)
        return error(p_free, t, data)
    end

    # Adjust so fitting indicies and parameters can be specified
    result = optimize(p -> objective(p, lb, ub), p0, Newton(),
     Optim.Options(time_limit = 100.0, iterations = 10, g_tol=1e-5, store_trace=true,
     show_trace = true, allow_f_increases=true, extended_trace=true))
end