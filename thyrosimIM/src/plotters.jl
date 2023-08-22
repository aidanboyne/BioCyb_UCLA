"""
Emulate data_plotter.py plots used for data selection to double check import success.
"""
function plot_patient(patient; title::String="Patient Labs")
    labs = patient[1]
    meds = patient[2]

    T4T3Labels = ["FT4 (ng/dL)", "FT3 (pg/dL)"]
    T4T3Ranges = [[0.9, 2.3], [260, 480]]
    xlims = (min(minimum(labs.t)/24, minimum(meds.t)/24), max(maximum(labs.t)/24, maximum(meds.t))/24)

    t = labs.t/24
    
    p1 = plot(t, labs.T4, label="T4", xlim=xlims, ylabel=T4T3Labels[1])
    any(meds.Levothyroxine .!= 0) && (p1 = plot!(twinx(), meds.Levothyroxine, xlim=xlims, ylabel="Levothyroxine (mcg/day)", label="Dose", color = "red", xticks=:none))
    p1 = hline!(T4T3Ranges[1], label= "")
    p1 = scatter!(t, labs.T4, alpha = 0.3, label="")
    
    p2 = plot(t, labs.T3, label="", xlim=xlims, ylabel=T4T3Labels[2])
    any(meds.Lyothyronine .!= 0) && (p2 = plot!(twinx(), meds.Lyothyronine, xlim=xlims, ylabel="Liothyronine (mcg/day)", label = "Dose", color = "red", xticks=:none))
    p2 = hline!(T4T3Ranges[2], label= "")
    p2 = scatter!(t, labs.T3, alpha = 0.3, label="")
    
    p3 = plot(t, labs.TSH, label="", xlim=xlims, ylabel="TSH (mU/L)")
    p3 = hline!([0.45, 4.5], label= "")
    p3 = scatter!(t, labs.TSH, alpha = 0.3)

    p4 = plot(t, labs.Lymphocytes, xlim=xlims, ylabel="Lymphocytes (cell/mL)", label="")
    p4 = scatter!(t, labs.Lymphocytes, alpha = 1, label="")

    p5 = plot(t, labs.Ab, ylabel="TPOAb (IU/mL)", xlim=xlims, xlabel="time [days]", label="")
    p5 = scatter!(t, labs.Ab, alpha = 1, label="")

    layout = @layout [a b ; c ; d e]
    plot!(size=(900,900))
    p_main = plot(p1, p2, p3, p4, p5, layout=layout)
    p_main[:plot_title] = title
    plot(p_main)

end

"""
Default plotter for ThyrosimIM optimization results. Displays FT4, FT3, TSH, Lymphocytes (with breakdown), and TPOAb levels.
"""
function output_plotIM(sol, data; title::AbstractString = "ThyrosimIM simulation", automargins::Bool=true, save_to_file = false, plot_no::Int = 1)

    # can definitley clean up this code quite a bit later (looped plotting, get rid of unnecessary plots...)
    # parameters to adjust figure limits
    p = sol.prob.p 
    t4lim, t3lim, tshlim = 140, 4, 10
    T4 = FT4.(sol[1, :])
    T3 = FT3.(sol[4, :])
    T4T3Labels = ["FT4 (ng/dL)", "FT3 (pg/dL)"]
    T4T3Ranges = [[0.9, 2.3], [260, 480]]
    TSH = TSHfit.(sol[7, :])
    Bcell = sol[20, :]
    Pcell = sol[21, :]
    Tcell = sol[22, :]
    Lymphocytes = sol[20, :] + sol[21, :] + sol[22, :]
    TPOAb = TPOConvert.(sol[25, :]) # convert to pM from molecules/mL

    xlim=(0,sol.t[end]/24) 
    if automargins
        t4lim = max(1.2*maximum(T4), T4T3Ranges[1][2]*1.1)
        t3lim = max(1.2*maximum(T3), T4T3Ranges[2][2]*1.1)
        tshlim = max(1.2maximum(TSH), 5.5)
        Llim = 1.2maximum(Lymphocytes)
        Alim = 1.2maximum(TPOAb)
    end

    p1 = plot(sol.t / 24.0, T4, ylim=(0, t4lim), xlim=xlim, label="", ylabel=T4T3Labels[1])
    p1 = hline!(T4T3Ranges[1], label= "")
    p1 = scatter!(data.t./24, data.T4, alpha = 0.3, label= "")
    
    p2 = plot(sol.t / 24.0, T3, ylim=(0, t3lim), xlim=xlim, label="", ylabel=T4T3Labels[2])
    p2 = hline!(T4T3Ranges[2], label= "")
    p2 = scatter!(data.t./24, data.T3, alpha = 0.3, label= "")
    
    p3 = plot(sol.t / 24.0, TSH, ylim=(0, tshlim), xlim=xlim, label="", ylabel="TSH (mU/L)")
    p3 = hline!([0.45, 4.5], label= "")
    p3 = scatter!(data.t./24, data.TSH, alpha = 0.3, label= "")

    p4 = plot(sol.t / 24.0, Lymphocytes, ylabel="Lymphocytes (cell/mL)", ylim=(0,Llim), xlim=xlim, label="")
    p4 = scatter!(data.t./24, data.Lymphocytes, alpha = 0.3, label= "")
    p4 = plot!(sol.t /24.0, Bcell, label = "B-Cells")
    p4 = plot!(sol.t /24.0, Tcell, label = "T-Cells")
    p4 = plot!(sol.t /24.0, Pcell, label = "Plasma Cells")

    p5 = plot(sol.t / 24.0, TPOAb, ylabel="TPOAb (IU/mL)", xlabel="time [days]", ylim=(0,Alim), xlim=xlim, label="")
    p5 = scatter!(data.t./24, data.Ab, alpha = 0.3, label= "")

    #p6 = plot(dose_df.t / 24.0, dose_df.Levothyroxine, ylabel="Oral Levothyroxine (mcg/day)", xlabel = "time [days]", label="")
    
    l = @layout [a b ; c ; d e]
    plot!(size=(900,900))
    p_main = plot!(p1, p2, p3, p4, p5, layout=l)
    p_main[:plot_title] = title
    if save_to_file
        saveplot = plot(p_main)
        savefig(saveplot, "img/optimization_prog_$(plot_no)")
    else
        plot(p_main)
    end
end

"""
Plotter for in-optimization solution tracking.
"""
function optim_plotter(sol, data; title::AbstractString = "ThyrosimIM simulation", automargins::Bool=true, save_to_file = false, plot_no::Int = 1)

    # can definitley clean up this code quite a bit later (looped plotting, get rid of unnecessary plots...)
    # parameters to adjust figure limits
    T4 = FT4.(sol[1, :])
    T3 = FT3.(sol[4, :])
    TSH = TSHfit.(sol[7, :])
    Bcell = sol[20, :]
    Pcell = sol[21, :]
    Tcell = sol[22, :]
    Lymphocytes = sol[20, :] + sol[21, :] + sol[22, :]
    TPOAb = TPOConvert.(sol[25, :]) # convert to pM from molecules/mL
    T4T3Ranges = [[0.9, 2.3], [260, 480]]

    xlim=(0,sol.t[end]/24) 
    if automargins
        t4lim = max(1.2*maximum(T4), T4T3Ranges[1][2]*1.1)
        t3lim = max(1.2*maximum(T3), T4T3Ranges[2][2]*1.1)
        tshlim = max(1.2maximum(TSH), 5.5)
        Llim = 1.2maximum(Lymphocytes)
        Alim = 1.2maximum(TPOAb)
    end

    p1 = plot(sol.t / 24.0, T4, ylim=(0, t4lim), xlim=xlim, label="", ylabel="FT4 (ng/dL)")
    p1 = hline!(T4T3Ranges[1], label= "")
    p1 = scatter!(data.t./24, data.T4, alpha = 0.3, label= "")
    
    p2 = plot(sol.t / 24.0, T3, ylim=(0, t3lim), xlim=xlim, label="", ylabel="FT3 (pg/dL)")
    p2 = hline!(T4T3Ranges[2], label= "")
    p2 = scatter!(data.t./24, data.T3, alpha = 0.3, label= "")
    
    p3 = plot(sol.t / 24.0, TSH, ylim=(0, tshlim), xlim=xlim, label="", ylabel="TSH (mU/L)")
    p3 = hline!([0.45, 4.5], label= "")
    p3 = scatter!(data.t./24, data.TSH, alpha = 0.3, label= "")

    p4 = plot(sol.t / 24.0, Lymphocytes, ylabel="Lymphocytes (cell/mL)", ylim=(0,Llim), xlim=xlim, label="")
    p4 = scatter!(data.t./24, data.Lymphocytes, alpha = 0.3, label= "")
    p4 = plot!(sol.t /24.0, Bcell, label = "B-Cells")
    p4 = plot!(sol.t /24.0, Tcell, label = "T-Cells")
    p4 = plot!(sol.t /24.0, Pcell, label = "Plasma Cells")

    p5 = plot(sol.t / 24.0, TPOAb, ylabel="TPOAb (IU/mL)", xlabel="time [days]", ylim=(0,Alim), xlim=xlim, label="")
    p5 = scatter!(data.t./24, data.Ab, alpha = 0.3, label= "")

    #p6 = plot(dose_df.t / 24.0, dose_df.Levothyroxine, ylabel="Oral Levothyroxine (mcg/day)", xlabel = "time [days]", label="")
    
    l = @layout [a b ; c ; d e]
    plot!(size=(900,900))
    p_main = plot!(p1, p2, p3, p4, p5, layout=l)
    p_main[:plot_title] = title
    if save_to_file
        saveplot = plot(p_main)
        savefig(saveplot, "img/optimization_prog_$(plot_no)")
    else
        plot(p_main)
    end
end

