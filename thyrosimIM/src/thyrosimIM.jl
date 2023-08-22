module thyrosimIM

    export thyrosimIM!, fixed_parameters, initialize_free, initialize_all, ics # Aidan's model
    export meha_model!, meha_free, meha_ics # Meha's model
    export plot_patient, optim_plotter, output_plotIM # Plotting Functions

    using DifferentialEquations
    using Optim
    using Plots
    using DataFrames
    using LinearAlgebra
    using CSV
    using Statistics
    using MuladdMacro

    # Add includes here...
    include("initialize.jl")
    include("ODEs.jl")
    include("plotters.jl")

end # module thyrosimIM
