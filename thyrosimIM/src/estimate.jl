using DataFrames
import NaNMath # handles log of negative number that can occur due to adaptive ODE solvers

function demeaned_neg_logl(sol, time, data, variances)
    loss = zeros(size(variances, 1))
    data_columns = names(data)[2:10] # this selects T4,T3, TSH, and 6 immune SV... will need to change 
    sol_indicies = [1,4,7,20,21,22,23,24,25] # indicies for T4, T3, TSH, and 6 immune SV
    for (i, column) in enumerate(data_columns)
        column_mean = mean(skipmissing(data[!, column]))
        for (j, t) in enumerate(time)
            datapoint = data[!, column][j]
            if ismissing(datapoint)
                loss = loss
            else
                loss[i] += ((sol(t)[sol_indicies[i]] - datapoint)/column_mean)^2
            end
        end
        loss[i] /= 2variances[i]^2
        loss[i] += length(time) * log(2pi) / 2 + length(time) * NaNMath.log(variances[i])
    end
    return sum(loss)
end

function lsq_loss(sol, time, data, _)
    loss = 0
    data_columns = names(data)[2:10] # this selects T4,T3, TSH, and 6 immune SV... will need to change 
    sol_indicies = [1,4,7,20,21,22,23,24,25] # indicies for T4, T3, TSH, and 6 immune SV
    for (i, column) in enumerate(data_columns)
        column_mean = sum(collect(skipmissing(data[!, column])))/length(collect(skipmissing(data[!, column])))
        for (j, t) in enumerate(time)
            datapoint = data[!, column][j]
            if ismissing(datapoint)
                loss = loss
            else
                loss += ((sol(t)[sol_indicies[i]] - datapoint)/column_mean)^2
            end
        end
    end
    return loss
end

function objective(
    free_parameters::Vector,
    fixed_parameters::Vector,
    lb::Vector,
    ub::Vector,
    data::DataFrame;
    free_indicies::Vector=[],
    loss_function::Function=neg_logl,
    variances::Vector=free_parameters[end-8:end],
    solver=Rosenbrock23())

    # assign t from data, initialize thyrosim problem and compute solution
    t = data.t; tspan = (t[1], t[end])
    problem = ODEProblem(thyrosimIM!, ics(),
        tspan, combine_params(free_parameters[1:end-9], fixed_parameters, free_indicies=free_indicies))
    sol = solve(problem, solver)

    # For now, we are fitting 9 state variables (1,4,7, 20:25). Thus, last 9 parameters must be
    # corresponding variance of the state variables for MLE, can make variable SV number in future
    loss = loss_function(sol, t, data, variances)

    return loss
end

function fit_all(
    free_parameters_guess::Vector,
    fixed_parameters::Vector,
    lb::Vector,
    ub::Vector,
    data::DataFrame;
    free_indicies::Vector=[],
    loss_function::Function=demeaned_neg_logl)

    if isempty(free_indicies)
        fixed_indicies = collect(eachindex(fixed_parameters))
        free_indicies = setdiff(1:(size(fixed_indicies,1)+size(free_parameters_guess,1)), fixed_indicies)
    end

    return optimize(free_params -> objective(free_params, fixed_parameters, lb, ub, data,
        free_indicies=free_indicies, loss_function=loss_function),
        free_parameters_guess, NelderMead(), 
        Optim.Options(time_limit = 500.0, iterations = 100000, g_tol=1e-5,
        show_trace = false, allow_f_increases=true))
end

function execute_fitting(data, loss_function)
    free_params = initialize_free()
    fixed_params = fixed_parameters()
    lb = zeros(size(free_params)[1])
    ub = free_params*10

    solution = fit_all(free_params, fixed_params, lb, ub, data, loss_function=loss_function)
    return solution
end

function CV_estimation(
    free_parameters_guess::Vector,
    fixed_parameters::Vector,
    lb::Vector,
    ub::Vector,
    data::DataFrame;
    free_indicies::Vector=[],
    loss_function::Function=demeaned_neg_logl,
    variances::Vector=[])

    if isempty(free_indicies)
        fixed_indicies = collect(eachindex(fixed_parameters))
        free_indicies = setdiff(1:(size(fixed_indicies,1)+size(free_parameters_guess,1)), fixed_indicies)
    end

    optimal_sol = optimize(free_params -> objective(free_params, fixed_parameters,
        lb, ub, data, free_indicies=free_indicies, loss_function=lsq_loss),
        free_parameters_guess, Newton(), 
        Optim.Options(time_limit = 100.0, iterations = 0, g_tol=1e-5,
        show_trace = true, allow_f_increases=true, extended_trace=true))

    return optimal_sol
    #---------------------------------------------------#
    hessian = optimal_sol.trace[end].metadata["h(x)"]

    inv_hessian = try
        inv(hessian)
        return inv(hessian)
    catch problem
        if typeof(problem) == SingularException
            println("Inversion failed. Hessian is singular")
            return hessian
        else
            println(problem)
        end
    end
end