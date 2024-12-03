function single_gene_perturbations(; param_to_perturb::String, range_to_perturb::Tuple, range_size::Int, time_range::Tuple = (0, 1e5))

    # organize perturbation related info in a dictionary
    perturb_info = Dict(
        "param_to_perturb" => param_to_perturb,
        "range_to_perturb" => range_to_perturb,
        "range_size"       => range_size
    )

    # organize ODESystem related info in a dictionary
    model_info = Dict(
        "model_name"          => ha_het_model,
        "steady_state_values" => ha_ss,
        "model_parameters"    => ha_params,
        "integr_time"         => time_range
    )
    # organize solver related info in a dictionary
    solver_info = Dict(
        "solver_option" => Rodas4(),
        "jacobian_opt"  => true,
        "max_iters"      => 1e5,
        "progress"      => false,
        "abstol"        => 1e-8,
        "reltol"        => 1e-8,
        "isoutofdomain" => (m,p,t) -> any(x->x<0, m)
    )

    return perturb_info, model_info, solver_info
end

function repressilator_perturbations()
    # time to integrate the system
    time_integr = (0, 1e5)

    # use common keys for parameters and respective ranges

    # dictionary for perturbation-related info
    perturb_info =  Dict(
        "param_to_perturb" => Dict("gene_induction" => ["w_max_1", "w_max_2", "w_max_3"],
                                "binding_rate"      => ["kb_h_1", "kb_h_2", "kb_h_3"],
                                #"unbinding_rate"   => ["ku_h_1", "ku_h_2", "ku_h_3"]
                                ),
        "range_to_perturb" => Dict("gene_induction" => (0, 4),
                                "binding_rate"      => (-2, 0),
                                #"unbinding_rate"   => (-4, -2)
                                ),
        "range_size" => 5
    )

    # dictionary for ODESystem-related info
    model_info = Dict(
        "model_name"          => repr_model,
        "steady_state_values" => repr_ss,
        "model_parameters"    => repr_params,
        "integr_time"         => time_integr
    )

    # organize solver related info in a dictionary
    solver_info = Dict(
        "solver_option" => Rodas4(),
        "jacobian_opt"  => true,
        "max_iters"      => 1e7,
        "progress"      => false,
        "abstol"        => 1e-8,
        "reltol"        => 1e-8,
        "isoutofdomain" => (m,p,t) -> any(x->x<0, m)
    )

    return time_integr, perturb_info, model_info, solver_info
end

function perturb_one_param!(; perturb_info::Dict, model_info::Dict, solver_info::Dict)
    # initialize results vector and iterator
    rng_iterator = exp10.(range(perturb_info["range_to_perturb"][1],
                                perturb_info["range_to_perturb"][2], 
                                length = perturb_info["range_size"]))

    # initialize the format dictionary to return                                        
    format_dict = Dict(string(state)[1:end-3] => [] for state in states(model_info["model_name"]))

        # Main loop
        for val in rng_iterator 

            # update parameter value
            Main.GrowthModels._update!(model_info["model_parameters"], 
                                       perturb_info["param_to_perturb"], 
                                       val)

            # call to ODEProblem
            prob = ODEProblem(model_info["model_name"], 
                              Main.GrowthModels._format(model_info["steady_state_values"]), 
                              model_info["integr_time"], 
                              Main.GrowthModels._format(model_info["model_parameters"]); 
                              jac=solver_info["jacobian_opt"]);
            # call to solver
            sol = solve(prob, 
                        solver_info["solver_option"], 
                        abstol=solver_info["abstol"],
                        reltol=solver_info["reltol"], 
                        maxiters=solver_info["max_iters"], 
                        progress=solver_info["progress"], 
                        isoutofdomain = solver_info["isoutofdomain"]);

            # map species to solutions
            species_to_solutions = (Main.GrowthModels.map_trajectories(model_info["model_name"], sol))

            # concatenate results
            for key in keys(format_dict)
                format_dict[key] = push!(format_dict[key], species_to_solutions[key])
            end

        end

    return format_dict, rng_iterator
end