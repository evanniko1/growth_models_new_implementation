# Define dictionarties for system perturbations

# time to simulate the system
time_integr = (0, 1e5)

# use common keys for parameters and respective ranges

# dictionary for perturbation-related info
perturb_info =  Dict(
    "param_to_perturb" => Dict("gene_induction" => ["w_max_1", "w_max_2", "w_max_3"],
                               "binding_rate"   => ["kb_h_1", "kb_h_2", "kb_h_3"],
                               #"unbinding_rate" => ["ku_h_1", "ku_h_2", "ku_h_3"]
                               ),
    "range_to_perturb" => Dict("gene_induction" => (0, 4),
                               "binding_rate"   => (-2, 0),
                               #"unbinding_rate" => (-4, -2)
                               ),
    "range_size" => 10
)

# dictionary for ODESystem-related info
model_info = Dict(
    "model_name"          => repr_model,
    "steady_state_values" => repr_ss,
    "model_parameters"    => repr_params,
    "integr_time"         => time_integr
)
# dictionary for solver-related info
solver_info = Dict(
    "solver_option" => Rodas4(autodiff=false),
    "jacobian_opt"  => true,
    #"solver_conditions" -> nothing
)

######
# code to check

#time_integr = (0, 1e5)

# organize perturbation related info in a dictionary
#perturb_info = Dict(
#    "param_to_perturb" => "w_max",
#    "range_to_perturb" => (0, 3),
#    "range_size" => 25
#)

# organize ODESystem related info in a dictionary
#model_info = Dict(
#    "model_name" => ha_het_model,
#    "steady_state_values" => ha_ss,
#    "model_parameters" => ha_params,
#    "integr_time" => time_integr
#)
# organize solver related info in a dictionary
#solver_info = Dict(
#    "solver_option" => Rodas4(),
#    "jacobian_opt" => true,
    #"solver_conditions" -> nothing
#)

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
            sol = solve(prob, solver_info["solver_option"]);

            # map species to solutions
            species_to_solutions = (map_vals_TEST(model_info["model_name"], sol))

            # concatenate results
            for key in keys(format_dict)
                format_dict[key] = push!(format_dict[key], species_to_solutions[key])
            end

        end

    return format_dict, rng_iterator
end

function get_results(; equations_dict::Dict, solutions_dict::Dict)
    # initialize the fractions dictionary
    fractions_dict = Dict()

    for species_key in keys(equations_dict)
        fractions_dict[species_key] = Dict(string(equations_dict[species_key][idx].lhs)[17:end-4] => [] for idx in range(1, size(equations_dict[species_key])[1]))
    end

    for frac_key in keys(fractions_dict)
        for sol_key in keys(solutions_dict)
            if sol_key in keys(fractions_dict[frac_key])
                fractions_dict[frac_key][sol_key] = solutions_dict[sol_key]
            end
        end
    end

    return fractions_dict
end

# GENERATE AN EXTRA DICT THAT INCLUDES (STEADY STATE (END) VALUES, VALUE USED) (??)
# initialize fractions dict
fractions_dict = Dict()
# use the left hand side of the ODE system definition to assign each species to the correct fraction
for fraction_key in keys(ha_eqs)
    fractions_dict[fraction_key] = Dict(string(ha_eqs[fraction_key][idx].lhs)[17:end-4] => [] for idx in range(1, size(ha_eqs[fraction_key])[1]))
end

# go over the dictionary containing solutions of the ODE system and organize accordingly
for key_frac in keys(fractions_dict)
    for key_res in keys(res_dict)
        if key_res in keys(fractions_dict[key_frac])
            fractions_dict[key_frac][key_res] = res_dict[key_res]
        end
    end
end