# PKGs
include("GrowthModels.jl")
using ModelingToolkit, DifferentialEquations, Main.GrowthModels, Plots, DataFrames, CSV, Combinatorics, OrderedCollections
plotly()

##################
# generate differential equations for base model
base_eqs_dict, base_lam, base_params, base_ss = Main.GrowthModels.base_model_eqs();

# compose the base model from proteome fractions
base_model = Main.GrowthModels.compose_model(base_eqs_dict)

# solve the base model
prob = ODEProblem(base_model, Main.GrowthModels._format(base_ss), (0, 1e5), Main.GrowthModels._format(base_params); jac=true);
sol  = solve(prob, Rodas4());

# update with steady state values
Main.GrowthModels.update_multiple!(Main.GrowthModels.map_end_vals(base_model, sol), base_ss)

##################
# simulate host + heterologous
# generate differential equations for base model extended with a single heterologous protein
ha_eqs, ha_lam, ha_params, ha_ss = Main.GrowthModels.het_model_eqs(input_eqs_dict   = base_eqs_dict, 
                                                                   input_lam        = base_lam, 
                                                                   input_param_vals = base_params,
                                                                   input_ss_vals    = base_ss);
ha_het_model = Main.GrowthModels.compose_model(ha_eqs)

prob_h = ODEProblem(ha_het_model, Main.GrowthModels._format(ha_ss), (0, 1e5), Main.GrowthModels._format(ha_params); jac=true);
sol_h  = solve(prob_h, Rodas4());
# plot only the solution for the heterologous expression -- a bit boring 
plot(sol_h)

# update heterologous construct parameters from heter_model
# re-execute lines 28-29 to solve a new model
Main.GrowthModels._update!(ha_params, "w_max", 100.0)


#-----------#

# needs to be tested
#res_dict, perturb_vals = perturb_one_param!(perturb_info = perturb_info, model_info = model_info, solver_info = solver_info);

# Define a repressilator genetic circuit
# NOTE: there is definitely an issue with formalism when it comes to stacking
# I probably need to pass in a vector of variable names
# let uss first do it the stupid way

repr_eqs, repr_lam, repr_params, repr_ss = Main.GrowthModels.repr_model_eqs(input_eqs_dict   = base_eqs_dict,
                                                                            input_lam        = base_lam,
                                                                            input_param_vals = base_params,
                                                                            input_ss_vals    = base_ss);

repr_model = Main.GrowthModels.compose_model(repr_eqs)

# solve the model once
prob_repr = ODEProblem(repr_model, Main.GrowthModels._format(repr_ss), (0, 1e5), Main.GrowthModels._format(repr_params); jac=true);
sol_repr  = solve(prob_repr, Rodas4(autodiff=false));
# plot all species
plot(sol_repr)

# figure out why I did this
# line updates all steady state values
Main.GrowthModels.update_multiple!(Main.GrowthModels.map_end_vals(repr_model, sol_repr), repr_ss)

# exports full time trajectories for each species
test_map = Main.GrowthModels.map_trajectories(repr_model, sol_repr)

# test multi parameter perturbations
# WORKS for updates of same values for different parameters, e.g. max induction strength

include("perturbation_values.jl")

function perturb_multi_single_params!(; perturb_info::Dict, model_info::Dict, solver_info::Dict)
    # initialize results vector and iterator
    rng_iterator = exp10.(range(perturb_info["range_to_perturb"][1],
                                perturb_info["range_to_perturb"][2], 
                                length = perturb_info["range_size"]))

    # initialize the format dictionary to return                                        
    format_dict = Dict(string(state)[1:end-3] => [] for state in states(model_info["model_name"]))

        # Main loop
        for val in rng_iterator 

            if length(perturb_info["param_to_perturb"]) > 1
                for ii in range(1, length(perturb_info["param_to_perturb"]))
                    Main.GrowthModels._update!(model_info["model_parameters"],
                                               perturb_info["param_to_perturb"][ii],
                                               val
                                               )
                end
            else
                # update parameter value
                Main.GrowthModels._update!(model_info["model_parameters"], 
                                        perturb_info["param_to_perturb"], 
                                        val)
            end
            # call to ODEProblem
            prob = ODEProblem(model_info["model_name"], 
                              Main.GrowthModels._format(model_info["steady_state_values"]), 
                              model_info["integr_time"], 
                              Main.GrowthModels._format(model_info["model_parameters"]); 
                              jac=solver_info["jacobian_opt"]);
            # call to solver
            sol = solve(prob, solver_info["solver_option"]);

            # map species to solutions
            species_to_solutions = (Main.GrowthModels.map_trajectories(model_info["model_name"], sol))

            # concatenate results
            for key in keys(format_dict)
                format_dict[key] = push!(format_dict[key], species_to_solutions[key])
            end

        end

    return format_dict, rng_iterator
end

res_dict, perturb_vals = perturb_multi_single_params!(perturb_info = perturb_info, model_info = model_info, solver_info = solver_info);

function perturb_multi_params!(; perturb_info::Dict, model_info::Dict, solver_info::Dict)
    # initialize results vector and iterator
    range_iterator_dict = OrderedDict(key => exp10.(range(perturb_info["range_to_perturb"][key][1],
                                                        perturb_info["range_to_perturb"][key][2],
                                                        length = perturb_info["range_size"])) for key in keys(perturb_info["param_to_perturb"]))

    combinations_list = allcombinations_(values(range_iterator_dict)...)

    # initialize the format dictionary to return                                        
    format_dict = Dict(string(state)[1:end-3] => [] for state in states(model_info["model_name"]))

    for combi in combinations_list
        for (idx, key) in enumerate(keys(perturb_info["param_to_perturb"]))
            for idx_param in range(1, length(perturb_info["param_to_perturb"][key]))
                Main.GrowthModels._update!(model_info["model_parameters"],
                                        perturb_info["param_to_perturb"][key][idx_param],
                                        combi[idx])
            end
        end

        # call to ODEProblem
        prob = ODEProblem(model_info["model_name"], 
                        Main.GrowthModels._format(model_info["steady_state_values"]), 
                        model_info["integr_time"], 
                        Main.GrowthModels._format(model_info["model_parameters"]), 
                        jac=solver_info["jacobian_opt"]);
        # call to solver
        sol = solve(prob, solver_info["solver_option"]);

        # map species to solutions
        species_to_solutions = (Main.GrowthModels.map_trajectories(model_info["model_name"], sol))

        # concatenate results
        for key in keys(format_dict)
            format_dict[key] = push!(format_dict[key], species_to_solutions[key])
        end

    end

    return format_dict, combinations_list, range_iterator_dict
end

res_dict, combi_list, range_itr_dict = Main.GrowthModels.perturb_multi_params!(perturb_info = perturb_info, model_info = model_info, solver_info = solver_info);

# turn results dictionary into a dataframe
df = DataFrame(res_dict)

# save the dataframe into a CSV 
save_to_file = false
file_name = "results_perturbation_induction_RBS.csv"
if save_to_file
    CSV.write(file_name, df)
end

# plot
plot(df.p_h_1, df.p_h_2, df.p_h_3)