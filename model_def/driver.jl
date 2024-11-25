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

#res_dict, perturb_vals = perturb_multi_single_params!(perturb_info = perturb_info, model_info = model_info, solver_info = solver_info);

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