# growth_models_new_implementation

driver.jl is the main workflow for simulating a host-aware system (more below)

GrowthModel.jl contains model definitions along with utility functions for updating and formatting simulation values

perturbation_values.jl contains dictionaries for defining parameter ranges, model and solver specifications, for running multiple simulations of a host-aware repressilator

## driver.jl details

run lines 1-18 to simulate a single host model (without any heterologous circuit) -- important: line 18 updates all species with the steady-state values resulting from solving the base model; those values are used as starting conditions for any host-aware system
ignore lines 19-45 for now

run lines 46-57 to siimulate a single host-aware repressilator system -- line 57 should generate an interactive plotly figure

ignore lines 58-68 for now

run lines 69-76 to simulate the joint effects of changes in maximal induction strength, binding and unbinding ribosomal rate for the host-aware repressilator -- resulting time trajectories for all species should be stored in a dataframe

set line 79 to true to save the dataframe as a CSV

run line 86 to get an iteractive plot of the resulting limit cycles for all combinations of perturbations


## TO-DO
- list of packages
- better way to define and pass dirs with perturbation specifications
- update model definition to include the translation initiation efficiency parameter
- provide more options for the ODE solver