# growth_models_new_implementation

## List of packages

Note: the current setup is for Julia v1.8

Agents v5.6.2

CSV v0.10.9

CairoMakie v0.8.13

Catalyst v12.3.1

ColorSchemes v3.26.0

Combinatorics v1.0.2

DataFrames v1.4.4

DifferentialEquations v7.6.0

Distributions v0.25.79

GraphViz v0.2.0

IJulia v1.23.3

InteractiveDynamics v0.21.11

LSODA v0.7.3

LaTeXStrings v1.3.1

Latexify v0.15.17

MAT v0.10.4

ModelingToolkit v8.34.0

OrderedCollections v1.6.0

OrdinaryDiffEq v6.31.2

ParameterizedFunctions v5.15.0

PlotlyBase v0.8.19

Plots v1.36.3

PyPlot v2.11.0

StatsBase v0.33.21

Sundials v4.15.1

Symbolics v4.14.0

## project breakdown

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