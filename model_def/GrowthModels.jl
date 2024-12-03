module GrowthModels
    
using ModelingToolkit, DifferentialEquations, OrderedCollections, Combinatorics

# global base model parameters and variables
prm = @parameters thetar k_cm s0 gmax thetax Kt M we Km vm nx Kq Kp vt wr wq wp nq nr dm kb ku ns  
var= @variables t mq(t) rmq(t) q(t) mm(t) rmm(t) em(t) mt(t) rmt(t) et(t) mr(t) rmr(t) r(t) si(t) a(t)
#some global calculations
Kgamma = gmax/Kp
gamma  = gmax*a/(Kgamma + a)
nucat  = em*vm*si/(Km + si)

function get_prm_var()
    prm_dict = Dict([
        "thetar" => thetar 
        "k_cm" => k_cm 
        "s0" => s0 
        "gmax" => gmax 
        "thetax" => thetax 
        "Kt" => Kt 
        "M" => M 
        "we" => we 
        "Km" => Km 
        "vm" => vm 
        "nx" => nx
        "Kq" => Kq 
        "Kp" => Kp 
        "vt" => vt 
        "wr" => wr 
        "wq" => wq 
        "wp" => wp 
        "nq" => nq 
        "nr" => nr 
        "dm" => dm 
        "kb" => kb 
        "ku" => ku 
        "ns" => ns 
    ]);
    var_dict  = Dict([
        "rmr" => rmr
        "em"  => em
        "rmq" => rmq
        "rmt" => rmt
        "et"  => et
        "rmm" => rmm
        "mt"  => mt
        "mm"  => mm
        "q"   => q
        "si"  => si
        "mq"  => mq
        "mr"  => mr
        "r"   => r
        "a"   => a
    ]);
    return prm_dict, var_dict
end

function base_model_eqs(; R::Dict = Dict([
    "house_keeping" => 1,
    "metabolic"     => 1,
    "transporter"   => 1,
    "ribosome"      => 1
     ]))

    # unpack regulatory functions
    lam = ((rmq + rmr + rmt + rmm)*gamma)/M
    @variables t
    D = Differential(t)

    eqs_hk_pr = [
        D(mq)  ~ (wq*a/(thetax + a))*R["house_keeping"]+ku*rmq+gamma/nx*rmq-kb*r*mq-dm*mq-lam*mq,
        D(rmq) ~ kb*r*mq-ku*rmq-gamma/nx*rmq-lam*rmq,
        D(q)   ~ gamma/nx*rmq-lam*q
    ];
    eqs_m_pr = [
        D(mm)  ~ (we*a/(thetax + a))*R["metabolic"][1]+ku*rmm+gamma/nx*rmm-kb*r*mm-dm*mm-lam*mm,
        D(rmm) ~ kb*r*mm-ku*rmm-gamma/nx*rmm-lam*rmm,
        D(em)  ~ gamma/nx*rmm-lam*em
    ];
    eqs_t_pr = [
        D(mt)  ~ (we*a/(thetax + a))*R["transporter"][1]+ku*rmt+gamma/nx*rmt-kb*r*mt-dm*mt-lam*mt,
        D(rmt) ~ kb*r*mt-ku*rmt-gamma/nx*rmt-lam*rmt,
        D(et)  ~ gamma/nx*rmt-lam*et
    ];
    eqs_r_mc = [
        D(mr)  ~ (wr*a/(thetar + a))*R["ribosome"][1]+ku*rmr+gamma/nr*rmr-kb*r*mr-dm*mr-lam*mr,
        D(rmr) ~ kb*r*mr-ku*rmr-gamma/nr*rmr-lam*rmr,
    ];
    eqs_r_pr = [
        D(r)   ~ ku*rmr+ku*rmt+ku*rmm+ku*rmq+gamma/nr*rmr+gamma/nr*rmr+gamma/nx*rmt+gamma/nx*rmm+gamma/nx*rmq-kb*r*mr-kb*r*mt-kb*r*mm-kb*r*mq-lam*r
    ];          
    eqs_s = [
        D(si)  ~ (et*vt*s0/(Kt + s0))-nucat-lam*si,
    ];
    eqs_a = [
        D(a)   ~ ns*nucat-lam*M-lam*a
    ];

    # organize all protein fractions into a nice dictionary
    eqs_dict = Dict(
        "house_keeping" => eqs_hk_pr,
        "metabolic"     => eqs_m_pr,
        "transporter"   => eqs_t_pr,
        "ribo_mc"       => eqs_r_mc,
        "ribosomes"     => eqs_r_pr,
        "nutrient"      => eqs_s,
        "energy"        => eqs_a
    );

    ss_values = Dict([
        "rmr" => rmr => 0
        "em"  => em  => 0
        "rmq" => rmq => 0
        "rmt" => rmt => 0
        "et"  => et  => 0
        "rmm" => rmm => 0
        "mt"  => mt  => 0
        "mm"  => mm  => 0
        "q"   => q   => 0
        "si"  => si  => 0
        "mq"  => mq  => 0
        "mr"  => mr  => 0
        "r"   => r   => 10
        "a"   => a   => 1e3
    ]);

    param_values = Dict([
        "thetar" => thetar => 426.8693338968694
        "k_cm" => k_cm => 0.005990373118888
        "s0" => s0 => 10000
        "gmax" => gmax => 1260.0
        "thetax" => thetax => 4.379733394834643
        "Kt" => Kt => 1.0e3
        "M" => M => 1.0e8
        "we" => we => 4.139172187824451
        "Km" => Km => 1000
        "vm" => vm => 5800.0
        "nx" => nx => 300.0
        "Kq" => Kq => 1.522190403737490e+05
        "Kp" => Kp => 180.1378030928276
        "vt" => vt => 726.0
        "wr" => wr => 929.9678874564831
        "wq" => wq => 948.9349882947897
        "wp" => wp => 0.0
        "nq" => nq => 4
        "nr" => nr => 7549.0
        "dm" => dm => 0.1
        "kb" => kb => 0.0095
        "ku" => ku => 1
        "ns" => ns => 0.5
    ]);
    return eqs_dict, lam, param_values, ss_values
end;

function het_model_eqs(; input_eqs_dict::Dict, input_lam::Num, input_param_vals::Dict, input_ss_vals::Dict, R = 1)
    @parameters w_max dm_h dp_h kb_h ku_h kappa_ini
    @variables t m_h(t) c_h(t) p_h(t) r(t) a(t)
    D = Differential(t)

    het_param_vals = Dict([
        "w_max" => w_max => 150
        "dm_h"  => dm_h  => log(2)/2
        "dp_h"  => dp_h  => log(2)/4
        "kb_h"  => kb_h  => 0.0095
        "ku_h"  => ku_h  => 1
        "kappa_ini" => kappa_ini => 0.8
    ])
    het_ss_vals = Dict([
        "m_h" => m_h => 0.0
        "c_h" => c_h => 0.0
        "p_h" => p_h => 0.0
    ])

    # add contribution to growth
    lam = input_lam + (kappa_ini*c_h*gamma)/M

    # update all growth rates for the input system of equations
    upd_eqs_dict = update_eqs!(eqs_dict= input_eqs_dict, term = input_lam - lam);

    # define equations for heterologous species
    eqs_het = [
        D(m_h) ~ R*w_max*a/(thetax+a) - (lam + dm_h)*m_h + (gamma/nx)*kappa_ini*c_h - kb_h*r*m_h + ku_h*(1-kappa_ini)*c_h
        D(c_h) ~ -lam*c_h + kb_h*r*m_h - ku_h*(1-kappa_ini)*c_h - (gamma/nx)*kappa_ini*c_h
        D(p_h) ~ (gamma/nx)*kappa_ini*c_h - (lam + dp_h)*p_h
    ];
    # update energy and ribosomal contributions
    eqs_r_pr_ha = [
        D(r) ~ upd_eqs_dict["ribosomes"][1].rhs - kb_h*m_h*r + ku_h*(1-kappa_ini)*c_h + (gamma/nx)*kappa_ini*c_h
    ];
    eqs_a_ha = [
        D(a) ~ upd_eqs_dict["energy"][1].rhs - gamma*kappa_ini*c_h
    ];

    # organize host_aware model equations in a dictionary
    eqs_dict = Dict(
        "house_keeping" => upd_eqs_dict["house_keeping"],
        "metabolic"     => upd_eqs_dict["metabolic"],
        "transporter"   => upd_eqs_dict["transporter"],
        "ribo_mc"       => upd_eqs_dict["ribo_mc"],
        "ribosomes"     => eqs_r_pr_ha,
        "heterologous"  => eqs_het,
        "nutrient"      => upd_eqs_dict["nutrient"],
        "energy"        => eqs_a_ha
    )

    param_vals = merge(input_param_vals, het_param_vals)
    ss_vals    = merge(input_ss_vals, het_ss_vals)

    return eqs_dict, lam, param_vals, ss_vals
end;

function repr_model_eqs(; input_eqs_dict::Dict, input_lam::Num, input_param_vals::Dict, input_ss_vals::Dict)
    @parameters w_max_1 dm_h_1 dp_h_1 kb_h_1 ku_h_1 Kq_p_1 nq_p_1 w_max_2 dm_h_2 dp_h_2 kb_h_2 ku_h_2 Kq_p_2 nq_p_2 w_max_3 dm_h_3 dp_h_3 kb_h_3 ku_h_3 Kq_p_3 nq_p_3
    @variables t m_h_1(t) c_h_1(t) p_h_1(t) m_h_2(t) c_h_2(t) p_h_2(t) m_h_3(t) c_h_3(t) p_h_3(t) r(t) a(t)
    D = Differential(t)

    het_param_vals =  Dict([
        "w_max_1" => w_max_1 => 150
        "dm_h"  => dm_h_1  => log(2)/2
        "dp_h"  => dp_h_1  => log(2)/4
        "kb_h_1"  => kb_h_1  => exp10(-1.3335)
        "ku_h_1"  => ku_h_1  => exp10(-2.6575)
        "Kq_p_1"=> Kq_p_1  => 100
        "nq_p_1"=> nq_p_1  => 2
        "w_max_2" => w_max_2 => 150
        "dm_h_2"  => dm_h_2  => log(2)/2
        "dp_h_2"  => dp_h_2  => log(2)/4
        "kb_h_2"  => kb_h_2  => exp10(-1.3335)
        "ku_h_2"  => ku_h_2  => exp10(-2.6575)
        "Kq_p_2"  => Kq_p_2  => 100
        "nq_p_2"  => nq_p_2  => 2
        "w_max_3" => w_max_3 => 150
        "dm_h_3"  => dm_h_3  => log(2)/2
        "dp_h_3"  => dp_h_3  => log(2)/4
        "kb_h_3"  => kb_h_3  => exp10(-1.3335)
        "ku_h_3"  => ku_h_3  => exp10(-2.6575)
        "Kq_p_3"  => Kq_p_3  => 100
        "nq_p_3"  => nq_p_3  => 2
    ])

    het_ss_vals = Dict([
        "m_h_1" => m_h_1 => 10.0
        "c_h_1" => c_h_1 => 0.0
        "p_h_1" => p_h_1 => 0.0
        "m_h_2" => m_h_2 => 100.0
        "c_h_2" => c_h_2 => 0.0
        "p_h_2" => p_h_2 => 0.0
        "m_h_3" => m_h_3 => 1000.0
        "c_h_3" => c_h_3 => 0.0
        "p_h_4" => p_h_3 => 0.0
    ])

    # define regulatory functions
    R_1 = (1 / (1 + (p_h_3/Kq_p_3)^nq_p_3))
    R_2 = (1 / (1 + (p_h_1/Kq_p_1)^nq_p_1))
    R_3 = (1 / (1 + (p_h_2/Kq_p_2)^nq_p_2))

    # add contirbution to growth
    lam = input_lam + (c_h_1 + c_h_2 + c_h_3) * gamma/M

    # update all growth rates for the input system of equations
    upd_eqs_dict = update_eqs!(eqs_dict = input_eqs_dict, term = input_lam - lam);

    # define equations for heterologous species
    eqs_het = [
        # protein 1
        D(m_h_1) ~ R_1*w_max_1*a/(thetax+a) - (lam + dm_h_1)*m_h_1 + gamma/nx*c_h_1 - kb_h_1*r*m_h_1 + ku_h_1*c_h_1
        D(c_h_1) ~ -lam*c_h_1 + kb_h_1*r*m_h_1 - ku_h_1*c_h_1 - gamma/nx*c_h_1
        D(p_h_1) ~ gamma/nx*c_h_1 - (lam + dp_h_1)*p_h_1

        # protein_2
        D(m_h_2) ~ R_2*w_max_2*a/(thetax+a) - (lam + dm_h_2)*m_h_2 + gamma/nx*c_h_2 - kb_h_2*r*m_h_2 + ku_h_2*c_h_2
        D(c_h_2) ~ -lam*c_h_2 + kb_h_2*r*m_h_2 - ku_h_2*c_h_2 - gamma/nx*c_h_2
        D(p_h_2) ~ gamma/nx*c_h_2 - (lam + dp_h_2)*p_h_2

        # protein 3
        D(m_h_3) ~ R_3*w_max_1*a/(thetax+a) - (lam + dm_h_3)*m_h_3 + gamma/nx*c_h_3 - kb_h_3*r*m_h_3 + ku_h_3*c_h_3
        D(c_h_3) ~ -lam*c_h_3 + kb_h_3*r*m_h_3 - ku_h_3*c_h_3 - gamma/nx*c_h_3
        D(p_h_3) ~ gamma/nx*c_h_3 - (lam + dp_h_3)*p_h_3
    ];

    # update energy and ribosomal contributions
    eqs_r_pr_ha = [
        D(r) ~ upd_eqs_dict["ribosomes"][1].rhs - kb_h_1*m_h_1*r + ku_h_1*c_h_1 - kb_h_2*m_h_2*r + ku_h_2*c_h_2 - kb_h_3*m_h_3*r + ku_h_3*c_h_3 + gamma/nx*(c_h_1 + c_h_2 + c_h_3)
    ];
    eqs_a_ha = [
        D(a) ~ upd_eqs_dict["energy"][1].rhs - gamma*(c_h_1 + c_h_2 + c_h_3)
    ];

    # organize host_aware model equations in a dictionary
    eqs_dict = Dict(
        "house_keeping" => upd_eqs_dict["house_keeping"],
        "metabolic"     => upd_eqs_dict["metabolic"],
        "transporter"   => upd_eqs_dict["transporter"],
        "ribo_mc"       => upd_eqs_dict["ribo_mc"],
        "ribosomes"     => eqs_r_pr_ha,
        "heterologous"  => eqs_het,
        "nutrient"      => upd_eqs_dict["nutrient"],
        "energy"        => eqs_a_ha
    )

    param_vals = merge(input_param_vals, het_param_vals)
    ss_vals    = merge(input_ss_vals, het_ss_vals)

    return eqs_dict, lam, param_vals, ss_vals
end;

##################
function update_eqs!(; eqs_dict::Dict, term::Num)
    # create a deep copy of the dictionary
    dict_to_updt = deepcopy(eqs_dict)

    for key in keys(eqs_dict)
        for idx in range(1, length(eqs_dict[key]))
            dict_to_updt[key][idx] = Equation(eqs_dict[key][idx].lhs, eqs_dict[key][idx].rhs + term)
        end
    end

    return dict_to_updt
end;

function compose_model(eqs_dict::Dict)
    D= Differential(t)
    connected = compose(ODESystem(
        [eqs_dict[key][idx] for key in keys(eqs_dict) for idx in range(1, length(eqs_dict[key]))], t; name = :connected));
    return connected
end;

# utility function for dictionary formating and updating
function _format(dict::Dict)
    return vec(collect(values(dict)))
end;

function _update!(dict::Dict, param::String, val::Float64)
    dict[param] = dict[param][1] => val
end;

function update_multiple!(vect::Vector, dict::Dict)
    [Main.GrowthModels._update!(dict, key, val) for (elm, val) in vect for key in keys(dict) if isequal(elm, dict[key][1])];
end;

function map_trajectories(model::ODESystem, sol::ODESolution)
    return Dict(string(states(model)[idx])[1:end-3] => sol[idx,:] for idx in range(1, length(sol[:, end])));
end;

function map_end_vals(model::ODESystem, sol::ODESolution)
    # for each species of a model, returns 
    # the end value -- corresponds to steady state if the system reaches steady state
    return [states(model)[idx] => sol[:, end][idx] for idx in range(1, length(sol[:, end]))];
end;

##################
# function definition for combinatorial assembly
allcombinations_(v...) = vec(collect(Iterators.product(v...)))

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
end;

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

    return format_dict, combinations_list, range_iterator_dict
end;

end