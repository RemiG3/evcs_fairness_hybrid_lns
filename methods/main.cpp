#include "utils.h"
#include "mip.h"
#include "lns.h"

void print_help() {
    cout << "\n=== EVCS (Electric Vehicle Charging Scheduling) System ===\n\n";
    cout << "USAGE:\n";
    cout << "  ./evcs_optimizer [OPTIONS]\n\n";
    
    cout << "ALGORITHM MODES:\n";
    cout << "  Constructive Heuristic-Only (Fast):\n";
    cout << "    --timeout_lns 0 --timeout_miqp 0\n";
    cout << "  LNS-Only Algorithm (Local Search):\n";
    cout << "    --timeout_lns 600 --timeout_miqp 0\n";
    cout << "  MIQP-Only Algorithm (Exact):\n";
    cout << "    --timeout_lns 0 --timeout_miqp 3600\n";
    cout << "  Hybrid LNS+MIQP (Recommended):\n";
    cout << "    --timeout_lns 600 --timeout_miqp 3000\n\n";
    
    cout << "SCENARIO PARAMETERS:\n";
    cout << "  --scenario N              Overall scenario ID (sets all sub-scenarios)\n";
    cout << "  --ev_scenario N           EV arrival/departure scenario\n";

    cout << "  --pv_scenario N           Solar production scenario\n";
    cout << "  --station_scenario N      Charging station configuration\n";
    cout << "  --nb_ev N                 Number of EVs (overrides scenario file)\n\n";
    
    cout << "OBJECTIVE FUNCTIONS:\n";
    cout << "  --objective_index N       Objective function to optimize:\n";
    cout << "    0 = Utilitarian (minimize total shortfall)\n";
    cout << "    1 = Quadratic (minimize sum of squared shortfalls)\n";
    cout << "    2 = Proportional Fairness (minimize sum of log shortfalls)\n";
    cout << "    3 = Harmonic Fairness (minimize sum of pow shortfalls)\n";
    cout << "    4 = Max-Min (minimize maximum shortfall)\n";
    cout << "    5 = Convex_Comb (RMD, Utilitarian)\n";
    cout << "    6 = Convex_Comb (Gini, Utilitarian)\n";
    cout << "    7 = Convex_Comb (Jain, Utilitarian)\n";
    cout << "    8 = Convex_Comb (Envy, Utilitarian)\n\n";
    cout << "  --lambda L             Weight in [0,1] for convex combination objectives (5-8). Default: 0.5\n\n";
    
    cout << "ALGORITHM PARAMETERS:\n";
    cout << "  --timeout_lns SECONDS     LNS time limit (default: 600)\n";
    cout << "  --timeout_miqp SECONDS    MIQP time limit (default: 3000)\n";
    cout << "  --timeout_neighbor SECS   Neighborhood search time limit (default: 10)\n";
    cout << "  --start_temp TEMP         Initial temperature for annealing (default: 10.0)\n";
    cout << "  --end_temp TEMP           Final temperature (default: 0.00001)\n";
    cout << "  --max_trials N            Maximum trials per temperature (default: 1000)\n";
    cout << "  --max_generated N         Maximum solutions per trial (default: 100)\n";
    cout << "  --seed N                  Random seed (-1 for current time)\n\n";
    
    cout << "OUTPUT CONTROL:\n";
    cout << "  --final_output_filename FILE      Final solution output file\n";
    cout << "  --lns_output_filename FILE        LNS-only solution output file\n";
    cout << "  --objective_evolution_filepath FILE  Objective evolution tracking\n\n";
    
    cout << "DATA REQUIREMENTS:\n";
    cout << "  The following CSV files must be present in the working directory:\n";
    cout << "    ev_scenario-{scenario}.csv    - EV specifications (arrivals, departures, SoC, battery)\n";
    cout << "    station_{scenario}.csv        - Charging station configuration\n";
    cout << "    PV-{scenario}.csv             - Solar production profile\n\n";
    
    cout << "EXAMPLES:\n";
    cout << "  # Constructive heuristic only (fast heuristic)\n";
    cout << "  ./evcs_optimizer --scenario 1 --objective_index 0 --timeout_lns 0 --timeout_miqp 0\n\n";
    cout << "  # LNS optimization (local search)\n";
    cout << "  ./evcs_optimizer --scenario 1 --objective_index 0 --timeout_lns 600 --timeout_miqp 0\n\n";
    cout << "  # Gurobi only (exact optimization)\n";
    cout << "  ./evcs_optimizer --scenario 1 --objective_index 0 --timeout_lns 0 --timeout_miqp 3600\n\n";
    cout << "  # LNS+Gurobi (hybrid optimization)\n";
    cout << "  ./evcs_optimizer --scenario 1 --objective_index 0 --timeout_lns 600 --timeout_miqp 3000\n\n";
    cout << "  # Convex combination objective with custom lambda\n";
    cout << "  ./evcs_optimizer --scenario 1 --objective_index 6 --lambda 0.3 --timeout_lns 600 --timeout_miqp 3000\n\n";
    cout << "  # Custom parameters\n";
    cout << "  ./evcs_optimizer --ev_scenario 1 --pv_scenario 1 \\\n";
    cout << "    --station_scenario 1 --objective_index 6 --seed 42 --start_temp 200.0\n\n";
    
    cout << "For more information, see README.md\n";
}


void run_model(LargeNeighboroodSearch lns, ProblemParams params, string lns_output_filename, string final_output_filename, double timeout_lns, double timeout_miqp){
    try {
        cout << "\n\nMODEL STATE \n" << endl;
        AnnealingState state;
        ModelResult* model_result = nullptr;

        if(lns.seed >= 0) {
            std::mt19937 rng(lns.seed);
            state.rng = rng;
        } else {
            std::mt19937 rng(static_cast<unsigned int>(std::time(nullptr)));
            state.rng = rng;
        }

        auto start_chrono = chrono::high_resolution_clock::now();
        lns.optimize(params, state, timeout_lns);
        auto end_chrono = chrono::high_resolution_clock::now();
        double elapsed_lns = chrono::duration_cast<chrono::duration<double>>(end_chrono - start_chrono).count();

        cout << "BEST OBJECTIVE FOUND: " << lns.objective(state, params) << endl;
        
        // Write LNS-only log if requested
        if (!lns_output_filename.empty()) {
            ModelResult* lns_result = new ModelResult();
            lns_result->solved = true;
            lns_result->elapsed_time = elapsed_lns;
            lns_result->gap = -1.0;
            lns_result->obj_val = lns.objective(state, params);
            lns_result->yij_assigned = new vector<vector<float>>(params.m, vector<float>(params.n, 0.0f));
            for (int j = 0; j < params.n; j++) {
                if (state.charger[j] >= 0) {
                    (*lns_result->yij_assigned)[state.charger[j]][j] = 1.0f;
                }
            }
            lns_result->soc_jf_assigned = new vector<double>(params.n, 0.0);
            for (int j = 0; j < params.n; j++) {
                int i = state.charger[j];
                if (i >= 0) {
                    double sum_power = 0.0;
                    for (int t = params.r_j[j]; t < params.d_j[j]; t++) {
                        int l_idx = state.solution[j][t-params.r_j[j]];
                        if (l_idx >= 0) {
                            sum_power += params.w_il[i][l_idx] * params.eta_il[i][l_idx];
                        }
                    }
                    (*lns_result->soc_jf_assigned)[j] = params.soc_0[j] + (params.tau * sum_power) / params.b_j[j];
                } else {
                    (*lns_result->soc_jf_assigned)[j] = params.soc_0[j];
                }
            }
            lns_result->x_ijt = new vector<vector<vector<float>>>(params.m, vector<vector<float>>(params.n, vector<float>(params.T, 0.0f)));
            for (int j = 0; j < params.n; j++) {
                int i = state.charger[j];
                if (i >= 0) {
                    for (int t = params.r_j[j]; t < params.d_j[j]; t++) {
                        int l_idx = state.solution[j][t-params.r_j[j]];
                        if (l_idx >= 0) {
                            (*lns_result->x_ijt)[i][j][t] = params.w_il[i][l_idx];
                        }
                    }
                }
            }
            write_log_file(params, lns_result, lns_output_filename, elapsed_lns);
        }
        
        // If MIQP timeout is <= 0, skip MIQP solving and use LNS solution
        if (timeout_miqp <= 0) {
            cout << "Skipping MIQP solving as requested (timeout_miqp <= 0)" << endl;
            
            // Create a ModelResult for the LNS solution
            model_result = new ModelResult();
            model_result->solved = true;
            model_result->elapsed_time = elapsed_lns;
            model_result->gap = -1.0; // Unknown gap without MIQP
            model_result->obj_val = lns.objective(state, params);
            
            // Create yij_assigned from charger vector
            model_result->yij_assigned = new vector<vector<float>>(params.m, vector<float>(params.n, 0.0f));
            for (int j = 0; j < params.n; j++) {
                if (state.charger[j] >= 0) {
                    (*model_result->yij_assigned)[state.charger[j]][j] = 1.0f;
                }
            }
            
            // Create soc_jf_assigned vector
            model_result->soc_jf_assigned = new vector<double>(params.n, 0.0);
            for (int j = 0; j < params.n; j++) {
                int i = state.charger[j];
                if (i >= 0) {
                    double sum_power = 0.0;
                    for (int t = params.r_j[j]; t < params.d_j[j]; t++) {
                        int l_idx = state.solution[j][t-params.r_j[j]];
                        if (l_idx >= 0) {
                            sum_power += params.w_il[i][l_idx] * params.eta_il[i][l_idx];
                        }
                    }
                    (*model_result->soc_jf_assigned)[j] = params.soc_0[j] + (params.tau * sum_power) / params.b_j[j];
                } else {
                    (*model_result->soc_jf_assigned)[j] = params.soc_0[j];
                }
            }
            
            // Create x_ijt matrix for the solution
            model_result->x_ijt = new vector<vector<vector<float>>>(
                params.m, vector<vector<float>>(params.n, vector<float>(params.T, 0.0f)));
            
            for (int j = 0; j < params.n; j++) {
                int i = state.charger[j];
                if (i >= 0) {
                    for (int t = params.r_j[j]; t < params.d_j[j]; t++) {
                        int l_idx = state.solution[j][t-params.r_j[j]];
                        if (l_idx >= 0) {
                            (*model_result->x_ijt)[i][j][t] = params.w_il[i][l_idx];
                        }
                    }
                }
            }
        } else {
            model_result = solve_mip(params, state, timeout_miqp);
            model_result->elapsed_time += elapsed_lns;
        }

        if (model_result != nullptr && model_result->solved){

            cout << "\n\nSolution found!\n" << endl;

            ////////////////////
            /// Post-process ///
            ////////////////////
            write_log_file(params, model_result, final_output_filename, timeout_lns+timeout_miqp);
        } else {
            cerr << "\n\nGurobi error!\n";
        }
    } catch(const GRBException e) {
        cerr << "\n\nGUROBI Raised an exception:\n";
        cerr << e.getMessage() << "\n";
        throw;
    }
}


int main(int argc, char* argv[]) {
    int scenario = 0, nb_ev = -1, ev_scenario = -1, pv_scenario = -1, station_scenario = -1;
    double start_temp = 10., end_temp = 0.00001;
    int max_trials = 1000, max_generated = 100, seed = -1;
    double timeout_lns = 600, timeout_miqp = 3000, timeout_neighbor = 10;
    int objective_index = -1;
    double lambda = 0.5;
    double tau = -1.0;  // -1 indicates not set, will use value from config file
    int power_allocation_heuristic = 0;
    string objective_evolution_filepath = "";
    string final_output_filename = "";
    string lns_output_filename = "";

    // Check for help first
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h") {
            print_help();
            return 0;
        }
    }

    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--scenario") {
            if (i+1 < argc) {
                scenario = atoi(argv[++i]);
            }
        }
        if (std::string(argv[i]) == "--ev_scenario") {
            if (i+1 < argc) {
                ev_scenario = atoi(argv[++i]);
            }
        }
        if (std::string(argv[i]) == "--nb_ev") {
            if (i+1 < argc) {
                nb_ev = atoi(argv[++i]);
            }
        }

        if (std::string(argv[i]) == "--pv_scenario") {
            if (i+1 < argc) {
                pv_scenario = atoi(argv[++i]);
            }
        }
        if (std::string(argv[i]) == "--station_scenario") {
            if (i+1 < argc) {
                station_scenario = atoi(argv[++i]);
            }
        }
        if (std::string(argv[i]) == "--objective_index") {
            if (i+1 < argc) {
                objective_index = atoi(argv[++i]);
            }
        }
        if (string(argv[i]) == "--lambda") {
            if (i+1 < argc) {
                lambda = atof(argv[++i]);
                if (lambda < 0.0 || lambda > 1.0) {
                    cerr << "Error: --lambda must be between 0 and 1." << endl;
                    return 1;
                }
            }
        }
        if (string(argv[i]) == "--final_output_filename") {
            final_output_filename = argv[++i];
        }
        if (string(argv[i]) == "--objective_evolution_filepath") {
            objective_evolution_filepath = argv[++i];
        }
        if (string(argv[i]) == "--lns_output_filename") {
            lns_output_filename = argv[++i];
        }
        if (string(argv[i]) == "--start_temp") {
            if (i+1 < argc) {
                start_temp = atof(argv[++i]);
            }
        }
        if (string(argv[i]) == "--end_temp") {
            if (i+1 < argc) {
                end_temp = atof(argv[++i]);
            }
        }
        if (string(argv[i]) == "--max_trials") {
            if (i+1 < argc) {
                max_trials = atoi(argv[++i]);
            }
        }
        if (string(argv[i]) == "--max_generated") {
            if (i+1 < argc) {
                max_generated = atoi(argv[++i]);
            }
        }

        if (string(argv[i]) == "--seed") {
            if (i+1 < argc) {
                seed = atoi(argv[++i]);
            }
        }

        if (string(argv[i]) == "--timeout_lns") {
            if (i+1 < argc) {
                timeout_lns = atof(argv[++i]);
            }
        }
        if (string(argv[i]) == "--timeout_miqp") {
            if (i+1 < argc) {
                timeout_miqp = atof(argv[++i]);
            }
        }
        if (string(argv[i]) == "--timeout_neighbor") {
            if (i+1 < argc) {
                timeout_neighbor = atof(argv[++i]);
            }
        }
        if (string(argv[i]) == "--tau") {
            if (i+1 < argc) {
                tau = atof(argv[++i]);
                if (tau <= 0.0 || tau > 1.0) {
                    cerr << "Error: --tau must be in (0, 1]." << endl;
                    return 1;
                }
            }
        }
        if (string(argv[i]) == "--power_allocation_heuristic") {
            if (i+1 < argc) {
                power_allocation_heuristic = atoi(argv[++i]);
                if (power_allocation_heuristic < 0 || power_allocation_heuristic > 1) {
                    cerr << "Error: --power_allocation_heuristic must be 0 (default) or 1 (fairness)." << endl;
                    return 1;
                }
            }
        }
    }

    if(final_output_filename == "")
        final_output_filename = "complete_solution_" + to_string(scenario) + ".log";
    if (lns_output_filename == "")
        lns_output_filename = "lns_solution_" + to_string(scenario) + ".log";

    if(scenario >= 0){
        ev_scenario = (ev_scenario >= 0) ? ev_scenario : scenario;
        pv_scenario = (pv_scenario >= 0) ? pv_scenario : scenario;
        station_scenario = (station_scenario >= 0) ? station_scenario : scenario;
    }

    ProblemParams params = buildProblemParams(nb_ev, ev_scenario, pv_scenario, station_scenario, objective_index, timeout_neighbor, lambda, power_allocation_heuristic, tau);
    LargeNeighboroodSearch lns(start_temp, end_temp, max_trials, max_generated, (objective_evolution_filepath != ""), seed);
    
    run_model(lns, params, lns_output_filename, final_output_filename, timeout_lns, timeout_miqp);

    return 0;
}

