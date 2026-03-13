#include "lns.h"



// Earliest Deadline First
void solve_power_allocation(ProblemParams const&params, AnnealingState &state, double const&timeout, bool const&log) {
    assert(state.solution.size() == params.n);
    assert(state.charger.size() == params.n);
    // Initialize the solution with no charging
    for (int j = 0; j < params.n; ++j) {
        assert(state.solution[j].size() == params.d_j[j]-params.r_j[j]);
        for (int t = params.r_j[j]; t < params.d_j[j]; ++t) {
            state.solution[j][t-params.r_j[j]] = -1;
        }
    }
    // Keep track of used capacity at each time step
    vector<double> used_capacity(params.T, 0.0);

    // Compute current_soc
    vector<double> current_soc(params.n, 0.0);
    for(int j = 0; j < params.n; ++j){
        current_soc[j] = params.soc_0[j];
    }
    // Compute energy need of each demand
    vector<double> energy_need(params.n, 0.0);
    for(int j = 0; j < params.n; ++j){
        if(state.charger[j] == -1) continue;
        energy_need[j] = (params.soc_d[j] - current_soc[j]) * params.b_j[j];
    }
    // Compute power_demand
    vector<double> power_demand(params.T, 0.0);
    vector<int> s_j(params.n, 0);
    for(int j = 0; j < params.n; ++j){
        s_j[j] = params.r_j[j];
        if(state.charger[j] == -1) continue;
        for(int j_ = 0; j_ < params.n; ++j_){
            if(j != j_ && state.charger[j_] == state.charger[j] && params.r_j[j_] < params.d_j[j] && params.d_j[j_] > params.r_j[j] && params.d_j[j_] <= params.d_j[j]){
                s_j[j] = max(s_j[j], params.d_j[j_]);
            }
        }
    }
    // Compute energy max available for each demand
    vector<double> energy_max_available(params.n, 0.0);
    for(int j = 0; j < params.n; ++j){
        if(state.charger[j] == -1) continue;
        int i = state.charger[j];
        assert(i >= 0 && i < params.m);
        
        for(int t = s_j[j]; t < params.d_j[j]; ++t){
            assert(t >= 0 && t < params.T);
            assert(params.levels[i]-1 >= -1 && params.levels[i]-1 < (int)params.w_il[i].size() && params.levels[i]-1 < (int)params.eta_il[i].size());
            
            for(int l = params.levels[i] - 1; l >= 0; --l){
                if(power_demand[t] + params.w_il[i][l] <= params.w_G + params.pv_t[t]){
                    energy_max_available[j] += params.w_il[i][l] * params.eta_il[i][l] * params.tau;
                    break;
                }
            }
        }
    }
    assert(s_j.size() >= params.n);
    assert(params.d_j.size() >= params.n);
    assert(energy_need.size() >= params.n);
    assert(energy_max_available.size() >= params.n);
    // Sort the demands in sort_demands by their ratio (energy need / energy max available) in decreasing order
    vector<tuple<int, int, int, double, double, double>> sort_demands;
    for(int j = 0; j < params.n; ++j){
        if(state.charger[j] == -1) continue;
        sort_demands.push_back(make_tuple(j, s_j[j], params.d_j[j], energy_need[j], energy_max_available[j], (energy_max_available[j] > 0.0) ? energy_need[j] / energy_max_available[j] : 0.0));
    }
    sort(sort_demands.begin(), sort_demands.end(), [](const tuple<int, int, int, double, double, double>& a, const tuple<int, int, int, double, double, double>& b) {
        if(get<2>(a) != get<2>(b))
            return get<2>(a) < get<2>(b);
        else
            return get<3>(a) > get<3>(b);
    });
    for (auto &dmd : sort_demands) {
        int j = get<0>(dmd);
        int i = state.charger[j]; // Pre-assigned charger
        if(i == -1) continue;
        assert(i >= 0 && i < params.m);
        for (int t = s_j[j]; t < params.d_j[j]; ++t) {
            assert(t >= 0 && t < params.T);
            assert(t-params.r_j[j] >= 0 && t-params.r_j[j] < (int)state.solution[j].size());
            assert(t < (int)power_demand.size() && t < (int)params.pv_t.size());

            int current_level = state.solution[j][t-params.r_j[j]];
            assert(current_level >= -1 && current_level < params.levels[i] && current_level < (int)params.w_il[i].size() && current_level < (int)params.eta_il[i].size());
            if(current_level == params.levels[i] - 1) continue;
            int new_level = params.levels[i] - 1;
            assert(new_level >= 0 && new_level < params.levels[i] && new_level < (int)params.w_il[i].size() && new_level < (int)params.eta_il[i].size());
            double current_power = 0.0;
            double current_soc_gain = 0.0;
            if(current_level >= 0){
                current_power = params.w_il[i][current_level];
                current_soc_gain = (params.w_il[i][current_level] * params.eta_il[i][current_level] * params.tau) / params.b_j[j];
            }
            while(new_level >= 0 && ((power_demand[t] + params.w_il[i][new_level] - current_power > params.w_G + params.pv_t[t]) || \
                 (current_soc[j] + params.w_il[i][new_level] * params.eta_il[i][new_level] * params.tau / params.b_j[j] - current_soc_gain > params.soc_d[j]))){
                assert(new_level >= 0 && new_level < params.levels[i] && new_level < (int)params.w_il[i].size() && new_level < (int)params.eta_il[i].size());
                new_level--;
            }
            if(new_level >= 0){
                assert(t-params.r_j[j] >= 0 && t-params.r_j[j] < params.d_j[j]-params.r_j[j]);
                state.solution[j][t-params.r_j[j]] = new_level;
                current_soc[j] -= current_soc_gain;
                current_soc[j] += (params.w_il[i][new_level] * params.eta_il[i][new_level] * params.tau / params.b_j[j]);
                power_demand[t] += params.w_il[i][new_level] - current_power;
            }
        }
    }
    // Assert that current_soc if actually the soc of the current demand
    for(int j = 0; j < params.n; ++j) {
        double sumPower = 0.;
        int i = state.charger[j];
        if(i != -1){
            for(int t = params.r_j[j]; t < params.d_j[j]; ++t) {
                int level = state.solution[j][t-params.r_j[j]];
                if(level >= 0)
                    sumPower += params.w_il[i][level] * params.eta_il[i][level];
            }
        }
        assert(std::abs((params.soc_0[j] + (sumPower * params.tau) / params.b_j[j]) - current_soc[j]) < 0.001);
        assert(current_soc[j] >= 0.0);
        assert(current_soc[j] <= params.soc_d[j]);
    }
}

// Helper function to compute objective value given SOC vector
// This is used in the greedy post-optimization of solve_power_allocation_fairness
static double compute_objective_from_soc(const std::vector<double>& soc_jf, const ProblemParams& params) {
    double obj = 0.;
    if(params.objective_index == 0) { // Utilitarian
        for(int j = 0; j < params.n; ++j){
            obj += (params.soc_d[j] - soc_jf[j]) / (params.soc_d[j] - params.soc_0[j]);
        }
    } else if(params.objective_index == 1) { // Quadratic
        for(int j = 0; j < params.n; ++j){
            obj += (((params.soc_d[j] - soc_jf[j]) / (params.soc_d[j] - params.soc_0[j]))*((params.soc_d[j] - soc_jf[j]) / (params.soc_d[j] - params.soc_0[j])));
        }
    } else if(params.objective_index == 2) { // Proportional Fairness
        for(int j = 0; j < params.n; ++j){
            if(params.soc_d[j] - soc_jf[j] > 0)
                obj += log((params.soc_d[j] - soc_jf[j]) / (params.soc_d[j] - params.soc_0[j]));
        }
    } else if(params.objective_index == 3) { // Harmonic Fairness
        double alpha = 2.;
        double exponent = 1. - alpha;
        for(int j = 0; j < params.n; ++j) {
            obj += pow((params.soc_d[j] - soc_jf[j]) / (params.soc_d[j] - params.soc_0[j]), exponent);
        }
        obj = 1. / (1. - alpha) * obj;
    } else if(params.objective_index == 4) { // Max-Min
        for(int j = 0; j < params.n; ++j){
            if((params.soc_d[j] - soc_jf[j]) / (params.soc_d[j] - params.soc_0[j]) > obj)
                obj = (params.soc_d[j] - soc_jf[j]) / (params.soc_d[j] - params.soc_0[j]);
        }
    } else if(params.objective_index == 5) { // Convex_Comb (RMD, Utilitarian)
        const double epsilon = 1e-6;
        int N = params.n;
        vector<double> u(N, 0.0);
        double sum_u = 0.0;
        double sum_s = 0.0;
        for (int j = 0; j < params.n; ++j) {
            double shortfall = (params.soc_d[j] - soc_jf[j]) / (params.soc_d[j] - params.soc_0[j]);
            u[j] = shortfall;
            sum_u += shortfall;
            double rel_shortfall = shortfall;
            sum_s += rel_shortfall;
        }
        double u_bar = sum_u / N;
        if (abs(u_bar) < epsilon) {
            return sum_s / N;
        }
        double sum_abs_deviation = 0.0;
        for (int j = 0; j < N; ++j) {
            sum_abs_deviation += abs(u[j] - u_bar);
        }
        double MAD = sum_abs_deviation / N;
        double RMD = MAD / u_bar;
        double S = sum_s / N;
        const double weight = params.lambda;
        obj = weight * RMD + (1. - weight) * S;
    } else if(params.objective_index == 6) { // Convex_Comb (Gini, Utilitarian)
        const double epsilon = 1e-6;
        int N = params.n;
        vector<double> u;
        u.reserve(N);
        double sum_u = 0.0;
        double sum_rel_shortfall = 0.0;
        for (int j = 0; j < params.n; ++j) {
            double shortfall = (params.soc_d[j] - soc_jf[j]) / (params.soc_d[j] - params.soc_0[j]);
            u.push_back(shortfall);
            sum_u += shortfall;
            double rel_shortfall = shortfall;
            sum_rel_shortfall += rel_shortfall;
        }
        double numerator = 0.0;
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                numerator += abs(u[j] - u[k]);
            }
        }
        double denominator = 2.0 * N * sum_u;
        double G = (denominator > epsilon) ? (numerator / denominator) : 0.0;
        double S = sum_rel_shortfall / N;
        const double weight = params.lambda;
        obj = weight * G + (1. - weight) * S;
    } else if(params.objective_index == 7) { // Convex_Comb (Jain, Utilitarian)
        const double epsilon = 1e-6;
        int N = params.n;
        double sum_u = 0.0;
        double sum_u_sq = 0.0;
        double sum_rel_shortfall = 0.0;
        for (int j = 0; j < params.n; ++j) {
            double shortfall = (params.soc_d[j] - soc_jf[j]) / (params.soc_d[j] - params.soc_0[j]);
            sum_u += shortfall;
            sum_u_sq += shortfall * shortfall;
            double rel_shortfall = shortfall;
            sum_rel_shortfall += rel_shortfall;
        }
        double S = sum_rel_shortfall / N;
        double jain_index = 1.0;
        if (sum_u_sq > epsilon) {
            jain_index = (sum_u * sum_u) / (N * sum_u_sq);
        }
        double unfairness = 1.0 - jain_index;
        const double weight = params.lambda;
        obj = weight * unfairness + (1. - weight) * S;
    } else if(params.objective_index == 8) { // Convex_Comb (Envy, Utilitarian)
        const double epsilon = 1e-6;
        int N = params.n;
        vector<double> u;
        vector<double> s;
        u.reserve(N);
        s.reserve(N);
        double sum_s = 0.0;
        for (int j = 0; j < N; ++j) {
            double shortfall = (params.soc_d[j] - soc_jf[j]) / (params.soc_d[j] - params.soc_0[j]);
            u.push_back(shortfall);
            double rel_shortfall = shortfall;
            s.push_back(rel_shortfall);
            sum_s += rel_shortfall;
        }
        double S = sum_s / N;
        double sum_envy = 0.0;
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                if (j == k) continue;
                sum_envy += max(0.0, s[j] - s[k]);
            }
        }
        double E = (sum_envy / (N * (N - 1)) );
        const double weight = params.lambda;
        obj = weight * E + (1. - weight) * S;
    } else {
        throw std::runtime_error("Invalid objective index");
    }
    return obj;
}

// Water-filling EDF
void solve_power_allocation_fairness(
    ProblemParams const& params,
    AnnealingState      & state,
    double timeout,
    bool   log = false) {
    using clock = std::chrono::high_resolution_clock;
    auto start_time = clock::now();

    const int N = params.n, T = params.T, M = params.m;

    /* ------------------------------------------------------------------
    * 1. initialisation
    * ----------------------------------------------------------------*/
    for (int j = 0; j < N; ++j)
        std::fill(state.solution[j].begin(), state.solution[j].end(), -1);

    std::vector<double> soc (N);          // current SoC (absolute, kWh)
    std::vector<double> need(N);          // remaining energy need  (Wh)
    std::vector<double> allocated(N, 0.0);

    for (int j = 0; j < N; ++j) {
        soc [j] = params.soc_0[j];
        need[j] = (state.charger[j] != -1)
                ? (params.soc_d[j] - soc[j]) * params.b_j[j]
                : 0.0;
    }

    /* ------------------------------------------------------------------
    * 2. serialisation order per charger (EDF by departure)
    * ----------------------------------------------------------------*/
    std::vector<int> s_j(N, 0);                   // earliest eligible slot
    std::vector<std::vector<int>> by_charger(M);
    for (int j = 0; j < N; ++j)
        if (state.charger[j] != -1)
            by_charger[state.charger[j]].push_back(j);

    for (int i = 0; i < M; ++i) {
        auto &L = by_charger[i];
        std::sort(L.begin(), L.end(),
                [&](int a,int b){ return params.d_j[a] < params.d_j[b]; });
        int last_d = 0;
        for (int j : L) {
            s_j[j]  = std::max(params.r_j[j], last_d);
            last_d  = params.d_j[j];
        }
    }

    /* ------------------------------------------------------------------
    * 3. helpers and pre-computed data
    * ----------------------------------------------------------------*/
    std::vector<double> power_demand(T, 0.0);     // Σ grid draw per slot

    /*  minimum non-zero level per charger (used to skip pointless checks) */
    std::vector<double> min_lvl_pow(M, 0.0);
    for (int i = 0; i < M; ++i)
        for (int l = 0; l < params.levels[i]; ++l)
            if (params.w_il[i][l] > 1e-9) { min_lvl_pow[i] = params.w_il[i][l]; break; }

    /* ------------------------------------------------------------------
    * 4. demand-first allocation loop
    *     – iterate over EVs ordered by remaining need (largest first)
    *       and greedily fill their timeline within grid/charger limits
    * ----------------------------------------------------------------*/
    while (true) {
        /* --- timeout guard ------------------------------------------- */
        if (timeout > 0 &&
            std::chrono::duration_cast<std::chrono::duration<double>>(clock::now() - start_time).count() >= timeout) {
            if (log) std::cout << "solve_power_allocation: timeout reached during demand-first loop" << std::endl;
            return;
        }

        /* build list of still-needing demands */
        std::vector<int> order;
        for (int j = 0; j < N; ++j)
            if (state.charger[j] != -1 && need[j] > 1e-6)
                order.push_back(j);
        if (order.empty()) break;                 // all satisfied

        std::sort(order.begin(), order.end(), [&](int a,int b){
            // return params.soc_0[a] < params.soc_0[b];
            // return (params.soc_d[a] - params.soc_0[a])/(params.d_j[a] - s_j[a]) > (params.soc_d[b] - params.soc_0[b])/(params.d_j[b] - s_j[b]);
            return (params.soc_d[a] - params.soc_0[a]) > (params.soc_d[b] - params.soc_0[b]);
            // return (params.soc_d[a] - soc[a]) > (params.soc_d[b] - soc[b]);
        });

        bool any_progress = false;

        for (int j : order) {
            int i = state.charger[j];
            if (i < 0) continue;                 // safety

            for (int t = s_j[j]; t < params.d_j[j] && need[j] > 1e-6; ++t) {
                /* --- timeout guard every 1k iterations -------------- */
                if ((t & 0x3FF) == 0 && timeout > 0 &&
                    std::chrono::duration_cast<std::chrono::duration<double>>(clock::now() - start_time).count() >= timeout) {
                    if (log) std::cout << "solve_power_allocation: timeout reached inside demand-first loop" << std::endl;
                    return;
                }

                /* skip inactive timeslot */
                if (t < params.r_j[j]) continue;

                double room = params.w_G + params.pv_t[t] - power_demand[t];
                if (room < min_lvl_pow[i] - 1e-6) continue;          // no space even for the cheapest level

                /* pick best (highest) feasible charger level */
                int lvl = -1;
                for (int l = params.levels[i] - 1; l >= 0; --l)
                    if (params.w_il[i][l] <= room + 1e-6) { lvl = l; break; }
                if (lvl < 0) continue;

                /* energy delivered by this level */
                double P = params.w_il[i][lvl];
                double E = P * params.eta_il[i][lvl] * params.tau;   // Wh

                /* shrink level if it would overshoot target SoC */
                if (soc[j] + E / params.b_j[j] > params.soc_d[j] + 1e-6) {
                    bool found = false;
                    for (int l = lvl - 1; l >= 0; --l) {
                        P = params.w_il[i][l];
                        E = P * params.eta_il[i][l] * params.tau;
                        if (soc[j] + E / params.b_j[j] <= params.soc_d[j] + 1e-6 && P <= room + 1e-6) {
                            lvl = l; found = true; break;
                        }
                    }
                    if (!found) continue;         // cannot place energy without overshooting
                }

                /* commit allocation */
                int pos = t - params.r_j[j];
                if (state.solution[j][pos] != -1) continue;          // slot already filled (shouldn’t happen)

                state.solution[j][pos] = lvl;
                power_demand[t] += P;
                soc [j]        += E / params.b_j[j];
                need[j]         = std::max(0.0, need[j] - E);
                allocated[j]   += E;
                any_progress    = true;
            }
        }

        if (!any_progress) break;               // grid saturated – stop
    }

    /* ------------------------------------------------------------------
    * 5. sanity check: no EV exceeds target SoC
    * ----------------------------------------------------------------*/
    for (int j = 0; j < N; ++j) {
        int i = state.charger[j];  if (i == -1) continue;
        double sumWh = 0.0;
        for (int t = params.r_j[j]; t < params.d_j[j]; ++t) {
            int l = state.solution[j][t - params.r_j[j]];
            if (l >= 0)
                sumWh += params.w_il[i][l] * params.eta_il[i][l] * params.tau;
        }
        double final_soc = params.soc_0[j] + sumWh / params.b_j[j];
        assert(final_soc <= params.soc_d[j] + 1e-6);
    }

    /* ------------------------------------------------------------------
    * 6. greedy post-optimisation
    * ----------------------------------------------------------------*/
    {
        /* compute initial objective value */
        double best_obj = compute_objective_from_soc(soc, params);

        /* ---- scan EVs, remove slots greedily ------------ */
        for (int j = 0; j < N; ++j) {
            if (timeout > 0 &&
                std::chrono::duration_cast<std::chrono::duration<double>>(clock::now() - start_time).count() >= timeout) {
                if (log) std::cout << "solve_power_allocation: timeout reached inside greedy removal" << std::endl;
                return;
            }
            int i = state.charger[j];            if (i < 0) continue;

            for (int t = s_j[j]; t < params.d_j[j]; ++t) {
                int idx = t - params.r_j[j];
                int lvl = state.solution[j][idx]; if (lvl < 0) continue;

                double E   = params.w_il[i][lvl] * params.eta_il[i][lvl] * params.tau;
                double d_u =  E / params.b_j[j];        // ↑ shortfall

                /* compute objective with this slot removed */
                std::vector<double> test_soc = soc;
                test_soc[j] -= d_u;
                double obj_est = compute_objective_from_soc(test_soc, params);

                if (obj_est + 1e-6 < best_obj) {          // improvement – commit removal
                    state.solution[j][idx] = -1;
                    soc[j]       -=  d_u;
                    allocated[j] -=  E;
                    best_obj = obj_est;
                } else {
                    break;                                  // first miss – stop for this EV
                }
            }
        }
    }
}



void LargeNeighboroodSearch::constructiveOperator(AnnealingState& state, const ProblemParams& params) {
    // If no rejected demands, no insertion is possible
    if (state.rejected_demands.empty()) return;
    bool anyInsertionDone = false;
    double currentObjective = objective(state, params);
    auto op_start = chrono::high_resolution_clock::now();

    // Keep trying until no rejected demands or no insertion found
    while (!state.rejected_demands.empty()) {
        auto now = chrono::high_resolution_clock::now();
        double elapsedOp = chrono::duration_cast<chrono::duration<double>>(now - op_start).count();
        struct DemandInsertInfo {
            int demand;
            double bestImprovement;
            int bestCharger;
            double secondBestImprovement;
            // If no second scenario, secondBestImprovement = bestImprovement or 0
            bool feasible;
            DemandInsertInfo() : demand(-1), bestImprovement(-1e9), bestCharger(-1), secondBestImprovement(-1e9), feasible(false) {}
        };

        // Compute insertion scenarios for all currently rejected demands
        vector<DemandInsertInfo> insertionData;
        insertionData.reserve(state.rejected_demands.size());
        sort(state.rejected_demands.begin(), state.rejected_demands.end(), [&](int j1, int j2) {
            return params.d_j[j1] < params.d_j[j2];
        });
        for (int j : state.rejected_demands) {
            // Check elapsed time and break early if the operator timeout is reached
            now = chrono::high_resolution_clock::now();
            elapsedOp = chrono::duration_cast<chrono::duration<double>>(now - op_start).count();
            if (elapsedOp >= params.timeout_neighbor) {
                cout << "[WARNING] Operator timeout reached during constructiveOperatorRegretRepair" << endl;
                break;
            }
            DemandInsertInfo info;
            info.demand = j;
            auto originalSolution = state.solution;
            int originalCharger = state.charger[j];
            vector<double> improvements;
            vector<int> chargersTried;

            // Try each charger as a scenario
            for (int i = 0; i < params.m; ++i) {
                // Ensure that charger i does not already have a request
                bool conflict = false;
                for (int k = 0; k < params.n; ++k) {
                    // Skip checking the current demand j.
                    if (k != j && state.charger[k] == i && params.d_j[k] == params.d_j[j]) {
                        conflict = true;
                        break;
                    }
                }
                if (conflict)
                    continue;
                state.charger[j] = i;
                if (params.power_allocation_heuristic == 0) {
                    solve_power_allocation(params, state, params.timeout_neighbor, false);
                } else {
                    solve_power_allocation_fairness(params, state, params.timeout_neighbor, false);
                }
                repairingSolution(state, params);
                double newObjective = objective(state, params);
                double improvement = currentObjective - newObjective;
                if(improvement != 0.){
                    // Check if the assignment is not in the tabu list or if it satisfies the aspiration criterion
                    auto tabuEndIter = state.tabu_list.find(state.charger);
                    bool isTabu = false;
                    
                    if (tabuEndIter != state.tabu_list.end()) {
                        if (state.iteration < tabuEndIter->second) {
                            isTabu = true;
                            if (newObjective < state.best_objective_so_far) {
                                isTabu = false;
                            }
                        }
                    }
                    
                    if (!isTabu) {
                        improvements.push_back(improvement);
                        chargersTried.push_back(i);
                    }
                }
                // revert
                state.charger[j] = originalCharger;
                state.solution = originalSolution;
            }

            // Identify best and second best improvements
            if (!improvements.empty()) {
                vector<pair<double,int>> scenarioResults;
                for (size_t idx = 0; idx < improvements.size(); ++idx) {
                    scenarioResults.push_back({improvements[idx], chargersTried[idx]});
                }
                // Sort by improvement descending
                sort(scenarioResults.begin(), scenarioResults.end(),
                          [](auto &a, auto &b){return a.first > b.first;});

                double bestImp = scenarioResults[0].first;
                int bestC = scenarioResults[0].second;
                double secondBestImp = (scenarioResults.size() > 1) ? scenarioResults[1].first : 0.0;

                info.bestImprovement = bestImp;
                info.bestCharger = bestC;
                info.secondBestImprovement = secondBestImp;
                info.feasible = true;
            }
            // else: no improvements means no feasible scenario, info remains not feasible
            insertionData.push_back(info);
        }

        // Select the demand with highest regret
        // regret = bestImprovement - secondBestImprovement if second scenario exists
        // else regret = bestImprovement
        double bestRegret = -1e9;
        int chosenDemand = -1;
        int chosenCharger = -1;
        double chosenBestImp = -1e9;
        bool foundInsertion = false;
        for (auto &info : insertionData) {
            if (!info.feasible || info.demand < 0 || info.bestCharger < 0) continue;
            double regret = info.bestImprovement - info.secondBestImprovement;
            if (info.secondBestImprovement < 0) {
                // Means no second scenario was recorded, use bestImprovement as regret
                regret = info.bestImprovement;
            }
            // Pick the one with the highest regret
            if (regret > bestRegret) {
                bestRegret = regret;
                chosenDemand = info.demand;
                chosenCharger = info.bestCharger;
                chosenBestImp = info.bestImprovement;
                foundInsertion = true;
            }
        }

        if (!foundInsertion) {
            // No improving insertion found: Pick one demand from rejected (e.g. the one with largest 'improvement' even if negative)
            if (!state.rejected_demands.empty()) {
                int fallbackDemand = -1;
                double bestFallback = -1e9;
                int fallbackCharger = -1;
                // Try to force an insertion of some demand even if not improving
                for (auto &info : insertionData) {
                    if (info.demand < 0) continue;
                    // Even if not feasible scenario is found, we try to pick the "least bad" scenario
                    if (!info.feasible) continue; // no scenario at all, skip
                    // We have at least a best scenario
                    if (info.bestImprovement > bestFallback) {
                        bestFallback = info.bestImprovement;
                        fallbackDemand = info.demand;
                        fallbackCharger = info.bestCharger;
                    }
                }
                if (fallbackDemand >= 0 && fallbackCharger >= 0) {
                    // Insert this fallback demand anyway
                    state.rejected_demands.erase(remove(state.rejected_demands.begin(), state.rejected_demands.end(), fallbackDemand), state.rejected_demands.end());
                    state.accepted_demands.push_back(fallbackDemand);
                    state.charger[fallbackDemand] = fallbackCharger;
                    anyInsertionDone = true;
                    // continue to next iteration
                    continue;
                }
            }
            // If we reach here, no fallback insertion found either
            // Break out since we can't insert anymore
            break;
        }
        // Perform the chosen insertion
        int j = chosenDemand;
        state.rejected_demands.erase(remove(state.rejected_demands.begin(), state.rejected_demands.end(), j), state.rejected_demands.end());
        state.accepted_demands.push_back(j);
        state.charger[j] = chosenCharger;
        if (params.power_allocation_heuristic == 0) {
            solve_power_allocation(params, state, params.timeout_neighbor, false);
        } else {
            solve_power_allocation_fairness(params, state, params.timeout_neighbor, false);
        }
        repairingSolution(state, params);
        currentObjective = objective(state, params);
        anyInsertionDone = true;

        if (elapsedOp >= params.timeout_neighbor)
            break;
    }
}

void LargeNeighboroodSearch::destructiveOperator(AnnealingState& state, const ProblemParams& params) {
    // Consider Top-K least served demands
    // If d_j - s_j < S_min
    //   Select the current charger with demand j as a starting point for substring
    // Else
    //   Select a demand k that overlaps with demand j on another charger as a starting point for substring
    //   (use random distribution weighted by s_k - r_k to select k)
    //
    // ==> Instead of removing the whole substring, consider removing only one demand every two demands

    if(state.accepted_demands.empty()) return;
    
    // Compute s_j for each demand (effective start time considering overlaps)
    vector<int> s_j(params.n, 0);
    for(int j = 0; j < params.n; ++j){
        s_j[j] = params.r_j[j];
        if(state.charger[j] == -1) continue;
        for(int j_ = 0; j_ < params.n; ++j_){
            if(j != j_ && state.charger[j_] == state.charger[j] && params.r_j[j_] < params.d_j[j] && params.d_j[j_] > params.r_j[j] && params.d_j[j_] <= params.d_j[j]){
                s_j[j] = max(s_j[j], params.d_j[j_]);
            }
        }
    }

    // Compute SoC difference for each accepted demand
    vector<double> socDiff(params.n, 0.0);
    for (int j : state.accepted_demands) {
        // Compute current SoC for demand j
        double currentSoc = params.soc_0[j];
        int i = state.charger[j];
        assert(i != -1);
        for(int t = s_j[j]; t < params.d_j[j]; ++t){
            int l = state.solution[j][t-params.r_j[j]];
            if(l >= 0)
                currentSoc += (params.w_il[i][l] * params.eta_il[i][l] * params.tau) / params.b_j[j];
        }
        socDiff[j] = params.soc_d[j] - currentSoc;
    }

    // Select a demand with a probability proportional to its SoC difference
    discrete_distribution<int> distribution(socDiff.begin(), socDiff.end());
    int j = distribution(state.rng);
    int charger = state.charger[j];

    // Compute the number of demands on the charger
    int nb_demands_charger = 0;
    for(int j_ = 0; j_ < params.n; ++j_){
        if(state.charger[j_] == charger)
            nb_demands_charger++;
    }

    int maxLen = max(1, (int)ceil(nb_demands_charger / 2));
    double currentObjective = objective(state, params);
    int length = state.lengthControl.adapt(currentObjective);
    length = min(length, maxLen);

    // List all the previous demands on the charger
    vector<int> list_demands;
    for(int j_ = 0; j_ < params.n; ++j_){
        if(state.charger[j_] == charger && params.d_j[j_] < params.d_j[j])
            list_demands.push_back(j_);
    }
    list_demands.push_back(j);
    // Sort list demands by their departure time
    sort(list_demands.begin(), list_demands.end(), [&](int j1, int j2) {
        return params.d_j[j1] < params.d_j[j2];
    });
    int nb_removed = 0;
    // Remove length demands from the solution by beginning from the end
    for(int j_idx = list_demands.size() - 1; j_idx >= 0 && nb_removed < length/2; --j_idx){
        int j_ = list_demands[j_idx];
        // Remove the demand from the solution
        for(int t = params.r_j[j_]; t < params.d_j[j_]; ++t){
            state.solution[j_][t-params.r_j[j_]] = -1;
        }
        state.charger[j_] = -1;
        state.accepted_demands.erase(remove(state.accepted_demands.begin(), state.accepted_demands.end(), j_), state.accepted_demands.end());
        nb_removed++;
    }
    
    // Compute a list of demands k that overlap with demand j on another charger
    vector<double> weights_demands;
    for(int k = 0; k < params.n; ++k){
        if (state.charger[k] != -1 && state.charger[k] != charger &&
                ((params.r_j[k] <= params.r_j[j] && params.r_j[j] <= params.d_j[k]) ||
                (params.r_j[k] <= params.d_j[j] && params.d_j[j] <  params.d_j[k]) ||
                (params.r_j[j] <= params.r_j[k] && params.r_j[k] <= params.d_j[j]) ||
                (params.r_j[j] <= params.d_j[k] && params.d_j[k] <  params.d_j[j]))) {
            weights_demands.push_back(s_j[k] - params.r_j[k]);
        } else {
            weights_demands.push_back(0.0);
        }
    }
    // Select a demand k with a probability proportional to s_k - r_k
    discrete_distribution<int> distribution_k(weights_demands.begin(), weights_demands.end());
    int k = distribution_k(state.rng);
    int charger_k = state.charger[k];
    if(charger_k == -1) return;
    // List all the previous demands on the charger
    list_demands.clear();
    for(int j_ = 0; j_ < params.n; ++j_){
        if(state.charger[j_] == charger_k && params.d_j[j_] < params.d_j[k])
            list_demands.push_back(j_);
    }
    list_demands.push_back(k);
    // Sort list demands by their departure time
    sort(list_demands.begin(), list_demands.end(), [&](int j1, int j2) {
        return params.d_j[j1] < params.d_j[j2];
    });
    nb_removed = 0;
    // Remove length demands from the solution by beginning from the end
    for(int j_idx = list_demands.size() - 1; j_idx >= 0 && nb_removed < length/2; --j_idx){
        int j_ = list_demands[j_idx];
        // Remove the demand from the solution
        for(int t = params.r_j[j_]; t < params.d_j[j_]; ++t){
            state.solution[j_][t-params.r_j[j_]] = -1;
        }
        state.charger[j_] = -1;
        state.accepted_demands.erase(remove(state.accepted_demands.begin(), state.accepted_demands.end(), j_), state.accepted_demands.end());
        nb_removed++;
    }
}


LargeNeighboroodSearch::LargeNeighboroodSearch(double initialTemp, double finalTemp, int trials, int generated, bool save_obj_evol, int seed)
    : initialTemperature(initialTemp), finalTemperature(finalTemp), maxTrials(trials), maxGenerated(generated), save_objective_evolution(save_obj_evol), seed(seed) {}

void LargeNeighboroodSearch::initializeSearch(AnnealingState &state, ProblemParams const &params) {
    // Initialize tabu parameters
    state.tabu_tenure = params.n;
    state.iteration = 0;
    state.best_objective_so_far = std::numeric_limits<double>::max();
    
    // Initialization based on the Interval Graph representation (inspiring from the algorithm proposed by Faigle and Nawjin, Bouzina and Emmons)
    state.charger = std::vector<int>(params.n, -1);
    Solution initial_solution(params.n, vector<int>());
    for(int j = 0; j < params.n; ++j){
        initial_solution[j] = vector<int>(params.d_j[j]-params.r_j[j], -1);
    }
    // Combine the arrival time r_j, departure time d_j, and demand index into a vector of tuples
    vector<tuple<int, int, float, int>> demands; // tuple: (arrival_time, departure_time, demand_index)
    for (int j = 0; j < params.n; ++j) {
        demands.push_back(make_tuple(params.r_j[j], params.d_j[j], params.soc_0[j], j));
    }
    // Sort demands by their arrival time in non-decreasing order
    sort(demands.begin(), demands.end());
    // Initialize an empty set S (it will store the indices of the demands currently assigned to chargers)
    set<int> S;
    // Initialize a vector s_j which is initially the same as r_j
    vector<int> s_j = params.r_j;
    vector<int> next_charger_availability(params.m, -1);
    vector<unordered_set<int>> charger_departures(params.m);

    // Loop through all the demands
    while (!demands.empty()) {
        // Get the next demand j in sorted order of arrival time
        int j = get<3>(demands.front());
        demands.erase(demands.begin());
        // Find if there is a charger available for demand j
        int charger_available = -1;
        for (int i = 0; i < params.m; ++i) {
            if((s_j[j] >= next_charger_availability[i]) && (charger_departures[i].find(params.d_j[j]) == charger_departures[i].end())){
                charger_available = i;
                break;
            }
        }
        // If the charger is available, assign demand j to this charger
        if (charger_available >= 0) {
            state.charger[j] = charger_available; // Assign demand j to charger i
            S.insert(j);               // Add demand j to the set S
            next_charger_availability[charger_available] = params.d_j[j];
            charger_departures[charger_available].insert(params.d_j[j]);
        } else {
            // If no charger was available, we need to reassign a demand
            int max_departure_k = -1;
            int min_departure_k = params.T;
            int k_to_remove = -1;
            int i_to_remove = -1;
            int i_assigned = -1;
            // Find the assigned demand k with the highest departure time in S
            for (int i = 0; i < params.m; ++i) {
                for (int k : S) {
                    if ((state.charger[k] == i) && (params.d_j[k] == next_charger_availability[i]) && (params.d_j[j] < params.d_j[k]) && (params.d_j[k] > max_departure_k)) {
                        max_departure_k = params.d_j[k];
                        k_to_remove = k;
                        i_to_remove = i;
                    } else if((state.charger[k] == i) && (params.d_j[k] == next_charger_availability[i]) && (params.d_j[j] > params.d_j[k]) && (params.d_j[k] < min_departure_k)) {
                        min_departure_k = params.d_j[k];
                    }
                }
            }
            if (min_departure_k != params.T) {
                s_j[j] = min_departure_k;
                demands.push_back(make_tuple(s_j[j], params.d_j[j], params.soc_0[j], j)); // Add k back to the demand queue
                sort(demands.begin(), demands.end()); // Keep demands sorted after reinserting k
            } else if ((min_departure_k == params.T) && (max_departure_k != -1)) {
                // If we found a k to reassign, remove k and assign j
                S.erase(k_to_remove);                  // Remove k from S
                state.charger[k_to_remove] = -1; // Unassign k from the charger
                s_j[k_to_remove] = params.d_j[j];             // Update k's start time to j's departure
                demands.push_back(make_tuple(s_j[k_to_remove], params.d_j[k_to_remove], params.soc_0[k_to_remove], k_to_remove)); // Add k back to the demand queue
                charger_departures[i_to_remove].erase(params.d_j[k_to_remove]);
                sort(demands.begin(), demands.end()); // Keep demands sorted after reinserting k
                // Assign j to this charger
                state.charger[j] = i_to_remove; // Assign j to charger i
                S.insert(j);                          // Add j to S
                next_charger_availability[i_to_remove] = params.d_j[j];
                charger_departures[i_to_remove].insert(params.d_j[j]);
            } else {
                // The demand j is implicitly rejected
                state.charger[j] = -1;
            }
        }
    }
    state.solution = initial_solution;
}

void LargeNeighboroodSearch::searchFunction(AnnealingState& state, const ProblemParams& params) {
    setupDemands(state, params);

    // cout << "destructiveOperator" << endl;
    destructiveOperator(state, params);

    setupDemands(state, params);

    // cout << "constructiveOperator" << endl;
    constructiveOperator(state, params);

    repairingSolution(state, params);

    // cout << "Solve power allocation" << endl;
    if (params.power_allocation_heuristic == 0) {
        solve_power_allocation(params, state, params.timeout_neighbor, false);
    } else {
        solve_power_allocation_fairness(params, state, params.timeout_neighbor, false);
    }

    repairingSolution(state, params);

    // Update tabu list with current solution
    state.iteration++;
    
    // Calculate the objective for the current solution
    double currentObjective = objective(state, params);
    
    // Update best objective if current is better
    if (currentObjective < state.best_objective_so_far) {
        state.best_objective_so_far = currentObjective;
    }
    
    // Add current assignment to tabu list with an expiration
    state.tabu_list[state.charger] = state.iteration + state.tabu_tenure;
    
    // Clean expired entries from tabu list
    auto it = state.tabu_list.begin();
    while (it != state.tabu_list.end()) {
        if (it->second <= state.iteration) {
            it = state.tabu_list.erase(it);
        } else {
            ++it;
        }
    }
}

void LargeNeighboroodSearch::optimize(ProblemParams const &params, AnnealingState &state, double timeout) {
    initializeSearch(state, params);
    if (params.power_allocation_heuristic == 0) {
        solve_power_allocation(params, state, params.timeout_neighbor, false);
    } else {
        solve_power_allocation_fairness(params, state, params.timeout_neighbor, false);
    }
    setupDemands(state, params);
    repairingSolution(state, params);

    Solution currentSolution = state.solution;
    Solution bestSolution = state.solution;
    vector<int> bestCharging = state.charger;

    double bestObjective = objective(state, params);
    state.best_objective_so_far = bestObjective; // Initialize best known objective
    
    if(save_objective_evolution)
        objectiveEvolution.push_back(bestObjective);
    double temperature = initialTemperature;
    default_random_engine generator;
    uniform_real_distribution<double> distribution(0.0, 1.0);

    if(timeout == 0.) {
        cout << "NOTE: Only initialization heuristic was ran (timeout == 0)." << endl;
        return;
    }

    int trial = 0;
    bool accepted = false;
    auto start = chrono::high_resolution_clock::now();

    while (trial < maxTrials){
        int generated = 0;
        accepted = false;

        while (generated < maxGenerated){
            auto now = chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed = chrono::duration_cast<chrono::duration<double>>(now - start);
            if (elapsed.count() >= timeout) {
                cerr << "WARNING: Timeout reached, stopping optimization." << endl;
                state.solution = bestSolution;
                state.charger = bestCharging;
                return; // Exit the function early due to timeout
            }

            searchFunction(state, params);
            double newObj = objective(state, params);
            cout << "Current objective: " << newObj << endl;
            double delta = newObj - bestObjective;
            if(save_objective_evolution)
                objectiveEvolution.push_back(newObj);

            if (delta < 0 || distribution(generator) <= exp(-delta / temperature)) {
                currentSolution = state.solution;
                if (newObj < bestObjective) {
                    accepted = true;
                    // Compute the number of accepted ev demand and print it
                    int acceptedEVCount = 0;
                    for (int j = 0; j < params.n; ++j) {
                        double sumPower = 0.;
                        bool hasCharge = false;
                        int i  = state.charger[j];
                        if(i != -1) {
                            for (int t = params.r_j[j]; t < params.d_j[j]; ++t) {
                                int l_idx = state.solution[j][t-params.r_j[j]];
                                if (l_idx >= 0) {
                                    sumPower += params.w_il[i][l_idx];
                                }
                            }
                            if (sumPower > 0.){
                                hasCharge = true;
                            }
                        }
                        if (hasCharge)
                            acceptedEVCount++;
                    }
                    cout << "New best solution with " << acceptedEVCount << " accepted EV demands - Objective function: " << newObj << endl;
                    //////////////////////////////////////////////////////////////

                    bestSolution = state.solution;
                    bestCharging = state.charger;
                    bestObjective = newObj;
                    state.best_objective_so_far = newObj; // Update best objective
                }
            }
            generated++;
        }
        cout << "\n\n ========= NEW TRIAL ========= \n\n" << endl;
        if(accepted) {
            double r = maxTrials/maxGenerated;
            double b = (temperature - finalTemperature) / (temperature * r * finalTemperature);
            temperature = temperature / (1 + b * temperature);
            trial++;
        }
    }

    state.solution = bestSolution;
    state.charger = bestCharging;
}

void LargeNeighboroodSearch::saveObjectiveEvolutionToFile(const string &filepath) {
    ofstream outFile(filepath);
    // Check if the file was successfully opened
    if (!outFile.is_open()) {
        cerr << "Error: Could not open file " << filepath << " for writing." << endl;
        return;
    }
    // Write each objective function value to the file
    for (const double &objValue : objectiveEvolution) {
        outFile << objValue << endl;
    }
    // Close the file
    outFile.close();
    // Notify that the file was written successfully
    cout << "Objective evolution saved to " << filepath << endl;
}

void repairingSolution(AnnealingState& state, const ProblemParams& params) {
    // Seed with the current time to ensure different outcomes in different runs
    vector<vector<int>> new_x_ijt = state.solution;

    // Calculate w_max for each charger
    vector<double> w_i_max, eta_i_max;
    for(int i = 0; i < params.m; ++i) {
        double max_w_il = 0;
        int idx_max = -1;
        for(int l = 0; l < params.levels[i]; ++l){
            if(params.w_il[i][l] > max_w_il) {
                max_w_il = params.w_il[i][l];
                idx_max = l;
            }
        }
        w_i_max.push_back( max_w_il );
        eta_i_max.push_back( params.eta_il[i][idx_max] );
    }

    // Charger Association
    vector<vector<bool>> y_ij(params.m, vector<bool>(params.n, false));
    vector<int> charger(params.n, -1);
    for(int j = 0; j < params.n; ++j){
        double sumPower = 0.;
        int i = state.charger[j];
        if(i != -1) {
            for(int t = params.r_j[j]; t < params.d_j[j]; ++t){
                int l_idx = new_x_ijt[j][t-params.r_j[j]];
                if(l_idx >= 0) {
                    if(l_idx >= params.levels[i]) {
                        cout << "[WARNING] Charging level " << l_idx << " is >= than the maximum level " << params.levels[i] << " for demand " << j << " and charger " << i << endl;
                        l_idx = params.levels[i] - 1;
                        state.solution[j][t-params.r_j[j]] = l_idx;
                    }
                    sumPower += params.w_il[i][l_idx];
                }
            }
            y_ij[i][j] = (bool)(sumPower > 0.);
            if(y_ij[i][j]){
                charger[j] = i;
            }
        } else {
            for(int t = params.r_j[j]; t < params.d_j[j]; ++t){
                int l_idx = new_x_ijt[j][t-params.r_j[j]];
                if(l_idx >= 0) {
                    cout << "[WARNING] Charging unassigned demand " << j << " charging at time " << t << endl;
                    new_x_ijt[j][t-params.r_j[j]] = -1;
                }
            }
        }
    }

    // Overlap
    for(int i = 0; i < params.m; ++i){
        for(int j = 0; j < params.n; ++j){
            if(y_ij[i][j]){
                for(int k = 0; k < params.n; ++k){
                    if((k != j) && (charger[j] == charger[k]) && (charger[j] != -1) && (params.d_j[k] <= params.d_j[j]) && (params.r_j[j] < params.d_j[k])){
                        for(int t = params.r_j[j]; t < params.d_j[k]; ++t) {
                            if(t >= params.r_j[j] && t < params.d_j[j] && new_x_ijt[j][t-params.r_j[j]] > 0)
                                cout << "[WARNING] Overlap" << endl;
                            new_x_ijt[j][t-params.r_j[j]] = -1;
                        }
                    }
                }
            }
        }
    }

    // Compute s_j for each demand
    vector<int> s_j(params.n, 0);
    for (int j = 0; j < params.n; ++j) {
        s_j[j] = params.r_j[j];
        if (charger[j] == -1) continue;
        for (int j_ = 0; j_ < params.n; ++j_) {
            if (j != j_ && charger[j_] == charger[j] &&
                params.r_j[j_] < params.d_j[j] &&
                params.d_j[j_] > params.r_j[j] &&
                params.d_j[j_] <= params.d_j[j]) {
                s_j[j] = max(s_j[j], params.d_j[j_]);
            }
        }
    }

    for(int j = 0; j < params.n; ++j){
        if(charger[j] != -1){
            for(int t = params.r_j[j]; t < s_j[j]; ++t){
                if(new_x_ijt[j][t-params.r_j[j]] > 0)
                    cout << "[WARNING] Charging inside of [r_j,s_j] for demand " << j << " at time " << t << " below " << s_j[j] << endl;
                new_x_ijt[j][t-params.r_j[j]] = -1;
            }
        }
    }

    // Grid limit
    vector<vector<int>> ev_in_excess_time(params.T);
    vector<double> power_time_excess;
    for(int t = 0; t < params.T; ++t){
        double sumPower = 0.;
        vector<int> list_ev;
        for(int j = 0; j < params.n; ++j){
            int i = state.charger[j];
            if(i == -1) continue;
            if(t >= params.r_j[j] && t < params.d_j[j]){
                int l_idx = new_x_ijt[j][t-params.r_j[j]];
                if(l_idx >= 0) {
                    sumPower += params.w_il[i][l_idx];
                    list_ev.push_back(j);
                }
            }
        }
        if(params.tau * (sumPower - params.pv_t[t]) > params.w_G * params.tau){
            ev_in_excess_time[t] = list_ev;  // Directly assign the list_ev to the appropriate index
        }
        power_time_excess.push_back( params.tau * (sumPower - params.pv_t[t]) - params.w_G * params.tau );
        if(power_time_excess[t] > 0.){
            cout << "[WARNING] grid limit excess:   " << sumPower << "  Limit:  " << params.w_G+params.pv_t[t] << endl;
        }
    }
    for(int t = 0; t < params.T; ++t){
        while(power_time_excess[t] > 0. && !ev_in_excess_time[t].empty()){
            uniform_int_distribution<int> dist(0, ev_in_excess_time[t].size() - 1);
            int random_index = dist(state.rng);
            int j = ev_in_excess_time[t][random_index]; // Select an EV randomly
            for(int i = 0; i < params.m; ++i){
                int power_idx = new_x_ijt[j][t-params.r_j[j]];
                new_x_ijt[j][t-params.r_j[j]] = -1;
                if(power_idx >= 0) {
                    power_time_excess[t] -= params.tau * params.w_il[i][power_idx];
                    ev_in_excess_time[t].erase(ev_in_excess_time[t].begin() + random_index);
                }
            }
        }
    }

    // soc_jf <= soc_jd
    double EPSILON = 0.01;
    for(int j = 0; j < params.n; ++j){
        int i = state.charger[j];
        if(i != -1){
            double sumPower = 0.;
            vector<int> valid_timesteps;
            for(int t = params.r_j[j]; t < params.d_j[j]; ++t){
                int l_idx = new_x_ijt[j][t-params.r_j[j]];
                assert(l_idx >= -1 && l_idx < params.levels[i]);
                if(l_idx >= 0) {
                    sumPower += params.w_il[i][l_idx] * params.eta_il[i][l_idx];
                    valid_timesteps.push_back(t);
                }
            }

            double soc_jf = params.soc_0[j] + (params.tau * sumPower) / params.b_j[j];
            assert(soc_jf <= params.soc_d[j]+EPSILON);
            while(soc_jf > params.soc_d[j]+EPSILON && !valid_timesteps.empty()){
                cout << "[WARNING] soc_jd" << endl;
                uniform_int_distribution<int> dist(0, valid_timesteps.size() - 1); // Select a timstep with power allocated randomly
                int random_index = dist(state.rng);
                int t = valid_timesteps[random_index];
                int l_idx = new_x_ijt[j][t-params.r_j[j]];
                new_x_ijt[j][t-params.r_j[j]] = -1;
                soc_jf -= (params.tau * params.eta_il[i][l_idx] * params.w_il[i][l_idx]) / params.b_j[j];
                valid_timesteps.erase(valid_timesteps.begin() + random_index);
            }
        }
    }

    // Adjusting for soc_jf <= 1 constraint
    for (int j = 0; j < params.n; ++j) {
        int i = state.charger[j];
        if(i != -1){
            double sumPower = 0.0;
            vector<int> chargeTimes;
            // Calculate initial state of charge based on charging
            for (int t = params.r_j[j]; t < params.d_j[j]; ++t) {
                int l_idx = new_x_ijt[j][t-params.r_j[j]];
                if (l_idx >= 0) {
                    sumPower += params.w_il[i][l_idx] * params.eta_il[i][l_idx] * params.tau;
                    chargeTimes.push_back(t);
                }
            }
            double soc_jf = params.soc_0[j] + (sumPower / params.b_j[j]);
            // Reduce charging if soc_jf exceeds 1, starting from the first charging time
            while (soc_jf > 1 && !chargeTimes.empty()) {
                cout << "[WARNING] soc_jf" << endl;
                int t = chargeTimes.front(); // Start reducing from the first charging time
                int l_idx = new_x_ijt[j][t-params.r_j[j]];
                if (l_idx >= 0) {
                    double reducedPower = params.w_il[i][l_idx] * params.eta_il[i][l_idx] * params.tau;
                    soc_jf -= reducedPower / params.b_j[j]; // Adjust soc_jf
                    new_x_ijt[j][t-params.r_j[j]] = -1; // Remove this charging time
                    chargeTimes.erase(chargeTimes.begin()); // Remove this time from vector
                }
            }
        }
    }
    state.solution = new_x_ijt;
}

double LargeNeighboroodSearch::objective(AnnealingState& state, const ProblemParams& params) {
    vector<double> soc_jf(params.n, 0.);
    for(int j = 0; j < params.n; ++j){
        int i = state.charger[j];
        if(i != -1){
            double sumPower = 0.;
            for(int t = params.r_j[j]; t < params.d_j[j]; ++t){
                int l_idx = state.solution[j][t-params.r_j[j]];
                if(l_idx >= 0) {
                    sumPower += params.w_il[i][l_idx] * params.eta_il[i][l_idx];
                }
            }
            soc_jf[j] = params.soc_0[j] + (params.tau * sumPower) / params.b_j[j];
        }
    }

    double obj = 0.;
    if(params.objective_index == 0) { // Utilitarian
        for(int j = 0; j < params.n; ++j){
            obj += (params.soc_d[j] - soc_jf[j]) / (params.soc_d[j] - params.soc_0[j]);
        }
    } else if(params.objective_index == 1) { // Quadratic
        for(int j = 0; j < params.n; ++j){
            obj += (((params.soc_d[j] - soc_jf[j]) / (params.soc_d[j] - params.soc_0[j]))*((params.soc_d[j] - soc_jf[j]) / (params.soc_d[j] - params.soc_0[j])));
        }
    } else if(params.objective_index == 2) { // Proportional Fairness
        for(int j = 0; j < params.n; ++j){
            if(params.soc_d[j] - soc_jf[j] > 0)
                obj += log((params.soc_d[j] - soc_jf[j]) / (params.soc_d[j] - params.soc_0[j]));
        }
    } else if(params.objective_index == 3) { // Harmonic Fairness
        double alpha = 2.;
        double exponent = 1. - alpha;
        for(int j = 0; j < params.n; ++j) {
            obj += pow((params.soc_d[j] - soc_jf[j]) / (params.soc_d[j] - params.soc_0[j]), exponent);
        }
        obj = 1. / (1. - alpha) * obj;
    } else if(params.objective_index == 4) { // Max-Min
        for(int j = 0; j < params.n; ++j){
            if((params.soc_d[j] - soc_jf[j]) / (params.soc_d[j] - params.soc_0[j]) > obj)
                obj = (params.soc_d[j] - soc_jf[j]) / (params.soc_d[j] - params.soc_0[j]);
        }
    } else if(params.objective_index == 5) { // Convex_Comb (RMD, Utilitarian)
        // New objective: convex combination of Relative Mean Deviation and average relative shortfall
        // using soc_jf for the final state-of-charge and remaining_energy for the shortfall.
        // Define: u_j = remaining_energy[j] = soc_d[j] - soc_jf[j] for each EV j.
        // Let:
        //   u_bar = (1/n)*sum_{j} remaining_energy[j],
        //   MAD = (1/n)*sum_{j} |remaining_energy[j] - u_bar|,
        //   RMD = MAD / u_bar,
        //   S   = (1/n)*sum_{j} (remaining_energy[j] / soc_d[j]).
        // The overall objective is: minimize 0.5 * RMD + 0.5 * S.

        const double epsilon = 1e-6; // safeguard for division by zero
        int N = params.n;
        vector<double> u(N, 0.0);  // u[j] = soc_d[j] - soc_jf for each EV j
        double sum_u = 0.0;             // to compute mean shortfall
        double sum_s = 0.0;             // to compute average relative shortfall S
        // Compute soc_jf, u_j and relative shortfall for each EV j
        for (int j = 0; j < params.n; ++j) {
            // Compute shortfall u_j = soc_d[j] - soc_jf (assumes soc_jf <= soc_d[j])
            double shortfall = (params.soc_d[j] - soc_jf[j]) / (params.soc_d[j] - params.soc_0[j]);
            u[j] = shortfall;
            sum_u += shortfall;
            // Compute relative shortfall for j; safeguard against division by zero if soc_d[j] is 0
            double rel_shortfall = shortfall;
            sum_s += rel_shortfall;
        }
        // Compute the mean shortfall u_bar
        double u_bar = sum_u / N;
        // Avoid division by zero in RMD computation
        if (abs(u_bar) < epsilon) {
            // If all EVs are nearly fully charged, return S (which should also be near zero)
            return sum_s / N;
        }
        // Compute Mean Absolute Deviation (MAD) of u values
        double sum_abs_deviation = 0.0;
        for (int j = 0; j < N; ++j) {
            sum_abs_deviation += abs(u[j] - u_bar);
        }
        double MAD = sum_abs_deviation / N;
        double RMD = MAD / u_bar;  // Relative Mean Deviation
        // Compute average relative shortfall S
        double S = sum_s / N;
        // Form the convex combination with equal weights (0.5 each)
        const double weight = params.lambda;
        obj = weight * RMD + (1. - weight) * S;
    } else if(params.objective_index == 6) { // Convex_Comb (Gini, Utilitarian)
        // New objective: convex combination of Gini coefficient and average relative energy shortfall.
        // For each EV j, define u_j = remaining_energy[j] = soc_d[j] - soc_jf[j].
        // Let N = n.
        // Define:
        //    u_bar = (1/N) * sum_{j=0}^{N-1} remaining_energy[j],
        //    and for each pair (j,k), introduce diff[j][k] >= |remaining_energy[j] - remaining_energy[k]|.
        // Then the Gini coefficient is approximated by:
        //    G = (sum_{j,k} diff[j][k]) / (2 * N * u_bar).
        // Also, define the average relative shortfall:
        //    S = (1/N) * sum_{j=0}^{N-1} (remaining_energy[j] / soc_d[j]).
        // The overall objective is: minimize 0.5 * G + 0.5 * S.

        const double epsilon = 1e-6; // safeguard for division by zero
        int N = params.n;
        vector<double> u;  // shortfall for each EV: u[j] = soc_d[j] - soc_jf
        u.reserve(N);
        double sum_u = 0.0;            // to compute the denominator of the Gini coefficient
        double sum_rel_shortfall = 0.0; // to compute the average relative shortfall S
        // For each EV j, compute soc_jf, u[j] and relative shortfall
        for (int j = 0; j < params.n; ++j) {
            // Compute the shortfall: soc_d[j] - soc_jf (assuming soc_jf <= soc_d[j])
            double shortfall = (params.soc_d[j] - soc_jf[j]) / (params.soc_d[j] - params.soc_0[j]);
            u.push_back(shortfall);
            sum_u += shortfall;
            // Compute the relative shortfall for EV j (guard against division by zero)
            double rel_shortfall = shortfall;
            sum_rel_shortfall += rel_shortfall;
        }
        // Compute the Gini coefficient over the shortfalls
        double numerator = 0.0;
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                numerator += abs(u[j] - u[k]);
            }
        }
        double denominator = 2.0 * N * sum_u;
        double G = (denominator > epsilon) ? (numerator / denominator) : 0.0;
        // Compute average relative shortfall S
        double S = sum_rel_shortfall / N;
        // Form the convex combination with equal weights (0.5 each)
        const double weight = params.lambda;
        obj = weight * G + (1. - weight) * S;
    } else if(params.objective_index == 7) { // Convex_Comb (Jain, Utilitarian)
        // New objective: convex combination of Jain's fairness (on shortfalls) and average relative energy shortfall.
        // For each EV j, define u_j = remaining_energy[j] = soc_d[j] - soc_jf[j].
        // Let:
        //    sum_u = sum_{j=0}^{N-1} remaining_energy[j],
        //    sum_u_sq = sum_{j=0}^{N-1} (remaining_energy[j])^2.
        // Jain's fairness index is: jain = (sum_u)^2 / (N * sum_u_sq),
        // and we define unfairness = 1 - jain.
        // Also, S = (1/N) * sum_{j=0}^{N-1} (remaining_energy[j] / soc_d[j]).
        // The overall objective is: minimize 0.5 * (1 - jain) + 0.5 * S.

        const double epsilon = 1e-6; // safeguard against division by zero
        int N = params.n;
        vector<double> u;  // shortfall for each EV: u[j] = soc_d[j] - soc_jf
        double sum_u = 0.0;       // Sum of shortfalls (for Jain's index)
        double sum_u_sq = 0.0;    // Sum of squares of shortfalls (for Jain's index)
        double sum_rel_shortfall = 0.0; // Sum of relative shortfalls (for S)
        // Compute soc_jf, u[j] and relative shortfall for each EV j
        for (int j = 0; j < params.n; ++j) {
            // Compute the shortfall: soc_d[j] - soc_jf (assumes soc_jf <= soc_d[j])
            double shortfall = (params.soc_d[j] - soc_jf[j]) / (params.soc_d[j] - params.soc_0[j]);
            sum_u += shortfall;
            sum_u_sq += shortfall * shortfall;
            // Compute the relative shortfall for EV j (guard against division by zero)
            double rel_shortfall = shortfall;
            sum_rel_shortfall += rel_shortfall;
        }
        // Compute average relative shortfall S
        double S = sum_rel_shortfall / N;
        // Compute Jain's fairness index based on shortfalls
        double jain_index = 1.0;
        if (sum_u_sq > epsilon) {
            jain_index = (sum_u * sum_u) / (N * sum_u_sq);
        }
        // Define an "unfairness" measure so that lower is better
        double unfairness = 1.0 - jain_index;
        // Form the convex combination (with equal weight 0.5 on both terms)
        const double weight = params.lambda;
        obj = weight * unfairness + (1. - weight) * S;
    } else if(params.objective_index == 8) { // Convex_Comb (Envy, Utilitarian)
        // New objective: convex combination of envy-freeness and energy shortfall.
        // For each EV j:
        //    u_j = remaining_energy[j] = soc_d[j] - soc_jf[j]
        //    s[j] = u_j / soc_d[j]
        // For each pair (j,k) with j != k, introduce e[j][k] >= max(0, s[j] - s[k]).
        // Define average envy: E = (1/(N*(N-1)))*sum_{j != k} e[j][k]
        // and average relative shortfall: S = (1/N)*sum_{j} s[j].
        // Overall objective: minimize 0.5 * E + 0.5 * S.

        const double epsilon = 1e-6; // safeguard against division by zero
        int N = params.n;
        vector<double> u; // u[j] = soc_d[j] - soc_jf (shortfall for each EV)
        vector<double> s; // s[j] = relative shortfall = u[j] / soc_d[j]
        u.reserve(N);
        s.reserve(N);
        double sum_s = 0.0; // for computing average relative shortfall S
        // Compute soc_jf, shortfall u[j], and relative shortfall s[j] for each EV j.
        for (int j = 0; j < N; ++j) {
            // Compute shortfall: (assumes soc_jf <= soc_d[j])
            double shortfall = (params.soc_d[j] - soc_jf[j]) / (params.soc_d[j] - params.soc_0[j]);
            u.push_back(shortfall);
            // Compute relative shortfall for EV j; guard against division by zero
            double rel_shortfall = shortfall;
            s.push_back(rel_shortfall);
            sum_s += rel_shortfall;
        }
        // Compute average relative shortfall S.
        double S = sum_s / N;
        // Compute the envy measure E
        // For each pair (j, k), compute e[j,k] = max(0, s[j] - s[k])
        double sum_envy = 0.0;
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                if (j == k) continue;
                sum_envy += max(0.0, s[j] - s[k]);
            }
        }
        double E = (sum_envy / (N * (N - 1)) );
        // Composite objective: convex combination of envy (E) and average relative shortfall (S)
        const double weight = params.lambda;
        obj = weight * E + (1. - weight) * S;
    } else {
        throw std::runtime_error("Invalid objective index");
    }

    return obj;
}

