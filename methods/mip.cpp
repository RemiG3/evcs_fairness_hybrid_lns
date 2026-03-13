#include "mip.h"

template <typename WarnCb>
static std::vector<int> sanitizeWarmStartChargerAssignments(
    ProblemParams const &params,
    AnnealingState const &state,
    WarnCb warn_cb) {
    std::vector<int> chargers = state.charger;
    for (int j = 0; j < params.n; ++j) {
        int i = chargers[j];
        if (i < 0 || i >= params.m) {
            continue;
        }
        double sumPower = 0.0;
        int r = params.r_j[j];
        int d = params.d_j[j];
        const auto &schedule = state.solution[j];
        int window = d - r;
        int max_len = std::min<int>(static_cast<int>(schedule.size()), window);
        for (int offset = 0; offset < max_len; ++offset) {
            int l_idx = schedule[offset];
            if (l_idx >= 0 && l_idx < params.levels[i]) {
                sumPower += params.w_il[i][l_idx];
            }
        }
        if (sumPower <= 0.0) {
            warn_cb(j, i);
            chargers[j] = -1;
        }
    }
    return chargers;
}


class StartCheckCB : public GRBCallback
{
    bool startRejected = false;          // flag we export

  protected:
    void callback() override
    {
        if (where == GRB_CB_MESSAGE) {                       // message hook
            std::string line = getStringInfo(GRB_CB_MSG_STRING);
            if (line.find("User MIP start did not produce a new incumbent solution")
                                                      != std::string::npos) {
                startRejected = true;                       // remember
                abort();                                    // stop Gurobi
            }
        }
    }

  public:
    bool failed() const { return startRejected; }           // public getter
};


class DebugCallback : public GRBCallback {
  std::ofstream &log;
public:
  DebugCallback(std::ofstream &logfile)
    : GRBCallback(), log(logfile) {}

protected:
  void callback() override {
    if (where == GRB_CB_MESSAGE) {
      // Capture every raw log line
      std::string msg = getStringInfo(GRB_CB_MSG_STRING);
      log << msg;
    }
    else if (where == GRB_CB_MIP) {
      // Periodic summary of best incumbent vs. best bound
      double nodecnt = getDoubleInfo(GRB_CB_MIP_NODCNT);
      double objbst  = getDoubleInfo(GRB_CB_MIP_OBJBST);
      double objbnd  = getDoubleInfo(GRB_CB_MIP_OBJBND);
      log << "Node " << nodecnt
          << "  Incumbent=" << objbst
          << "  Bound="     << objbnd << "\n";
    }
  }
};



ModelResult* solve_mip(ProblemParams const &params, AnnealingState const &state, double timeout) {
    try {
        GRBEnv env = GRBEnv(true);
        env.start();
        GRBModel model = GRBModel(env);
        stringstream name;
        auto warn_zero_power = [](int j, int i) {
            cerr << "[WARNING] Warm start: demand " << j
                 << " assigned to charger " << i
                 << " has zero total power; removing assignment for MIP warm start."
                 << endl;
        };
        auto warm_charger = sanitizeWarmStartChargerAssignments(params, state, warn_zero_power);

        ////////////////////////////
        //////// Parameters ////////
        ////////////////////////////

        float epsilon = .00001;
        float obj_coef = 0.;
        vector<vector<bool>> o_jt;
        for(int j = 0; j < params.n; ++j) {
            o_jt.push_back( *new vector<bool>() );
            for(int t = 0; t < params.T; ++t){
                o_jt[j].push_back( (t >= params.r_j[j]) && (t < params.d_j[j]) );
            }
        }

        double w_mean = 0., eta_mean = 0.;
        int nb = 0;
        for(int i = 0; i < params.m; ++i) {
            for(int l = 0; l < params.levels[i]; ++l) {
                w_mean += params.w_il[i][l];
                eta_mean += params.w_il[i][l];
                nb++;
            }
        }
        w_mean = w_mean / nb;
        eta_mean = eta_mean / nb;
        double max_min_power = 0; // Initialize to the lowest possible value
        vector<double> w_i_max, eta_i_max;
        for(int i = 0; i < params.m; ++i) {
            double max_w_il = 0.;
            double min_w_il = std::numeric_limits<double>::max(); // Initialize to maximum double value
            int idx_max = 0;
            for(int l = 0; l < params.levels[i]; ++l){
                if(params.w_il[i][l] > max_w_il) {
                    max_w_il = params.w_il[i][l];
                    idx_max = l;
                }
                if(params.w_il[i][l] < min_w_il) { // Find the minimum power
                    min_w_il = params.w_il[i][l];
                }
            }
            w_i_max.push_back( max_w_il );
            eta_i_max.push_back( params.eta_il[i][idx_max] );
            // Update the maximum of the minimum powers
            if (min_w_il > max_min_power) {
                max_min_power = min_w_il;
            }
        }
        vector<double> dur_min_j;
        vector<double> p_time_j;
        for(int j = 0; j < params.n; ++j){
            p_time_j.push_back( params.d_j[j] - params.r_j[j] );
            dur_min_j.push_back( ((params.soc_d[j] - params.soc_0[j]) * params.b_j[j]) / (w_mean * eta_mean * params.tau) );
        }

        
        ////////////////////////////
        //// Decision Variables ////
        ////////////////////////////

        cout << "Building decision variables" << endl;

        vector<GRBVar> soc_jf(params.n);
        vector<vector<vector<vector<GRBVar>>>> x_ijlt(params.m);
        for(int i = 0; i < params.m; ++i){
            x_ijlt[i] = vector<vector<vector<GRBVar>>>(params.n);
            for(int j = 0; j < params.n; ++j){
                x_ijlt[i][j] = vector<vector<GRBVar>>(params.levels[i]);
                for(int l = 0; l < params.levels[i]; ++l){
                    x_ijlt[i][j][l] = vector<GRBVar>(params.T);
                }
            }
        }
        vector<vector<GRBVar>> y_ij(params.m, vector<GRBVar>(params.n));
        vector<GRBVar> utility(params.n);
        vector<GRBVar> log_utility(params.n);
        vector<GRBVar> pow_utility(params.n);
        vector<GRBVar> log_remaining_energy_relative(params.n);
        
        // Create variables x_ijlt[0][0][0][0], ..., x_ijlt[m-1][n-1][levels[m-1]-1][T-1]
        for (int i = 0; i < params.m; ++i) {
            for (int j = 0; j < params.n; ++j) {
                for (int l = 0; l < params.levels[i]; ++l) {
                    for (int t = 0; t < params.T; t++) {
                        name << "x_ijlt_" << i << "_" << j << "_" << l << "_" << t;
                        if((t >= params.r_j[j]) && (t < params.d_j[j]))
                            x_ijlt[i][j][l][t] = model.addVar(0., 1., obj_coef, GRB_BINARY, name.str());
                        else
                            x_ijlt[i][j][l][t] = model.addVar(0., 0., obj_coef, GRB_BINARY, name.str());
                        name.str("");

                        // Warm start
                        double warm_val = 0.;
                        if((warm_charger[j] == i) && (t >= params.r_j[j]) && (t < params.d_j[j]) && (state.solution[j][t-params.r_j[j]] == l)) {
                            warm_val = 1.;
                        }
                        x_ijlt[i][j][l][t].set(GRB_DoubleAttr_Start, warm_val);
                    }
                }
            }
        }

        // Create variables soc_jf[0], ..., soc_jf[n-1]
        for(int j = 0; j < params.n; j++) {
            name << "soc_" << j << "_f";
            soc_jf[j] = model.addVar(params.soc_0[j], params.soc_d[j], obj_coef, GRB_CONTINUOUS, name.str()); // Constraint n°9
            name.str("");
        }

        // Create variables y_ij[0][0], ..., y_ij[m-1][n-1]
        for(int i = 0; i < params.m; i++) {
            for(int j = 0; j < params.n; j++) {
                name << "y_ij_" << "i" << "_" << j;
                y_ij[i][j] = model.addVar(0., 1., obj_coef, GRB_BINARY, name.str());
                name.str("");

                // Warm start from LNS
                double warm_val = (warm_charger[j] == i) ? 1. : 0.;
                y_ij[i][j].set(GRB_DoubleAttr_Start, warm_val);
            }
        }

        for(int j = 0; j < params.n; j++) {
            name << "utility_" << j;
            utility[j] = model.addVar(0.0, GRB_INFINITY, obj_coef, GRB_CONTINUOUS, name.str());
            name.str("");
            name << "log_utility_" << j;
            log_utility[j] = model.addVar(-GRB_INFINITY, GRB_INFINITY, obj_coef, GRB_CONTINUOUS, name.str());
            name.str("");
        }

        GRBLinExpr max_diff = model.addVar(0., GRB_INFINITY, obj_coef, GRB_CONTINUOUS, "max_diff");

        model.update();

        cout << "Building constraints" << endl;

        /////////////////////
        //// Constraints ////
        /////////////////////

        GRBLinExpr expr, expr1, expr2;

        for(int i = 0; i < params.m; ++i){
            for(int j = 0; j < params.n; ++j){
                bool is_first = true;
                for(int t = 0; t < params.T; ++t){
                    if ((is_first) && (t >= params.r_j[j]) && (t < params.d_j[j])) {
                        expr.clear();
                        for(int j_ = 0; j_ < params.n; ++j_){
                            for(int l_ = 0; l_ < params.levels[i]; ++l_) {
                                for(int t_ = params.r_j[j_]; t_ < min(params.d_j[j_], params.d_j[j]); ++t_) {
                                    assert(t < params.T);
                                    assert(t_ < params.T);
                                    if((j_ != j) && (((t_ >= params.r_j[j]) && (t_ < params.d_j[j]) && (params.d_j[j] <= params.d_j[j_])))){
                                        expr += x_ijlt[i][j_][l_][t_];
                                    }
                                }
                            }
                        }
                        expr1.clear();
                        for(int l = 0; l < params.levels[i]; ++l){
                            expr1 += x_ijlt[i][j][l][t];
                        }
                        model.addConstr( expr1 + 1./(params.n*params.T) * expr <= 1. );
                        is_first = false;
                    } else {
                        expr1.clear();
                        for(int l = 0; l < params.levels[i]; ++l){
                            expr1 += x_ijlt[i][j][l][t];
                        }
                        model.addConstr( expr1 + 1./(params.n*params.T) * expr <= 1. );
                    }
                }
            }

            for(int j = 0; j < params.n; j++){
                expr.clear();
                for(int t = 0; t < params.T; ++t){
                    for(int l = 0; l < params.levels[i]; ++l){
                        expr += x_ijlt[i][j][l][t];
                    }
                }
                model.addConstr(params.T*y_ij[i][j] >= expr);
                model.addConstr(y_ij[i][j] <= expr);
            }
        }

        for(int j = 0; j < params.n; ++j){
            expr.clear();
            for(int i = 0; i < params.m; ++i){
                for(int t = 0; t < params.T; ++t){
                    for(int l = 0; l < params.levels[i]; ++l){
                        expr += x_ijlt[i][j][l][t] * params.w_il[i][l] * params.eta_il[i][l];
                    }
                }
            }
            double epsilon_ = 0.02;
            model.addConstr( soc_jf[j] == (params.soc_0[j] + (params.tau * expr) / params.b_j[j]) );

            expr.clear();
            for(int i = 0; i < params.m; ++i)
                expr += y_ij[i][j];
            model.addConstr( expr <= 1 );

            for(int t = 0; t < params.T; ++t){
                expr.clear();
                for(int i = 0; i < params.m; ++i){
                    for(int l = 0; l < params.levels[i]; ++l){
                        expr += x_ijlt[i][j][l][t];
                    }
                }
                model.addConstr( expr <= 1 );
            }

            model.addConstr( utility[j] == (params.soc_d[j] - soc_jf[j]) / (params.soc_d[j] - params.soc_0[j]) );
            
            // If objective is proportional fairness, add constraint for logarithmic utility
            if(params.objective_index == 2)
                model.addGenConstrLog( utility[j], log_utility[j] );
            
            // If objective is Max-Min, add constraint for max_diff
            if(params.objective_index == 4)
                model.addConstr( max_diff >= utility[j] );
        }


        for(int t = 0; t < params.T; ++t){
            expr.clear();
            for(int j = 0; j < params.n; ++j){
                for(int i = 0; i < params.m; ++i){
                    for(int l = 0; l < params.levels[i]; ++l){
                        expr += x_ijlt[i][j][l][t] * params.w_il[i][l];
                    }
                }
            }
            model.addConstr( params.tau * (expr - params.pv_t[t]) <= params.tau * params.w_G );
        }

        model.update();

        cout << "Building objective function" << endl;


        ////////////////////////////
        //// Objective function ////
        ////////////////////////////

        if(params.objective_index == 0) { // Utilitarian
            GRBLinExpr objective; objective.clear();
            for(int j = 0; j < params.n; ++j)
                objective += utility[j];
            model.setObjective(objective, GRB_MINIMIZE);
        } else if(params.objective_index == 1) { // Quadratic
            GRBQuadExpr objective; objective.clear();
            for(int j = 0; j < params.n; ++j)
                objective += (utility[j]*utility[j]);
            model.setObjective(objective, GRB_MINIMIZE);
        } else if(params.objective_index == 2) { // Proportional Fairness
            GRBLinExpr objective; objective.clear();
            for(int j = 0; j < params.n; ++j)
                objective += log_utility[j];
            model.setObjective(objective, GRB_MINIMIZE);
        } else if(params.objective_index == 3) { // Harmonic Fairness
            double alpha = 2.;
            double exponent = 1.0 - alpha;
            GRBLinExpr objective; objective.clear();
            for(int j = 0; j < params.n; ++j) {
                string pvname = "pow_u_" + to_string(j);
                pow_utility[j] = model.addVar(0., GRB_INFINITY, 0., GRB_CONTINUOUS, pvname);
                model.addGenConstrPow( utility[j], pow_utility[j], exponent, "powu_constr_" + to_string(j) );
                objective += pow_utility[j];
            }
            objective = 1. / (1. - alpha) * objective;
            model.setObjective(objective, GRB_MINIMIZE);
        } else if(params.objective_index == 4) { // Max-Min
            model.setObjective(max_diff, GRB_MINIMIZE);
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

            const double weight = params.lambda;
            const int N = params.n;

            // Create a new variable for the average shortfall: u_bar = (1/N)*sum_{j} remaining_energy[j]
            GRBVar u_bar = model.addVar(0.0, GRB_INFINITY, obj_coef, GRB_CONTINUOUS, "u_bar");

            // For each EV j, remaining_energy[j] is already defined in the model (as soc_d[j]-soc_jf[j]).
            // Define auxiliary variables to linearize |remaining_energy[j] - u_bar|:
            vector<GRBVar> dev(N);
            for (int j = 0; j < N; ++j) {
                dev[j] = model.addVar(0.0, GRB_INFINITY, obj_coef, GRB_CONTINUOUS, "dev_" + to_string(j));
                // Enforce: dev[j] >= remaining_energy[j] - u_bar
                model.addConstr(dev[j] >= utility[j] - u_bar, "dev_pos_" + to_string(j));
                // Enforce: dev[j] >= u_bar - remaining_energy[j]
                model.addConstr(dev[j] >= u_bar - utility[j], "dev_neg_" + to_string(j));
            }

            // Define u_bar = (1/N)*sum_{j} remaining_energy[j]
            GRBLinExpr sum_remaining = 0;
            for (int j = 0; j < N; ++j)
                sum_remaining += utility[j];
            model.addConstr(u_bar * N == sum_remaining, "u_bar_def");

            // Define MAD = (1/N)*sum_{j} dev[j]
            GRBLinExpr MAD = 0;
            for (int j = 0; j < N; ++j)
                MAD += dev[j];
            MAD = MAD / N;

            // Introduce a new variable 'rmd' for the Relative Mean Deviation: RMD = MAD / u_bar.
            // To linearize the ratio as much as possible, we enforce the bilinear equality: rmd * u_bar == MAD.
            GRBVar rmd = model.addVar(0.0, GRB_INFINITY, obj_coef, GRB_CONTINUOUS, "rmd");
            model.addQConstr(rmd * u_bar == MAD, "rmd_def");

            // Define average relative shortfall S = (1/N)*sum_{j} (remaining_energy[j] / soc_d[j]).
            // Since soc_d[j] are parameters, these terms are linear.
            GRBLinExpr S = 0;
            for (int j = 0; j < N; ++j)
                S += utility[j];
            S = S / N;
            model.update();

            // Overall objective: convex combination with equal weights
            //   Objective = 0.5 * RMD + 0.5 * S.
            // Note: rmd is nonlinearly linked via the Q-constr above.
            GRBQuadExpr objective = weight * rmd + (1 - weight) * S;
            model.setObjective(objective, GRB_MINIMIZE);
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

            const double weight = params.lambda;
            const int N = params.n;

            // Sum of remaining_energy: sum_u = sum_{j} remaining_energy[j]
            GRBLinExpr sum_u = 0;
            for (int j = 0; j < N; ++j)
                sum_u += utility[j];

            // Define the average shortfall u_bar = (1/N)*sum_{j} remaining_energy[j]
            GRBVar u_bar = model.addVar(0.0, GRB_INFINITY, obj_coef, GRB_CONTINUOUS, "u_bar");
            model.addConstr(u_bar * N == sum_u, "u_bar_def");

            // For each pair (j,k), introduce an auxiliary variable diff[j][k] >= |remaining_energy[j] - remaining_energy[k]|
            vector<vector<GRBVar>> diff(N, vector<GRBVar>(N));
            for (int j = 0; j < N; ++j) {
                for (int k = 0; k < N; ++k) {
                    stringstream dname;
                    dname << "diff_" << j << "_" << k;
                    diff[j][k] = model.addVar(0.0, GRB_INFINITY, obj_coef, GRB_CONTINUOUS, dname.str());
                    // Linearize the absolute value:
                    model.addConstr(diff[j][k] >= utility[j] - utility[k],
                                    "diff_pos_" + to_string(j) + "_" + to_string(k));
                    model.addConstr(diff[j][k] >= utility[k] - utility[j],
                                    "diff_neg_" + to_string(j) + "_" + to_string(k));
                }
            }

            // Compute the sum of all differences: sum_diff = sum_{j,k} diff[j][k]
            GRBLinExpr sum_diff = 0;
            for (int j = 0; j < N; ++j) {
                for (int k = 0; k < N; ++k)
                    sum_diff += diff[j][k];
            }

            // Introduce a new variable for the Gini coefficient.
            // We impose the bilinear equality: Gini * (2 * N * u_bar) == sum_diff.
            GRBVar Gini = model.addVar(0.0, GRB_INFINITY, obj_coef, GRB_CONTINUOUS, "Gini");
            model.addQConstr(Gini * (2 * N * N * u_bar) == sum_diff, "Gini_def");

            // Define the average relative shortfall:
            // S = (1/N) * sum_{j=0}^{N-1} (remaining_energy[j] / soc_d[j])
            GRBLinExpr S = 0;
            for (int j = 0; j < N; ++j)
                S += utility[j];
            S = S / N;
            model.update();

            // Overall objective: minimize 0.5 * Gini + 0.5 * S.
            GRBQuadExpr objective = weight * Gini + (1 - weight) * S;
            model.setObjective(objective, GRB_MINIMIZE);
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

            const double weight = params.lambda;
            const int N = params.n;  // total number of EVs

            // Compute sum_u and sum_u_sq as (linear/quadratic) expressions in remaining_energy.
            GRBLinExpr sum_u = 0;
            GRBQuadExpr sum_u_sq = 0;
            for (int j = 0; j < N; j++) {
                sum_u += utility[j];
                sum_u_sq += utility[j] * utility[j];
            }

            // Introduce a new variable 'jain' representing Jain's fairness index (in [0,1]).
            GRBVar jain = model.addVar(0.0, 1.0, obj_coef, GRB_CONTINUOUS, "jain");
            GRBVar jain_sum_u_sq = model.addVar(0.0, GRB_INFINITY, obj_coef, GRB_CONTINUOUS, "jain_sum_u_sq");
            GRBVar su = model.addVar(0., GRB_INFINITY, 0.0, GRB_CONTINUOUS, "su");

            model.addConstr(su == sum_u);

            model.addQConstr(jain_sum_u_sq == (N * sum_u_sq), "jain_sum_u_sq_def");

            // Enforce the relationship: jain * (N * sum_u_sq) == (sum_u)^2
            model.addQConstr(jain * jain_sum_u_sq == (su * su), "jain_def");

            // Unfairness is defined as: 1 - jain.
            GRBQuadExpr unfairness = 1. - jain;

            // Define average relative shortfall S = (1/N)*sum_{j=0}^{N-1} (remaining_energy[j] / soc_d[j]).
            GRBLinExpr S = 0;
            for (int j = 0; j < N; j++) {
                S += utility[j];
            }
            S = S / N;
            model.update();

            // Overall objective: minimize 0.5 * (unfairness) + 0.5 * S.
            GRBQuadExpr objective = weight * unfairness + (1 - weight) * S;
            model.setObjective(objective, GRB_MINIMIZE);
        } else if(params.objective_index == 8) { // Convex_Comb (Envy, Utilitarian)
            // New objective: convex combination of Envy-freeness and energy shortfall.
            // For each EV j:
            //    u_j = remaining_energy[j] = soc_d[j] - soc_jf[j]
            //    s[j] = u_j / soc_d[j]
            // For each pair (j,k) with j != k, introduce e[j][k] >= max(0, s[j] - s[k]).
            // Define average envy: E = (1/(N*(N-1)))*sum_{j != k} e[j][k]
            // and average relative shortfall: S = (1/N)*sum_{j} s[j].
            // Overall objective: minimize 0.5 * E + 0.5 * S.

            const int N = params.n;
            const double weight = params.lambda; // equal weighting

            // For each EV j, define s[j] as a linear expression: s[j] = remaining_energy[j] / soc_d[j]
            vector<GRBLinExpr> s_expr(N);
            for (int j = 0; j < N; j++) {
                s_expr[j] = utility[j];
            }

            // For each ordered pair (j,k) with j != k, introduce auxiliary variable e[j][k] to capture max(0, s[j]-s[k])
            vector<vector<GRBVar>> e_vars(N, vector<GRBVar>(N));
            for (int j = 0; j < N; j++) {
                for (int k = 0; k < N; k++) {
                    if (j == k) continue;
                    stringstream ename;
                    ename << "e_" << j << "_" << k;
                    e_vars[j][k] = model.addVar(0., GRB_INFINITY, obj_coef, GRB_CONTINUOUS, ename.str());
                    // Enforce: e[j][k] >= s_expr[j] - s_expr[k] and e[j][k] >= 0 is automatic from the variable's lower bound.
                    model.addConstr(e_vars[j][k] >= s_expr[j] - s_expr[k], "envy_constr_" + to_string(j) + "_" + to_string(k));
                }
            }

            // Compute average envy: E = (1/(N*(N-1))) * sum_{j != k} e_vars[j][k]
            GRBLinExpr sum_envy = 0;
            for (int j = 0; j < N; j++) {
                for (int k = 0; k < N; k++) {
                    if (j == k) continue;
                    sum_envy += e_vars[j][k];
                }
            }
            GRBLinExpr E_expr = sum_envy / (N * (N - 1));

            // Compute average relative shortfall: S = (1/N) * sum_{j} s_expr[j]
            GRBLinExpr S_expr = 0;
            for (int j = 0; j < N; j++) {
                S_expr += s_expr[j];
            }
            S_expr = S_expr / N;
            model.update();

            // Overall objective: convex combination with equal weights: 0.5 * E + 0.5 * S.
            GRBLinExpr objective = weight * E_expr + (1 - weight) * S_expr;
            model.setObjective(objective, GRB_MINIMIZE);
        } else {
            throw std::runtime_error("Invalid objective index: " + to_string(params.objective_index));
        }


        /////////////////
        //// Solving ////
        /////////////////

        cout << "Solving model" << endl;

        model.update();
        model.set(GRB_DoubleParam_TimeLimit, timeout);

        double start_time = model.get(GRB_DoubleAttr_Runtime);
        StartCheckCB cb;
        model.setCallback(&cb);
        try {
          model.optimize();
        } catch (GRBException &e) {
          std::cerr << "Gurobi error " << e.getErrorCode()
                    << ": " << e.getMessage() << "\n";
        }
        double elapsed_time = model.get(GRB_DoubleAttr_Runtime) - start_time;
        int status = model.get(GRB_IntAttr_Status);
        bool solved = (status != GRB_INFEASIBLE);

        if (cb.failed())
            throw std::runtime_error("Gurobi rejected warm-start – aborting run.");

        if (status == GRB_OPTIMAL) {
            cout << "Optimal solution found." << endl;
        } else if (model.get(GRB_IntAttr_SolCount) > 0) {
            cout << "Solution found within the time limit." << endl;
        } else {
            solved = false;
            cerr << "No feasible solution found." << endl;
        }


        //////////////////////////
        //// Display Solution ////
        //////////////////////////

        ModelResult *model_result = (ModelResult*)malloc(sizeof(ModelResult));;
        model_result->elapsed_time = elapsed_time;
        model_result->solved = solved;

        int solCount = model.get(GRB_IntAttr_SolCount);
        double gap = 0.;
        if (status == GRB_OPTIMAL || 
            (status == GRB_TIME_LIMIT && solCount > 0) || 
            (status == GRB_INTERRUPTED && solCount > 0) || 
            status == GRB_SUBOPTIMAL) { // GRB_SUBOPTIMAL is available in Gurobi 11.x
            try {
                gap = model.get(GRB_DoubleAttr_MIPGap);
            } catch (GRBException &e) {
                std::cerr << "Error retrieving MIPGap: " << e.getMessage() << std::endl;
            }
        }

        model_result->gap = gap;
        model_result->obj_val = model.getObjective().getValue();
        model_result->yij_assigned = new vector<vector<float>>(params.m, vector<float>(params.n, 0.));
        model_result->soc_jf_assigned = new vector<double>(params.n);
        model_result->x_ijt = new vector<vector<vector<float>>>(params.m, vector<vector<float>>(params.n, vector<float>(params.T, 0.)));

        if(solved){
            for(int i = 0; i < params.m; ++i) {
                for(int j = 0; j < params.n; ++j) {
                    for(int t = 0; t < params.T; ++t) {
                        float power = 0.;
                        for(int l = 0; l < params.levels[i]; ++l){
                            if(int(round(x_ijlt[i][j][l][t].get(GRB_DoubleAttr_X))) == 1)
                                power = params.w_il[i][l];
                        }
                        (*model_result->x_ijt)[i][j][t] = power;
                    }
                }
            }

            cout << "y_ij: " << endl;
            for(int i = 0; i < params.m; ++i) {
                cout << "i=" << i << ": ";
                for(int j = 0; j < params.n; ++j) {
                    if(int(round(y_ij[i][j].get(GRB_DoubleAttr_X))) == 1) {
                        (*model_result->yij_assigned)[i][j] = 1;
                        cout << j << ", ";
                    }
                }
                cout << endl;
            }
            cout << endl;

            cout << "soc_jf: ";
            for(int j = 0; j < params.n; ++j){
                (*model_result->soc_jf_assigned)[j] = soc_jf[j].get(GRB_DoubleAttr_X);
                cout << (*model_result->soc_jf_assigned)[j];
                if(j < params.n-1)
                    cout << ", ";
            }
            cout << endl;
            

            double soc_sum = 0.0;
            int accepted_count = 0;
            for (int j = 0; j < params.n; ++j) {
                bool is_accepted = false;
                for (int i = 0; i < params.m; ++i) {
                    if (int(round(y_ij[i][j].get(GRB_DoubleAttr_X))) == 1) {
                        is_accepted = true;
                        break;
                    }
                }
                if (is_accepted) {
                    soc_sum += (*model_result->soc_jf_assigned)[j];
                    ++accepted_count;
                }
            }
            if (accepted_count > 0) {
                double soc_mean = soc_sum / accepted_count;
                double variance_sum = 0.0;
                for (int j = 0; j < params.n; ++j) {
                    bool is_accepted = false;
                    for (int i = 0; i < params.m; ++i) {
                        if (int(round(y_ij[i][j].get(GRB_DoubleAttr_X))) == 1) {
                            is_accepted = true;
                            break;
                        }
                    }
                    if (is_accepted) {
                        double diff = (*model_result->soc_jf_assigned)[j] - soc_mean;
                        variance_sum += diff * diff;
                    }
                }
                double standard_deviation = sqrt(variance_sum / accepted_count);
                cout << "Standard Deviation of soc_jf among accepted EVs: " << standard_deviation << endl;
            } else {
                cout << "No EVs accepted. Standard Deviation is not defined." << endl;
            }

        } else {
            cerr << "\tStatus: " << model.get(GRB_IntAttr_Status) << "\n";
        }

        return model_result;

    } catch(const GRBException& e) {
        cerr << "\n\nGUROBI Raised an exception:\n";
        cerr << e.getMessage() << "\n";
        throw;
    }
}
