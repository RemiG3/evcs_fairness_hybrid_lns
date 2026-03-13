#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <gurobi_c++.h>
#include <filesystem>
#include <experimental/filesystem>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <limits>
#include <regex>
#include <algorithm>
#include <cmath>
//#include <direct.h>
#include <stdexcept>
#include <chrono>
#include <random>
#include <unordered_set>
#include <set>
#include <unordered_map>



using namespace std;


template<typename T>
struct FixedVectorHash {
    std::size_t operator()(const std::vector<T>& vec) const {
        // Safely handle vectors of any size
        std::size_t hash = 0;
        // Use the vector's actual size
        std::size_t size = vec.size();
        // If vector is empty, return a default hash
        if (size == 0) {
            return hash;
        }
        // Start with the first element
        hash = std::hash<T>()(vec[0]);
        // Process remaining elements
        for (std::size_t i = 1; i < size; ++i) {
            hash = hash * 31 + std::hash<T>()(vec[i]); // 31 is a common prime multiplier
        }
        return hash;
    }
};

// Forward declaration of ModelState struct
typedef struct {
    bool solved;
    double elapsed_time;
    double gap;
    double obj_val;
    vector<vector<float>> *yij_assigned;
    vector<double> *soc_jf_assigned;
    vector<vector<vector<float>>> *x_ijt;
} ModelResult;



struct AdaptiveLengthControl {
    int minLen, maxLen;
    double currentMaxLen;
    int noImproveIters, stagnationThreshold;
    double lastObjective;

    AdaptiveLengthControl(int minL, int maxL, int stagnThresh)
      : minLen(minL), maxLen(maxL),
        currentMaxLen(minL),
        noImproveIters(0),
        stagnationThreshold(stagnThresh),
        lastObjective(std::numeric_limits<double>::infinity())
    {}

    // adjust the static bounds (e.g. when demList.size() changes)
    void setBounds(int newMin, int newMax) {
        minLen = newMin;
        maxLen = newMax;
        // clamp currentMaxLen into the new range
        if(currentMaxLen > maxLen)
            currentMaxLen = maxLen;
        else if(currentMaxLen < minLen)
            currentMaxLen = minLen;
    }

    // returns the adapted length bound for this iteration
    int adapt(double currentObj) {
        if (currentObj < lastObjective) {
            // improvement: intensify
            currentMaxLen = currentMaxLen / 1.5 < minLen ? minLen : currentMaxLen / 1.5;
            noImproveIters = 0;
            lastObjective = currentObj;
        } else {
            // stagnation
            if (++noImproveIters >= stagnationThreshold) {
                currentMaxLen = currentMaxLen * 1.25 > maxLen ? maxLen : currentMaxLen * 1.25;
                noImproveIters = 0;
            }
        }
        return (int)currentMaxLen;
    }
};


using Solution = vector<vector<int>>; //vector<vector<vector<int>>>;

struct AnnealingState {
    Solution solution;
    vector<int> accepted_demands, rejected_demands, charger;
    std::unordered_set<vector<int>, FixedVectorHash<int>> assignments_known;
    std::unordered_map<vector<int>, int, FixedVectorHash<int>> tabu_list; // Maps assignment to expiration iteration
    int iteration;                  // Current iteration counter
    double best_objective_so_far;   // Best objective value found so far
    int tabu_tenure;                // Number of iterations a solution remains tabu
    mt19937 rng;
    AdaptiveLengthControl lengthControl = AdaptiveLengthControl(8, /*initialMax=*/20, /*stagnThreshold=*/5); 
};

struct ProblemParams {
    int T, m, n;
    double w_G, tau;
    vector<double> soc_0, soc_d, b_j, pv_t, dur_min_j;
    vector<int> r_j, d_j, levels;
    vector<vector<double>> w_il, eta_il;
    vector<vector<bool>> o_jt;
    double timeout_neighbor;
    int objective_index;
    double lambda = 0.5; // weight for convex combination objectives
    int power_allocation_heuristic = 0; // 0 = solve_power_allocation (default), 1 = solve_power_allocation_fairness
    function<double(const AnnealingState&, const ProblemParams&)> objective;
};


// Utility function to get current working directory
string get_current_working_dir();

// Function to trim whitespace from a string
string trim(const string &s);

// Functions to read scenario data from files
void read_solar_production_signal(int scenario, int n, vector<double> &pv);
void read_ev_scenario(int scenario, int &n, vector<int> &r_j, vector<int> &d_j, 
                      vector<double> &soc_0, vector<double> &soc_d, vector<double> &b);
void read_station_config(int scenario, int &m, vector<int> &levels, 
                         vector<vector<double>> &w_il, 
                         vector<vector<double>> &eta_il, 
                         double &w_G, double &tau);

// Functions to convert timesteps when tau differs from default 0.1
void convert_pv_timesteps(vector<double>& pv_t, double tau_old, double tau_new);
void convert_ev_timesteps(vector<int>& r_j, vector<int>& d_j, double tau_old, double tau_new);

// Function to write solution to log file
void write_log_file(ProblemParams params, ModelResult *model_result, string output_filename, double timeout);


void setupDemands(AnnealingState& state, const ProblemParams& params);

ProblemParams buildProblemParams(int nb_ev, int ev_scenario, int pv_scenario, int station_scenario, int objective_index, double timeout_neighbor, double lambda, int power_allocation_heuristic = 0, double tau_override = -1.0);


#endif // UTILS_H 
