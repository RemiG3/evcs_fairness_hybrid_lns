#include "utils.h"


namespace fs = experimental::filesystem;

string get_current_working_dir() {
    string buff = fs::current_path().string();
    string current_working_dir(buff);
    return current_working_dir;
}

// Returns a string with leading and trailing whitespace removed
string trim(const string &s) {
    auto start = s.find_first_not_of(" \t\n\r");
    if (start == string::npos)
        return "";  // String is all whitespace
    auto end = s.find_last_not_of(" \t\n\r");
    return s.substr(start, end - start + 1);
}

void read_solar_production_signal(int scenario, int n, vector<double> &pv){
    string dir_path = get_current_working_dir();
    dir_path = ".";
    string scenario_file_path;
    if(n <= 0)
        scenario_file_path =  dir_path + "/PV-"+ to_string(scenario)+".csv";
    else
        scenario_file_path =  dir_path + "/PV-"+ to_string(n) + "-" + to_string(scenario) + ".csv";
    cout << "Opening file " << scenario_file_path <<  endl;
    fstream fout;
    string line;
    ifstream myFile(scenario_file_path);
    try {
        if (!myFile.is_open())
            throw  runtime_error("Could not open file " + scenario_file_path);
    } catch (const runtime_error &error) {
        cout << error.what()<<  endl;
    }

    string row;
    int idxRow =0;
    while ( getline(myFile, row)) {
        stringstream check1(row);
        string intermediate;
        int idxColumn = 0;

        while ( getline(check1, intermediate, ',')) {
            if (idxRow != 0) {
                switch(idxColumn) {
                    case 1:
                        pv.push_back( stod(intermediate) );
                        break;
                    default:
                        break;
                }

            }
            idxColumn++;
        }
        idxRow++;
    }
    myFile.close();
}

void read_ev_scenario(int scenario, int &n, vector<int> &r_j, vector<int> &d_j, vector<double> &soc_0, vector<double> &soc_d, vector<double> &b){
    string dir_path = get_current_working_dir();
    dir_path  = ".";
    string scenario_file_path;
    if (n <= 0)
        scenario_file_path =  dir_path + "/ev_scenario-"+ to_string(scenario)+".csv";
    else
        scenario_file_path =  dir_path + "/ev_scenario-"+to_string(n)+"-"+ to_string(scenario)+".csv";
    cout << "Opening file " << scenario_file_path <<  endl;
    fstream fout;
    string line;
    ifstream myFile(scenario_file_path);
    try {
        if (!myFile.is_open())
            throw  runtime_error("Could not open file " + scenario_file_path);
    } catch (const runtime_error &error) {
        cout << error.what()<<  endl;
    }

    string row;
    int idxRow = 0;
    vector<string> columnNames;
    while ( getline(myFile, row)) {
        stringstream check1(row);
        string intermediate;
        int idxColumn = 0;

        while ( getline(check1, intermediate, ',')) {
            if (idxRow == 0) {
                columnNames.push_back(intermediate);
            } else {
                if(columnNames[idxColumn] == "arrivals")
                    r_j.push_back( stoi(intermediate) );
                else if(columnNames[idxColumn] == "departures")
                    d_j.push_back( stoi(intermediate) );
                else if(columnNames[idxColumn] == "soc_0")
                    soc_0.push_back( stod(intermediate) );
                else if(columnNames[idxColumn] == "soc_d")
                    soc_d.push_back( stod(intermediate) );
                else if(columnNames[idxColumn] == "battery")
                    b.push_back( stod(intermediate) );
            }
            idxColumn++;
        }
        idxRow++;
    }
    n = --idxRow;
    myFile.close();
}

void read_station_config(int scenario, int &m, vector<int> &levels, vector<vector<double>> &w_il, vector<vector<double>> &eta_il, double &w_G, double &tau) {
    string dir_path = ".";
    string scenario_file_path = dir_path + "/station_" + to_string(scenario) + ".csv";
    ifstream myFile(scenario_file_path);

    try {
        if (!myFile.is_open())
            throw runtime_error("Could not open file " + scenario_file_path);
    } catch (const runtime_error &error) {
        cout << error.what() << endl;
        return;  // Exit the function if file could not be opened
    }

    string row;
    int idxRow = 0;
    vector<string> columnNames;
    tau = -1.;
    w_G = -0.123456789;
    m = -1;
    while (getline(myFile, row)) {
        stringstream check1(row);
        string intermediate;
        int idxColumn = 0;
        vector<double> current_w;  // Vectors start empty
        vector<double> current_eta;

        while (getline(check1, intermediate, ',')) {
            if (idxRow == 0) {
                // First row contains column names
                columnNames.push_back(trim(intermediate));
            } else {
                if (trim(intermediate) != "") {  // Check for non-empty strings
                    if (columnNames[idxColumn] == "tau") {
                        tau = stof(intermediate);
                    } else if (columnNames[idxColumn] == "w_G") {
                        w_G = stod(intermediate);
                    } else if (columnNames[idxColumn] == "m") {
                        m = stoi(intermediate);
                    } else if (columnNames[idxColumn].find("w_") != string::npos) {
                        // Handles w_i columns
                        int index = stoi(columnNames[idxColumn].substr(2)) - 1; // Extract index from name
                        if (index >= current_w.size()) {
                            current_w.resize(index + 1, 0.0);  // Resize and initialize new elements with 0.0
                        }
                        current_w[index] = stod(intermediate);
                    } else if (columnNames[idxColumn].find("eta_") != string::npos) {
                        // Handles eta_i columns
                        int index = stoi(columnNames[idxColumn].substr(4)) - 1; // Extract index from name
                        if (index >= current_eta.size()) {
                            current_eta.resize(index + 1, 0.0);  // Resize and initialize new elements with 0.0
                        }
                        current_eta[index] = stod(intermediate);
                    }
                }
            }
            idxColumn++;
        }
        if (idxRow > 0) {
            // Add the vectors to the main list, ensuring they align by index
            w_il.push_back(current_w);
            eta_il.push_back(current_eta);
        }
        idxRow++;
    }
    if(m == -1)
        m = w_il.size();

    assert(w_G != -0.123456789);
    assert(tau != -1.);
    assert(tau > 0. && tau <= 1.);  // tau must be in ]0,1] (0 excluded, 1 included)
    assert(w_il.size() == eta_il.size());
    assert(eta_il.size() == m);
    for(int i = 0; i< m; ++i) {
        levels.push_back(w_il[i].size());
        assert(w_il[i].size() == eta_il[i].size());
    }
    myFile.close();
}

void convert_pv_timesteps(vector<double>& pv_t, double tau_old, double tau_new) {
    if (tau_old == tau_new) {
        return;  // No conversion needed
    }
    
    size_t old_size = pv_t.size();
    size_t new_size = static_cast<size_t>(round(old_size * tau_old / tau_new));
    
    vector<double> new_pv_t;
    new_pv_t.reserve(new_size);
    
    if (tau_new < tau_old) {
        // Expanding timesteps: use linear interpolation
        for (size_t i = 0; i < new_size; ++i) {
            double old_time = i * tau_new / tau_old;
            size_t idx_low = static_cast<size_t>(floor(old_time));
            size_t idx_high = min(idx_low + 1, old_size - 1);
            double fraction = old_time - idx_low;
            
            if (idx_low >= old_size) {
                new_pv_t.push_back(0.0);  // Extend with zeros if beyond original data
            } else if (idx_low == idx_high || fraction < 1e-9) {
                new_pv_t.push_back(pv_t[idx_low]);
            } else {
                // Linear interpolation
                new_pv_t.push_back(pv_t[idx_low] * (1.0 - fraction) + pv_t[idx_high] * fraction);
            }
        }
    } else {
        // Reducing timesteps: average values over aggregated timesteps
        for (size_t i = 0; i < new_size; ++i) {
            double start_time = i * tau_new / tau_old;
            double end_time = (i + 1) * tau_new / tau_old;
            size_t idx_start = static_cast<size_t>(floor(start_time));
            size_t idx_end = min(static_cast<size_t>(ceil(end_time)), old_size);
            
            if (idx_start >= old_size) {
                new_pv_t.push_back(0.0);
            } else {
                double sum = 0.0;
                double count = 0.0;
                for (size_t j = idx_start; j < idx_end; ++j) {
                    sum += pv_t[j];
                    count += 1.0;
                }
                new_pv_t.push_back(count > 0 ? sum / count : 0.0);
            }
        }
    }
    
    pv_t = new_pv_t;
}

void convert_ev_timesteps(vector<int>& r_j, vector<int>& d_j, double tau_old, double tau_new) {
    if (tau_old == tau_new) {
        return;  // No conversion needed
    }
    
    double conversion_factor = tau_old / tau_new;
    
    for (size_t j = 0; j < r_j.size(); ++j) {
        // Use ceil to ensure time windows are not missed
        r_j[j] = static_cast<int>(ceil(r_j[j] * conversion_factor));
        d_j[j] = static_cast<int>(ceil(d_j[j] * conversion_factor));
        
        // Ensure d_j >= r_j (at least 1 timestep)
        if (d_j[j] <= r_j[j]) {
            d_j[j] = r_j[j] + 1;
        }
    }
}

void write_log_file(ProblemParams params, ModelResult *model_result, string output_filename, double timeout){
    string dir_path = get_current_working_dir();
    string log_filename = output_filename.substr(0, output_filename.find_last_of("."));
    string result_file_path =  dir_path + "/" + log_filename + ".log";
    fs::create_directory(dir_path);

    cout << "Writing file " << result_file_path <<  endl;
    ofstream ofs(result_file_path);

    ofs << "Parameters: \n";
    ofs << "n=" << params.n << "\n";
    ofs << "m=" << params.m << "\n";
    ofs << "T=" << params.T << "\n";
    ofs << "tau=" << params.tau << "\n";
    ofs << "w_G=" << params.w_G << "\n";
    for(unsigned int i = 0; i < params.m; ++i){
        ofs << "w_" + to_string(i+1) + "=";
        for(unsigned int l = 0; l < params.levels[i]; ++l) {
            ofs << params.w_il[i][l];
            if(l < params.levels[i]-1)
                ofs << ",";
        }
        ofs << "\n";
    }
    for(unsigned int i = 0; i < params.m; ++i){
        ofs << "eta_" + to_string(i+1) + "=";
        for(unsigned int l = 0; l < params.levels[i]; ++l) {
            ofs << params.eta_il[i][l];
            if(l < params.levels[i]-1)
                ofs << ",";
        }
        ofs << "\n";
    }

    string list_pv = "pv=";
    for(int t = 0; t < params.T; ++t)
        list_pv += to_string(params.pv_t[t]) + ",";
    list_pv.pop_back();
    ofs << list_pv << "\n";

    string list_soc = "soc_0=";
    for(int j = 0; j < params.n; ++j)
        list_soc += to_string(params.soc_0[j]) + ",";
    list_soc.pop_back();
    ofs << list_soc << "\n";

    string list_soc_d = "soc_d=";
    for(int j = 0; j < params.n; ++j)
        list_soc_d += to_string(params.soc_d[j]) + ",";
    list_soc_d.pop_back();
    ofs << list_soc_d << "\n";

    string list_bmax = "battery=";
    for(int j = 0; j < params.n; ++j)
        list_bmax += to_string(params.b_j[j]) + ",";
    list_bmax.pop_back();
    ofs << list_bmax << "\n";

    ofs << "arrivals=";
    string list_arr = "";
    for(int j = 0; j < params.n; ++j)
        list_arr += to_string(params.r_j[j]) + ",";
    list_arr.pop_back();
    ofs << list_arr << "\n";

    ofs << "departures=";
    string list_dep = "";
    for(int j = 0; j < params.n; ++j)
        list_dep += to_string(params.d_j[j]) + ",";
    list_dep.pop_back();
    ofs << list_dep << "\n";

    ofs << "Decision variables:\n";
    ofs << "X_ijt:\n";
    for(int i = 0; i < params.m; i++){
        ofs << i << ":\n";
        for(int j = 0; j < params.n; j++){
            string list_ev = "";
            ofs << j << ":";
            for(int t = 0; t < params.T; t++){
                list_ev += to_string((*model_result->x_ijt)[i][j][t]) + ",";
            }
            list_ev.pop_back();
            ofs << list_ev << "\n";
        }
    }

    ofs << "y_ij:\n";
    for(int i = 0; i < params.m; i++){
        ofs << i << ":";
        string list_ev = "";
        for(int j = 0; j < params.n; j++){
            if((*model_result->yij_assigned)[i][j] == 1)
                list_ev += to_string(j) + ",";
        }
        if(list_ev.length() > 0)
            list_ev.pop_back();
        ofs << list_ev << "\n";
    }

    string list_socf = "soc_f=";
    for(int j = 0; j < params.n; j++){
        list_socf += to_string((*model_result->soc_jf_assigned)[j]) + ",";
    }
    list_socf.pop_back();
    ofs << list_socf << "\n";

    bool solved = model_result->solved && (model_result->elapsed_time < timeout);
    double total_time = model_result->elapsed_time;

    ofs << "Additional infos: \n";
    ofs << "gap=" << to_string(model_result->gap) << "\n";
    ofs << "computation_time=" << to_string(total_time) << "\n";
    ofs << "obj=" << to_string(model_result->obj_val) << "\n";
    ofs << "solved=" << to_string(solved) << "\n";

    ofs.close();
}


void setupDemands(AnnealingState& state, const ProblemParams& params) {
    state.accepted_demands.clear();
    state.rejected_demands.clear();
    for(int j = 0; j < params.n; ++j){
        int i = state.charger[j];
        if(i == -1) {
            state.rejected_demands.push_back(j);
            for(int level : state.solution[j]){
                assert(level == -1);
            }
            continue;
        }

        double sumPower = 0.;
        for(int level : state.solution[j]){
            if(level != -1){
                sumPower += params.w_il[i][level];
            }
        }
        if(sumPower > 0.){
            state.accepted_demands.push_back(j);
        } else {
            state.rejected_demands.push_back(j);
        }
    }
}

ProblemParams buildProblemParams(int nb_ev, int ev_scenario, int pv_scenario, int station_scenario, int objective_index, double timeout_neighbor, double lambda, int power_allocation_heuristic, double tau_override){
    ProblemParams params;

    // Parameters of the problem
    params.tau = 1.f;
    params.n = nb_ev;
    params.m = 0;
    params.w_G = 0;
    params.objective_index = objective_index;
    params.timeout_neighbor = timeout_neighbor;
    params.lambda = lambda;
    params.power_allocation_heuristic = power_allocation_heuristic;

    read_solar_production_signal(pv_scenario, params.n, params.pv_t);
    read_station_config(station_scenario, params.m, params.levels, params.w_il, params.eta_il, params.w_G, params.tau);
    
    // Store the original tau from config file before potential override
    double tau_from_config = params.tau;
    
    // Override tau if provided via command line
    if (tau_override > 0.0) {
        params.tau = tau_override;
        cout << "Using tau=" << params.tau << " from command line (overriding config file value " << tau_from_config << ")" << endl;
    }
    
    read_ev_scenario(ev_scenario, params.n, params.r_j, params.d_j, params.soc_0, params.soc_d, params.b_j);

    // Convert timesteps if tau differs from the tau used in the data files
    // Data files are assumed to use DEFAULT_TAU = 0.1, but if tau was overridden,
    // we need to convert from the config file tau to the new tau
    const double DEFAULT_TAU = 0.1;
    double tau_source = (tau_override > 0.0) ? tau_from_config : DEFAULT_TAU;
    
    if (abs(params.tau - tau_source) > 1e-6) {
        cout << "Converting timesteps from tau=" << tau_source << " to tau=" << params.tau << endl;
        convert_pv_timesteps(params.pv_t, tau_source, params.tau);
        convert_ev_timesteps(params.r_j, params.d_j, tau_source, params.tau);
    }

    /////////////////////////////////
    // PARAMETERS ORDER PRE-TREATMENT
    std::vector<std::tuple<int, int, double, double, double>> combined;
    for (size_t j = 0; j < params.d_j.size(); ++j) {
        combined.push_back(std::make_tuple(params.d_j[j], params.r_j[j], params.soc_0[j], params.soc_d[j], params.b_j[j]));
    }
    std::sort(combined.begin(), combined.end());
    for (size_t j = 0; j < combined.size(); ++j) {
        params.d_j[j] = std::get<0>(combined[j]);
        params.r_j[j] = std::get<1>(combined[j]);
        params.soc_0[j] = std::get<2>(combined[j]);
        params.soc_d[j] = std::get<3>(combined[j]);
        params.b_j[j] = std::get<4>(combined[j]);
    }
    /////////////////////////////////

    ///////////////////////
    //// Print problem ////
    ///////////////////////

    cout <<"\tEV" << "\tArrival times" << "\tDeparture times" <<"\t\tBattery capacity"  << "\tInitial SOC" << "\t\tDesired SOC" << "\n";
    for(int j=0; j<params.n; j++){
        cout <<"\t" << j << "\t" << setw(sizeof("Arrival times")-1) <<   right <<  params.r_j[j]
            << "\t" << setw(sizeof("Departure times")-1) <<   right <<  params.d_j[j]
            << "\t\t" << setw(sizeof("Battery capacity")-1) <<   right <<  params.b_j[j]
            << "\t" << setw(sizeof("Initial SOC")-1) <<   right << params.soc_0[j]
            << "\t\t" << setw(sizeof("Desired SOC")-1) <<   right << params.soc_d[j] <<"\n";
    }

    cout << "\nNumber of chargers: " << params.m << "\n";
    cout << "Levels of power charging rate delivered by the chargers: " << endl;
    for(int i = 0; i < params.m; ++i) {
        cout << "Charger " << i << ": ";
        for(int l = 0; l < params.levels[i]; ++l) {
            cout << params.w_il[i][l] << " (eff: " << params.eta_il[i][l] << ")";
            if(l < params.levels[i]-1)
                cout << ", ";
        }
        cout << endl;
    }
    cout << "\nGrid limit: " << params.w_G << "\nTau (step time): " << params.tau << "\n";

    cout << "\nSolar production signal: " << endl;
    for (float e : params.pv_t){
        cout << e << ", ";
    }
    cout << "\n" << endl;


    // Find the maximum departure time
    int max_departure = 0;
    for(int j = 0; j < params.n; ++j) {
        max_departure = max(max_departure, params.d_j[j]);
    }
    // Set T to be max_departure + 1 to include the last departure time slot
    params.T = max_departure + 1;

    // Old calculation: params.T = int(round(24*1/tau)); // Number of time slots
    
    // Check if pv and pr arrays need to be extended to match T
    if (params.pv_t.size() < params.T) {
        size_t original_size = params.pv_t.size();
        cout << "Extending solar production signal from size " << original_size << " to " << params.T << " by adding zeros" << endl;
        params.pv_t.resize(params.T, 0.0); // Extend with zeros
    }
    
    assert(params.pv_t.size() == params.T);

    return params;
}
