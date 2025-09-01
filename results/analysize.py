import os
import argparse
import pandas as pd
import re
import math


def get_variables(filename):
    result = {}
    with open(filename, 'r') as f:
        current_var_type = ''
        current_decision_var = ''
        current_charger_num = ''
        possible_var_types = ['Hyper-parameters', 'Parameters', 'Decision variables', 'Additional infos', 'Metrics']
        decision_variables = ['X_ijt', 'y_ij', 'SOC_f', 'e_t']
        line = f.readline()
        while line != '':
            is_var_type_def = False
            for var_type in possible_var_types:
                if line.startswith(var_type):
                    result[var_type] = {}
                    current_var_type = var_type
                    is_var_type_def = True
                    current_decision_var = ''
                    break

            if not is_var_type_def:
                if '=' in line:
                    var_name, var_val = line.split('=')
                    var_name = var_name.strip()
                    # If the value contains commas, try to convert to a list
                    if ',' in var_val:
                        try:
                            var_val = list(map(eval, var_val.split(',')))
                        except:
                            var_val = [0.0 for _ in range(len(var_val.split(',')))]
                    else:
                        var_val = eval(var_val)
                    result[current_var_type][var_name] = var_val
                    current_decision_var = ''
                elif ':' in line:
                    var_name, var_val = line.split(':')
                    var_name = var_name.strip()
                    # Two cases: new decision variable or nested decision variable rows
                    if (var_val == '\n') and (current_decision_var == 'X_ijt') and (var_name not in decision_variables):
                        current_charger_num = eval(var_name)
                        result[current_var_type][current_decision_var][current_charger_num] = {}
                    elif (var_val == '\n') and (var_name in decision_variables):
                        current_decision_var = var_name
                        result[current_var_type][current_decision_var] = {}
                        current_charger_num = ''
                    else:
                        if (',' in var_val) or (current_decision_var in decision_variables):
                            if var_val == '\n':
                                var_val = []
                            else:
                                var_val = list(map(eval, var_val.split(',')))
                        else:
                            var_val = eval(var_val)
                        if current_charger_num == '':
                            key = eval(var_name) if var_name.isdigit() else var_name
                            if (current_decision_var == ''):
                                result[current_var_type][key] = var_val
                            else:
                                result[current_var_type][current_decision_var][key] = var_val
                        else:
                            result[current_var_type][current_decision_var][current_charger_num][eval(var_name)] = var_val
            line = f.readline()
    return result


def UTILITY_FCT(soc_d, soc_0, soc_jf):
    return (soc_d - soc_jf)/(soc_d - soc_0)


def metric_jain_fairness_index(result_state):
    params = result_state['Parameters']
    solution_state = result_state['Decision variables']['X_ijt']
    utilities = [0.0 for _ in range(params['n'])]
    for j in range(params['n']):
        sum_max_power = 0.0
        for i in range(params['m']):
            sum_power = 0.0
            for t in range(params['arrivals'][j], params['departures'][j]):
                power = solution_state[i][j][t]
                if power > 0.:
                   l_idx = params[f'w_{i+1}'].index(power)
                   sum_power += power * params[f'eta_{i+1}'][l_idx]
            sum_max_power = max(sum_max_power, sum_power)
        soc_jf = params['soc_0'][j] + (params['tau'] * sum_max_power) / params['battery'][j]
        utilities[j] = UTILITY_FCT(params['soc_d'][j], params['soc_0'][j], soc_jf)
    sum_utilities = sum(utilities)
    sum_utilities_sq = sum(u ** 2 for u in utilities)
    if sum_utilities_sq == 0.0:
        return 0.0
    fairness_index = (sum_utilities ** 2) / (params['n'] * sum_utilities_sq)
    return fairness_index


def metric_energy_shortfall(result_state):
    params = result_state['Parameters']
    solution_state = result_state['Decision variables']['X_ijt']
    total_energy_shortfall = 0.0
    for j in range(params['n']):
        sum_max_power = 0.0
        for i in range(params['m']):
            sum_power = 0.0
            for t in range(params['arrivals'][j], params['departures'][j]):
                power = solution_state[i][j][t]
                if power > 0.0:
                    l_idx = params[f'w_{i+1}'].index(power)
                    sum_power += power * params[f'eta_{i+1}'][l_idx]
            sum_max_power = max(sum_max_power, sum_power)
        soc_jf = params['soc_0'][j] + (params['tau'] * sum_max_power) / params['battery'][j]
        if soc_jf < params['soc_d'][j]:
            total_energy_shortfall += (params['soc_d'][j] - soc_jf) * params['battery'][j]
    return total_energy_shortfall


def metric_gini_coefficient(result_state):
    params = result_state['Parameters']
    solution_state = result_state['Decision variables']['X_ijt']
    utilities = []
    for j in range(params['n']):
        sum_max_power = 0.0
        for i in range(params['m']):
            sum_power = 0.0
            for t in range(params['arrivals'][j], params['departures'][j]):
                power = solution_state[i][j][t]
                if power > 0.0:
                    l_idx = params[f'w_{i+1}'].index(power)
                    sum_power += power * params[f'eta_{i+1}'][l_idx]
            sum_max_power = max(sum_max_power, sum_power)
        soc_jf = params['soc_0'][j] + (params['tau'] * sum_max_power) / params['battery'][j]
        utility_j = UTILITY_FCT(params['soc_d'][j], params['soc_0'][j], soc_jf)
        utilities.append(utility_j)
    sum_d = sum(utilities)
    N = params['n']
    denominator = 2.0 * N * sum_d
    if denominator == 0.0:
        return 0.0
    numerator = 0.0
    for j in range(N):
        for k in range(N):
            numerator += abs(utilities[j] - utilities[k])
    G = numerator / denominator
    return G


def metric_relative_mean_deviation(result_state):
    params = result_state['Parameters']
    solution_state = result_state['Decision variables']['X_ijt']
    u = []
    for j in range(params['n']):
        sum_max_power = 0.0
        for i in range(params['m']):
            sum_power = 0.0
            for t in range(params['arrivals'][j], params['departures'][j]):
                power = solution_state[i][j][t]
                if power > 0.0:
                    l_idx = params[f'w_{i+1}'].index(power)
                    sum_power += power * params[f'eta_{i+1}'][l_idx]
            sum_max_power = max(sum_max_power, sum_power)
        soc_jf = params['soc_0'][j] + (params['tau'] * sum_max_power) / params['battery'][j]
        u_j = UTILITY_FCT(params['soc_d'][j], params['soc_0'][j], soc_jf)
        u.append(u_j)
    sum_u = sum(u)
    n = params['n']
    u_bar = sum_u / n
    if u_bar == 0.0:
        return 0.0
    sum_abs_dev = sum(abs(value - u_bar) for value in u)
    MAD = sum_abs_dev / n
    RMD = MAD / u_bar
    return RMD


def metric_relative_range(result_state):
    """Relative Range: (max(u) - min(u)) / mean(u)"""
    params = result_state['Parameters']
    solution_state = result_state['Decision variables']['X_ijt']
    utilities = []
    for j in range(params['n']):
        sum_max_power = 0.0
        for i in range(params['m']):
            sum_power = 0.0
            for t in range(params['arrivals'][j], params['departures'][j]):
                power = solution_state[i][j][t]
                if power > 0.0:
                    l_idx = params[f'w_{i+1}'].index(power)
                    sum_power += power * params[f'eta_{i+1}'][l_idx]
            sum_max_power = max(sum_max_power, sum_power)
        soc_jf = params['soc_0'][j] + (params['tau'] * sum_max_power) / params['battery'][j]
        utilities.append(UTILITY_FCT(params['soc_d'][j], params['soc_0'][j], soc_jf))
    n = params['n']
    if n == 0:
        return 0.0
    u_bar = sum(utilities) / n
    if u_bar == 0.0:
        return 0.0
    return (max(utilities) - min(utilities)) / u_bar


def metric_coefficient_of_variation(result_state):
    """Coefficient of Variation: std(u) / mean(u)"""
    params = result_state['Parameters']
    solution_state = result_state['Decision variables']['X_ijt']
    utilities = []
    for j in range(params['n']):
        sum_max_power = 0.0
        for i in range(params['m']):
            sum_power = 0.0
            for t in range(params['arrivals'][j], params['departures'][j]):
                power = solution_state[i][j][t]
                if power > 0.0:
                    l_idx = params[f'w_{i+1}'].index(power)
                    sum_power += power * params[f'eta_{i+1}'][l_idx]
            sum_max_power = max(sum_max_power, sum_power)
        soc_jf = params['soc_0'][j] + (params['tau'] * sum_max_power) / params['battery'][j]
        utilities.append(UTILITY_FCT(params['soc_d'][j], params['soc_0'][j], soc_jf))
    n = params['n']
    if n == 0:
        return 0.0
    u_bar = sum(utilities) / n
    if u_bar == 0.0:
        return 0.0
    variance = sum((u - u_bar) ** 2 for u in utilities) / n
    std_dev = math.sqrt(variance)
    return std_dev / u_bar


def metric_hoover_index(result_state):
    """Hoover Index: (1 / (2 n u_bar)) * sum(|u - u_bar|)"""
    params = result_state['Parameters']
    solution_state = result_state['Decision variables']['X_ijt']
    utilities = []
    for j in range(params['n']):
        sum_max_power = 0.0
        for i in range(params['m']):
            sum_power = 0.0
            for t in range(params['arrivals'][j], params['departures'][j]):
                power = solution_state[i][j][t]
                if power > 0.0:
                    l_idx = params[f'w_{i+1}'].index(power)
                    sum_power += power * params[f'eta_{i+1}'][l_idx]
            sum_max_power = max(sum_max_power, sum_power)
        soc_jf = params['soc_0'][j] + (params['tau'] * sum_max_power) / params['battery'][j]
        utilities.append(UTILITY_FCT(params['soc_d'][j], params['soc_0'][j], soc_jf))
    n = params['n']
    if n == 0:
        return 0.0
    u_bar = sum(utilities) / n
    if u_bar == 0.0:
        return 0.0
    total_abs_dev = sum(abs(u - u_bar) for u in utilities)
    return total_abs_dev / (2 * n * u_bar)


def metric_envy_freeness(result_state):
    params = result_state['Parameters']
    solution_state = result_state['Decision variables']['X_ijt']
    utilities = []
    for j in range(params['n']):
        sum_max_power = 0.0
        for i in range(params['m']):
            sum_power = 0.0
            for t in range(params['arrivals'][j], params['departures'][j]):
                power = solution_state[i][j][t]
                if power > 0.0:
                    l_idx = params[f'w_{i+1}'].index(power)
                    sum_power += power * params[f'eta_{i+1}'][l_idx]
            sum_max_power = max(sum_max_power, sum_power)
        soc_jf = params['soc_0'][j] + (params['tau'] * sum_max_power) / params['battery'][j]
        utility_j = UTILITY_FCT(params['soc_d'][j], params['soc_0'][j], soc_jf)
        utilities.append(utility_j)
    total_envy = 0.0
    n = params['n']
    for j in range(n):
        for k in range(n):
            diff = utilities[k] - utilities[j]
            if diff > 0.0:
                total_envy += diff
    return 1/(n*(n-1)) * total_envy


def metric_utility_mean(result_state):
    params = result_state['Parameters']
    solution_state = result_state['Decision variables']['X_ijt']
    utilities = []
    for j in range(params['n']):
        sum_max_power = 0.0
        for i in range(params['m']):
            sum_power = 0.0
            for t in range(params['arrivals'][j], params['departures'][j]):
                power = solution_state[i][j][t]
                if power > 0.0:
                    l_idx = params[f'w_{i+1}'].index(power)
                    sum_power += power * params[f'eta_{i+1}'][l_idx]
            sum_max_power = max(sum_max_power, sum_power)
        soc_jf = params['soc_0'][j] + (params['tau'] * sum_max_power) / params['battery'][j]
        utility_j = UTILITY_FCT(params['soc_d'][j], params['soc_0'][j], soc_jf)
        utilities.append(utility_j)
    return sum(utilities) / params['n']


def flatten_result(result, director, filename):
    """
    Flatten the nested result dictionary into a single-level dictionary.
    Extracts metadata from the filename and includes it in the dictionary.
    """

    VAR_TYPE_FILTER = {'Hyper-parameters', 'Additional infos', 'Metrics'} # Save only this var types
    METRIC_FUNCTIONS = {
      'Envy Freeness': metric_envy_freeness,
      'Relative Mean Deviation': metric_relative_mean_deviation,
      'Gini Coefficient': metric_gini_coefficient,
      'Jain\'s Fairness Index': metric_jain_fairness_index,
      'Energy Shortfall': metric_energy_shortfall,
      'Relative Range': metric_relative_range,
      'Coefficient of Variation': metric_coefficient_of_variation,
      'Hoover Index': metric_hoover_index,
      'Utility': metric_utility_mean,
    }
    
    # Extracting metadata from the filename using regex
    pattern = r"solution_(\d+)_obj=(\d+).*\.log"
    match = re.match(pattern, filename)

    # Default values if the filename does not match the expected pattern
    instance_number, i_init, i_est, i_opt, i_obj = None, None, None, None, None

    if match:
        instance_number, i_obj = map(int, match.groups())

    flat_dict = {
        "Filename": filename,
        "Instance_Number": instance_number,
        "Objective_Index": i_obj
    }

    for var_type, vars_dict in result.items():
        if var_type in VAR_TYPE_FILTER:
            for key, value in vars_dict.items():
                if var_type == 'Metrics':
                    continue
                #    value = METRIC_FUNCTIONS[key](result)
                if isinstance(value, dict):
                    for sub_key, sub_value in value.items():
                        flat_key = f"{var_type}_{key}_{sub_key}"
                        flat_dict[flat_key] = sub_value
                else:
                    flat_key = f"{var_type}_{key}"
                    flat_dict[flat_key] = value

    for metric_name in METRIC_FUNCTIONS:
        flat_key = f"Metrics_{metric_name}"
        flat_dict[flat_key] = METRIC_FUNCTIONS[metric_name](result)

    return flat_dict


def process_directory(directory):
    """
    Process all `.csv` files in the given directory and compile the data into a Pandas DataFrame.
    """
    nb = 0
    for filename in os.listdir(directory):
        if filename.startswith('solution_') and filename.endswith(".log"):
            nb += 1
    
    data = []
    n = 0
    
    # Iterate over all `.txt` files in the directory
    for filename in os.listdir(directory):
        if filename.startswith('solution_') and filename.endswith(".log"):
            n += 1
            file_path = os.path.join(directory, filename)
            result = get_variables(file_path)
            flat_result = flatten_result(result, directory,filename)
            data.append(flat_result)
            print(f'{n}/{nb} files processed!', end='\r')
    
    # Convert list of dictionaries to a Pandas DataFrame
    df = pd.DataFrame(data)
    
    return df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process multiple log files in a directory and save the results as a CSV.")
    parser.add_argument("--dir", type=str, help="Path to the directory containing log files.")

    args = parser.parse_args()
    directory = args.dir

    if not os.path.isdir(directory):
        print(f"Error: The directory '{directory}' does not exist.")
        exit(1)

    # Process files and create DataFrame
    df = process_directory(directory)

    # Save DataFrame to CSV
    output_file = os.path.join(directory, "output_log_analysis.csv")
    df.to_csv(output_file, index=False)

    print(f"Processing complete. Results saved to {output_file}.")



# USAGE: python analysize.py --dir ./solution_h_lns/; mv ./solution_h_lns/output_log_analysis.csv output_h_lns.csv; python analysize.py --dir ./solution_h_lns_mip/; mv ./solution_h_lns_mip/output_log_analysis.csv output_h_lns_mip.csv; python analysize.py --dir ./solution_mip/; mv ./solution_mip/output_log_analysis.csv output_mip.csv; python analysize.py --dir ./solution_h/; mv ./solution_h/output_log_analysis.csv output_h.csv

