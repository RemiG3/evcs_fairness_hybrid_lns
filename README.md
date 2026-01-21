# Fairness in EV Charging Scheduling (EVCS): Hybrid LNS + MIP


# 🚧 Release Coming Soon...

## Overview

This repository accompanies the paper **“Fairness-based Optimization  in Electric Vehicle Charging Scheduling  with Heterogeneous Chargers: Comparative Analysis and Solution Approaches”**.
It provides code and experiment assets for a fairness-driven **Electric Vehicle Charging Scheduling (EVCS)** framework that:

- **Formalizes utility and fairness** for EVCS (energy-based utilities and multiple fairness metrics).
- **Implements several objective families**: Utilitarian, Egalitarian (e.g., Proportional & Harmonic Fairness, Max-Min, Quadratic), and **Hybrid welfare** (convex combinations of fairness metrics with total utility).
- **Solves EVCS** via:
  - A **Mixed-Integer** mathematical model (Gurobi).
  - A **Metaheuristic LNS** with destroy/repair operators and a tabu list.
  - A **Constructive Heuristic** with assignment and EDF-based power allocation.
  - An **LNS+Gurobi (warm start)** hybrid optimization for refinement.

> **Scope**: Offline EVCS with known arrivals/departures, heterogeneous chargers (multi-level power), grid power bound, PV supply, fairness objectives, and post-hoc fairness metrics.

**Paper**: **[Fairness-based Optimization in Electric Vehicle Charging Scheduling with Heterogeneous Chargers: Comparative Analysis and Solution Approaches](https://dx.doi.org/10.2139/ssrn.5437728)**

---

## What’s inside

- **Utilities** (energy-based):
  - $U^{(1)}$: Final SoC,
  - $U^{(2)}$: Shortfall,
  - $U^{(3)}$: Proportional shortfall (normalized by need).
- **Fairness metrics**:
  - Extremal: Relative Range (RR)
  - Dispersion: Relative Mean Deviation (RMD), Coefficient of Variation (CV), Jain’s fairness index, Gini coefficient, Hoover index
  - Pairwise: Envy-Freeness (aggregate positive pairwise differences)
  - Efficiency reference: Total utility
- **Objectives**:
  - **Utilitarian** Maximize total utility
  - **Egalitarian**: Max-Min, Quadratic, Proportional (Nash), Harmonic
  - **Hybrid welfare**: Convex combination `λx(utility) − (1−λ)x(fairness)` using RMD, (1−Jain), Gini, or Envy
- **Mathematical model (MIP)** with constraints for assignment, charging windows, grid limit, and SoC evolution.
- **Algorithms**:
  - **Constructive heuristic**:
    - Assignment inspired by interval scheduling with waiting,
    - EDF-style power allocation that respects grid/PV and desired SoC.
  - **LNS**:
    - Destroy: **SISR-based** string removal on overlapping requests (+ adaptive removal length),
    - Repair: **regret-based insertion** with tabu on EV-charger assignments,
    - **Metropolis acceptance** and cooling schedule for diversification/intensification.

---

## Build

**Requirements**
- **C++17**, **CMake ≥ 3.5**
- **Gurobi 11.x** (with C++ API configured)
- **Python 3.10+** (for analysis/plots), with `pandas`, `matplotlib`, `numpy`

**Configure & compile**

The build system automatically detects Gurobi installation. You have several options:

**Option 1: Automatic Detection (Recommended)**
```bash
cd methods
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j
```
This works if Gurobi is installed in standard locations:
- Linux: `/opt/gurobi*/linux64` or `~/gurobi*/linux64`
- macOS: `/Library/gurobi*/macos_universal2`
- Windows: `C:/gurobi*/win64`

**Option 2: Environment Variable**
```bash
export GUROBI_HOME=/path/to/your/gurobi/installation
cd methods
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j
```

**Option 3: CMake Parameter**
```bash
cd methods
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DGUROBI_HOME=/path/to/your/gurobi/installation ..
make -j
```

**For persistent setup (recommended):**
Add to your `~/.bashrc` or `~/.zshrc`:
```bash
export GUROBI_HOME=/path/to/your/gurobi/installation
```

**Getting Help**

The system includes comprehensive built-in help documentation:

```bash
cd methods/build
./evcs_optimizer --help    # Full help documentation
./evcs_optimizer -h        # Same as --help
```

**Usage Examples**

After successful compilation, the executable `evcs_optimizer` will be available in your build directory. The system supports four main algorithm modes:

### Algorithm Modes

**1. Constructive Heuristic-Only (Fast Heuristic)**
```bash
cd methods/build
./evcs_optimizer --scenario 1 --objective_index 0 --timeout_lns 0 --timeout_miqp 0
```
- Uses only the constructive heuristic with EDF-based power allocation
- Fastest execution, ideal for quick feasible solutions
- Good baseline for comparison with other methods
- Set both timeouts to `0` to use only the initialization heuristic

**2. LNS-Only Algorithm (Local Search)**
```bash
cd methods/build
./evcs_optimizer --scenario 1 --objective_index 0 --timeout_lns 600 --timeout_miqp 0
```
- Uses Large Neighborhood Search metaheuristic with constructive initialization
- Fast execution, good for large instances
- Set `--timeout_miqp 0` to skip exact optimization

**3. MIQP-Only Algorithm (Exact Optimization)**
```bash
cd methods/build
./evcs_optimizer --scenario 1 --objective_index 0 --timeout_lns 0 --timeout_miqp 3600
```
- Uses only Gurobi's exact mixed-integer optimization
- Optimal solutions for smaller instances
- Set `--timeout_lns 0` to skip heuristic initialization

**4. LNS+MIQP Algorithm (Hybrid Optimization)**
```bash
cd methods/build
./evcs_optimizer --scenario 1 --objective_index 0 --timeout_lns 600 --timeout_miqp 3000
```
- LNS provides warm start for Gurobi optimization
- Best of both worlds: speed + optimality
- Default approach for most experiments

### Objective Functions

The system supports multiple fairness and efficiency objectives:

| Index | Objective | Description |
|-------|-----------|-------------|
| 0 | **Utilitarian** | Minimize total relative shortfall |
| 1 | **Quadratic** | Minimize sum of squared shortfalls |
| 2 | **Proportional Fairness** | Minimize sum of log shortfalls |
| 3 | **Harmonic Fairness** | Minimize sum of log shortfalls |
| 4 | **Max-Min** | Minimize maximum shortfall |
| 5 | **RMD Hybrid** | Convex combination of RMD + utility |
| 6 | **Gini Hybrid** | Convex combination of Gini + utility |
| 7 | **Jain Hybrid** | Convex combination of Jain's fairness + utility |
| 8 | **Envy Hybrid** | Convex combination of envy-freeness + utility |
| - | **Egalitarian** | Baseline with 100% fairness and 0% utility |

### Data Requirements

The executable expects CSV data files in the working directory:
- `ev_scenario-{scenario}.csv` - Scenario configuration
- `station_{scenario}.csv` - Charging station configuration
- `PV-{scenario}.csv` - Solar production profile


### Advanced Parameters

**LNS Algorithm Tuning:**
```bash
./evcs_optimizer \
  --scenario 1 --objective_index 0 \
  --timeout_lns 600 --timeout_miqp 3000 \
  --start_temp 100.0 --end_temp 0.001 \
  --max_trials 10 --max_generated 10000 \
  --seed 1
```

**Output Control:**
```bash
./evcs_optimizer \
  --scenario 1 --objective_index 0 \
  --final_output_filename "my_solution.log" \
  --lns_output_filename "lns_only.log" \
  --objective_evolution_filepath "evolution.csv"
```

**Multi-Scenario Parameters:**
```bash
./evcs_optimizer \
  --ev_scenario 1 --pv_scenario 1 --station_scenario 1 \
  --objective_index 0 --timeout_lns 300 --timeout_miqp 1800
```

### Complete Parameter Reference

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--scenario` | int | (REQUIRED) | Overall scenario ID (sets all sub-scenarios) |
| `--ev_scenario` | int | [scenario] | EV arrival/departure scenario |
| `--pv_scenario` | int | [scenario] | Solar production scenario |
| `--station_scenario` | int | [scenario] | Charging station configuration |
| `--nb_ev` | int | -1 | Number of EVs (overrides scenario file) |
| `--objective_index` | int | (REQUIRED) | Objective function (see table above) |
| `--timeout_lns` | double | 300 | LNS time limit (seconds) |
| `--timeout_mip` | double | 3300 | MIP time limit (seconds) |
| `--timeout_neighbor` | double | 10 | Neighborhood search time limit |
| `--start_temp` | double | 100.0 | Initial temperature for simulated annealing |
| `--end_temp` | double | 0.001 | Final temperature |
| `--max_trials` | int | 10 | Maximum trials per temperature |
| `--max_generated` | int | 10000 | Maximum solutions generated per trial |
| `--seed` | int | -1 | Random seed (-1 for current time) |
| `--final_output_filename` | string | "complete_solution_[scenario].log" | Final solution output file |
| `--lns_output_filename` | string | "lns_solution_[scenario].log" | LNS-only solution output file |
| `--objective_evolution_filepath` | string | - | Objective evolution tracking file |

**Troubleshooting**

If you encounter Gurobi-related errors:

1. **Verify Gurobi installation:**
   ```bash
   ls $GUROBI_HOME/lib/  # Should show libgurobi*.so files
   ls $GUROBI_HOME/include/  # Should show gurobi_c++.h
   ```

2. **Check CMake output:**
   Look for the message `"Using GUROBI_HOME: /path/to/gurobi"` during cmake configuration.

3. **License issues:**
   Ensure your Gurobi license is properly configured:
   ```bash
   gurobi_cl --version  # Should show license info
   ```

---

## Reproducing experiments

* **Instance sets**: multiple groups (e.g., G1–G6) with varying number of EVs.
* **PV data**: hourly PV profiles resampled to 6-minute steps (τ = 0.1 h).
* **Time limits** (typical): Gurobi 3600s, LNS 600s, hybrid LNS→Gurobi (600s + 3000s).
* **Defaults**: proportional utility `U^(3)`, convex combinations with `λ = 0.5`.

> Exact configs, seeds, and figure scripts will be included to regenerate utility-fairness plots, GAP tables, and correlation matrices.

---

## How to cite

If you use this repo, please cite the paper:

```bibtex

```

---

## License
MIT License

Copyright (c) 2025 Rémi Gauchotte

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.





