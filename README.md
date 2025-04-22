# Genetic Algorithm Optimized Energy Utilization System

**Solar Energy Simulation and Load Shedding Optimization using MATLAB**

A MATLAB-based project for simulating solar-powered energy systems, managing battery storage, and optimizing load shedding using Genetic Algorithms (GA).

## 🌞 Project Overview

This project models a **solar-powered energy system** with battery storage and dynamic load management. It simulates solar energy generation using synthetic weather data, prioritizes critical loads during shortages, and applies a **Genetic Algorithm** to optimize load shedding schedules.

## 🔧 Key Features

- **Solar & Battery Modeling**: Simulates solar panel output and battery charge/discharge cycles.
- **Load Shedding Logic**: Prioritizes loads and sheds based on necessity (non‑essential first).
- **Genetic Algorithm**: Minimizes disruption by optimizing the load shedding schedule.
- **Visual Outputs**: Plots and tables showing solar power, battery SOC, and shed loads.

## 📦 Features

- ✅ Synthetic weather data generation (irradiance and temperature)
- ✅ Customizable system inputs (battery capacity, solar panels, load profiles)
- ✅ Interactive or static simulations
- ✅ Hourly visualization of load shedding and system performance

## 🚀 Installation

### 1. Clone the Repository

```bash
git clone https://github.com/your-username/solar-energy-simulation.git
cd genetic-algorithm-optimized-energy-utilization-system
```

### 2. Generate Weather Data

Run in MATLAB:

```matlab
generate_weather_data
```

This creates `weatherData.mat`, required for simulations.

## 💻 MATLAB Requirements

- MATLAB **R2020b** or newer
- No additional toolboxes required

## ⚙️ Usage

1. **Simple Simulation**

   ```matlab
   simple_main_script
   ```

   Executes a simulation using default system parameters.

2. **Dynamic Input Simulation**

   ```matlab
   dynamic_main_script
   ```

   Prompts for:

   - Battery capacity (e.g., 2000 Wh)
   - Number of solar panels (e.g., 5)
   - Solar panel rating (e.g., 300 W)
   - Load profiles (12×24 matrix: rows = circuits, columns = hours; 1 = on, 0 = off)

3. **Optimized Load Shedding (Genetic Algorithm)**

   ```matlab
   optimized_main_script
   ```

   Uses a Genetic Algorithm to compute the most efficient load‑shedding schedule.

## 📁 File Structure

| File                       | Purpose                                           |
| -------------------------- | ------------------------------------------------- |
| `dynamic_main_script.m`    | Interactive simulation with user-defined inputs   |
| `optimized_main_script.m`  | GA-optimized load shedding simulation             |
| `simple_main_script.m`     | Basic simulation using default values             |
| `simulate_load_shedding.m` | Core logic for battery and load management        |
| `generate_weather_data.m`  | Generates synthetic weather data                  |
| `house_specs.m`            | Defines household circuits, loads, and priorities |
| `load_profiles.m`          | Example 12×24 matrix for load usage               |
| `visualize_results.m`      | Plots and summarizes simulation outputs           |

## 📊 Outputs

### Plots

- 📈 Solar irradiance and temperature over 24 hours
- 🔥 Heatmap of load shedding (circuits vs. hour)
- 🔋 Solar power output and battery state of charge (SOC)

### Tables

- Hourly solar power generation
- Battery SOC over time
- Total and per-hour load shedding summaries

## 📄 Dependencies

- MATLAB (tested with **R2020b+**)
- `weatherData.mat` (generated via `generate_weather_data.m`)

## 🤝 Contributing

1. **Fork** the repository
2. **Create** a feature branch:
   ```bash
   git checkout -b feature/your-feature-name
   ```
3. **Commit** your changes:
   ```bash
   git commit -m "Add new feature"
   ```
4. **Push** the branch:
   ```bash
   git push origin feature/your-feature-name
   ```
5. **Open** a Pull Request and describe your changes

## 🪪 License

Distributed under the **MIT License**. See [`LICENSE`](LICENSE) for details.

## 📬 Contact

For questions or contributions, please contact: [[your-email@example.com](mailto\:your-email@example.com)]

