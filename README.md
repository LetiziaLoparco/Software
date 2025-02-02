# Quantum Simulation GUI for Register and Correlation Analysis

This software offers an intuitive graphical user interface (GUI) for simulating quantum registers and analyzing their correlation and Néel structure factors. Built with Qt Designer, the application enables users to configure simulation parameters, visualize quantum register arrangements, plot correlation matrices, and compute the Néel structure factor efficiently.

---

## Features

1. **Quantum Register Visualization**:
   - Set up and visualize a quantum register for a given interaction strength.

2. **Correlation Matrix Plotting**:
   - Simulate quantum states and compute correlation matrices.
   - Navigate through multiple time-sweep plots using the GUI.

3. **Néel Structure Factor Calculation**:
   - Compute and visualize the Néel structure factor for time sweeps.

---

## How to Download and Run the Software

### For Linux Users

1. Clone the repository:
   ```bash
   git clone https://github.com/LetiziaLoparco/Simulation.git
   ```

2. Install the required Python packages (Python 3.10 is required):
   ```bash
   pip install numpy matplotlib PyQt5 pulser qutip
   ```

3. Run the GUI:
   ```bash
   python Main.py
   ```

---

### For Anaconda Users

1. Clone the repository:
   ```bash
   git clone https://github.com/LetiziaLoparco/Simulation.git
   ```

2. Install the required packages in the Anaconda environment:
   ```bash
   conda install numpy matplotlib -c conda-forge
   conda install -c conda-forge pyqt pulser qutip
   ```

3. Open the `Main.py` file in Spyder or run it in the terminal:
   ```bash
   python Main.py
   ```

---

## Preparing an Antiferromagnetic State and Observing Correlations

This code and theory summary are based on the Quantum Computing course by Prof. Calarco, Prof. Ercolessi, and Prof. Bonacorsi at the University of Bologna. The code is inspired by a model experimentally implemented with a Rydberg-atom platform ([DOI: 10.1103/PhysRevX.8.021070]). It leverages Pulser, an open-source Python library for programming neutral atom devices at the pulse level ([arXiv:2104.15044v3]).

We aim to study a quantum many-body system, dynamically tuning the parameters to observe the buildup of antiferromagnetic order. The simulation emulates an Ising quantum antiferromagnet on a square lattice (pre-set to $$N_{side}=3$$ in `Config.json`), exploring the ground-state phase diagrams for the nearest-neighbor Ising model.

The GUI enables the observation of the influence of finite ramp speed on the extent of correlations, depending on the chosen parameters. The development of correlations during the ramp is visualized in both space and time.

The system is prepared in a product state with all atoms in the ground state. The Hamiltonian governing its evolution is:

$$ H = \sum_i \left[ \frac{\hbar \Omega(t)}{2} \sigma_i^x - \hbar \delta(t)n_i \right] + \frac{1}{2} \sum_{i \neq j} U_{ij} n_i n_j $$

The parameters are dynamically tuned as follows:

![Pulse_sequence](https://github.com/user-attachments/assets/4ce7d305-e77c-4545-aa4c-9bc6e5a4077c)


### Steps in Detail

#### 1. Parameters Initialization

Simulation parameters can be configured through the GUI. Default values are pre-set. Key parameters include:
- Interaction strength $$\( U \)$$
- Maximum Rabi frequency $$\( \Omega \)$$
- Initial and final detuning $$\( \Delta \)$$
- Rise and fall times
- Time sweep ranges (number of steps is pre-set to 100 in `Config.json`).

##### Note:
All limitations on pulse parameters can be found by printing the specifications from the Pulser library. Use the following command:

```python
from pulser.devices import AnalogDevice
AnalogDevice.print_specs()
```


#### 2. Quantum Register Preparation

The register is set up as a square lattice of $$\( N_{\text{side}} \times N_{\text{side}} \)$$ atoms, where $$\( N_{\text{side}} \)$$ is defined in `Config.json` ( default: $$\( N_{\text{side}} = 3 \)$$ ).

The interatomic distance (lattice spacing) is set equal to the Rydberg blockade radius, computed using Pulser’s `AnalogDevice.rydberg_blockade_radius(U)` function. The register is visualized using `reg.draw()`.

#### 3. Spin-Spin Correlation Function Matrix

The correlation function measures the quality of the resulting state after the simulation:

$$ g^{(2)}(k, l) = \frac{1}{N_{k, l}} \sum_{i, j} \left( \langle n_i n_j \rangle - \langle n_i \rangle \langle n_j \rangle \right) $$

This function evaluates the correlation between local operators $$\( n_i \)$$ and $$\( n_j \)$$ at different positions, averaged over all pairs of atoms separated by vector $$\( (ka, la) \)$$.

Pressing the **"Show 2D antiferromagnetic correlation"** button initiates a new simulation with a laser applied globally following the pulse sequence for different $$\( t_{\text{sweep}} \)$$ values. The correlation matrix is plotted for different time values and displayed using the **"Show Next"** button. Stronger correlations are indicated by more vivid colors.

#### 4. Néel Structure Factor

To understand the behavior further, the Néel structure factor is computed to detect antiferromagnetic correlations in the AF region of the phase diagram:

$$ S_{\text{Néel}} = 4 \sum_{(k, l) \neq (0, 0)} (-1)^{|k|+|l|} g^{(2)}(k, l) $$

In the AF region of the phase diagram, stronger correlations yield higher values of $$\( S_{\text{Néel}} \)$$.

---

## **How to Navigate the Software**

The software provides a **user-friendly interface** to visualize and analyze quantum systems. Follow these instructions to effectively utilize its features.


### **Before Starting**
The **GUI interface** allows you to select a `.json` file via the **"Select Configuration File"** option. This file defines:
- **`N_SIDE`**: The number of atoms per side in a square lattice.
- **`NUMBER_STEPS`**: The step size for the `t_sweep` range.

**If no file is provided, default values will be used:**  
- `N_SIDE = 3`  
- `NUMBER_STEPS = 100`


### **Main Features**
- **Visualize the Quantum Register**  
  Click **"Show Register"** to view the quantum register configuration.
  
- **Correlation Plot**  
  Generate and display the **correlation matrix** for the first time sweep by clicking **"Show 2D Antiferromagnetic Correlation"**.
  
- **Néel Structure Plot**  
  Compute and display the **Néel structure factor** over time by selecting **"Show Néel Structure Factor Depending on t"**.



### **Additional Options**
- **Navigate Time-Sweep Plots**  
  Use **"Show Next"** to cycle through different time-sweep plots.

- **Reset and Re-run Simulations**  
  Modify parameters or restart simulations using the **reset option**.


### Example GUI Interface
Below is an example screenshot showing the GUI after generating the antiferromagnetic correlation plot by clicking the corresponding button.

![GUInterface](https://github.com/user-attachments/assets/06792010-1db1-4d57-a4ea-8bc3aa1fb45d)


---

## **Technical Details**

### **Dependencies**
The software requires the following Python libraries:
- `numpy`
- `matplotlib`
- `PyQt5`
- `pulser`
- `qutip`

### **Core Modules**
- **`Main.py`**: The main entry point of the software.
- **`Handler.py`**: Manages the connections between the GUI layout and core functionalities.
- **`Functions.py`**: Contains the core simulation functions.
- **`Utilities.py`**: Provides helper functions for plotting and computation.
- **`GUI_setup.py`**: Defines and manages the **PyQt5-based** GUI layout.
- **`GUI_setup.ui`**: The GUI layout file in `.ui` format.
- **`Config.json`**: Stores configuration settings, including the number of atoms per side (`N_SIDE`) and the step size for the time sweep (`NUMBER_STEPS`).
- **`SharedVariables.py`**: Contains global variables and acts as a bridge between different Python modules.
- **`Test_Utilities.py`**: Includes utility functions for testing.
- **`Test_Functions.py`**: Contains specific test cases for validating core functionalities.


---

## Future Extensions

- Add support for 3D quantum register visualizations.
- Include additional interaction models.

---

## Acknowledgments

- Inspired by advanced quantum simulation techniques.
- Special thanks to contributors and the open-source community for invaluable resources.

