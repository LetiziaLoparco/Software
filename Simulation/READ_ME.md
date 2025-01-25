# Quantum Simulation GUI for Register and Correlation Analysis

This software provides a graphical user interface (GUI) to simulate quantum registers and analyze their correlation and Néel structure factors. Users can configure simulation parameters, visualize quantum register setups, plot correlation matrices, and compute the Néel structure factor using the GUI.

---

## Features

1. **Quantum Register Visualization**:
   - Set up and visualize a quantum register for a given interaction strength.
   - Customize parameters for the register size and spacing.

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

2. Install the required Python packages:
   ```bash
   pip install numpy matplotlib PyQt5 pulser qutip
   ```

4. Run the GUI:
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

## Parameters configuration

You can configure the simulation parameters through the GUI interface. Some correct values are pre-imposted by default. The parameters are:
- Interaction strength (U)
- Maximum Rabi frequency (Ω)
- Initial and final detuning (Δ)
- Rise and fall times
- Time sweep ranges

---

## How to Navigate in the Software

- **Main Features**:
  - Visualize the quantum register by clicking the **"Register Sequence"** button.
  - Plot the correlation matrix for the first time sweep using the **"Correlation Plot"** button.
  - Compute and display the Néel structure factor with the **"Néel Structure Plot"** button.

- **Additional Options**:
  - Use the **"Next"** button to navigate through time-sweep plots.
  - Reset parameters or re-run simulations as needed.

---
## Visualizing Correlation Functions and Néel Structure Factor

---

## Technical Details

- **Dependencies**:
  - Python libraries: `numpy`, `matplotlib`, `PyQt5`, `pulser`, `qutip`

- **Core Modules**:
  - `Functions.py`: Contains core simulation functions.
  - `Utilities.py`: Includes helper functions for plotting and computation.
  - `GUI_setup.py`: Handles the GUI layout and design.

---

## Future Extensions

- Add support for 3D quantum register visualizations.
- Include additional interaction models.

---

## License

This project is licensed under the MIT License. See the LICENSE file for details.

---

## Acknowledgments

- Inspired by advanced quantum simulation techniques.
- Special thanks to contributors and the open-source community for providing invaluable resources.

