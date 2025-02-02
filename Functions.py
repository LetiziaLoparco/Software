import numpy as np
import json
from pulser import Pulse, Sequence, Register
from pulser_simulation import QutipEmulator
from pulser.waveforms import RampWaveform
from pulser.devices import AnalogDevice
from Utilities import *
from SharedVariables import IMAGE_MATRIX_LIST, IMAGE_MATRIX_LIST_INDEX


def load_config(config_path="config.json"):
    """
    Load configuration settings from a JSON file.

    This function reads the specified JSON configuration file and returns its contents as a dictionary.
    If no file path is provided, it defaults to "config.json".

    Parameters
    ----------
    config_path : str, optional
        Path to the configuration file (default is "config.json").

    Returns
    -------
    dict
        A dictionary containing the configuration parameters read from the JSON file.

    Raises
    ------
    FileNotFoundError
        If the specified configuration file does not exist.
    json.JSONDecodeError
        If the file is not a valid JSON format.

    """
    with open(config_path, 'r') as file:
        return json.load(file)



def save_config(config, config_path="config.json"):
    """
    Save configuration settings to a JSON file.

    This function writes the given configuration dictionary to a JSON file, formatting it for readability.
    If no file path is specified, it defaults to "config.json".

    Parameters
    ----------
    config : dict
        A dictionary containing configuration parameters to be saved.
    config_path : str, optional
        Path to the configuration file (default is "config.json").

    Returns
    -------
    None

    Raises
    ------
    TypeError
        If `config` is not a dictionary.
    IOError
        If an error occurs while writing to the file.

    """
    with open(config_path, 'w') as file:
        json.dump(config, file, indent=4)



def prepare_register(U, N_SIDE):
    """
    Prepare a quantum register for a given interaction strength.

    Parameters
    ----------
    U : float
        Interaction strength used to calculate the Rydberg blockade radius.

    N_SIDE : int
        Number of qubits per side of the square register, defining a grid of (N_SIDE x N_SIDE) qubits.

    Returns
    -------
    tuple
        A tuple containing:
        - reg (Register): Quantum register with qubit positions.
        - R_interatomic (float): Interatomic distance scaling factor.
    """
    R_interatomic = AnalogDevice.rydberg_blockade_radius(U)
    reg = Register.square(N_SIDE, R_interatomic, prefix="q")
    return reg, R_interatomic



def run_single_simulation(reg, R_interatomic, Omega_max, delta_0, delta_f, t_rise, t_fall, t_sweep, N_SIDE):
    """
    Run a single quantum simulation for a given sweep time and calculate the final states and correlation values.

    Parameters
    ----------
    reg : Register
        Quantum register containing qubit positions.
    R_interatomic : float
        Interatomic distance scaling factor.
    Omega_max : float
        Maximum Rabi frequency for the pulse.
    delta_0 : float
        Initial detuning of the pulse.
    delta_f : float
        Final detuning of the pulse.
    t_rise : float
        Rise time of the pulse.
    t_fall : float
        Fall time of the pulse.
    t_sweep : float
        Sweep time to simulate.
    N_SIDE : int
        Number of qubits per side of the square register, defining a grid of (N_SIDE x N_SIDE) qubits.

    Returns
    -------
    tuple
        A tuple containing:
        - t_tot (float): Total duration of the sweep in microseconds.
        - sim_results_states (QutipEmulator.Result): Simulation results containing quantum states.
        - correlation_value (dict): Calculated correlation function values.
    """
    t_tot = (t_sweep + t_rise + t_fall) * 1e-3  # Convert to µs

    # Define the pulses
    rise = Pulse.ConstantDetuning(
        RampWaveform(t_rise, 0.0, Omega_max), delta_0, 0.0
    )
    sweep = Pulse.ConstantAmplitude(
        Omega_max, RampWaveform(t_sweep, delta_0, delta_f), 0.0
    )
    fall = Pulse.ConstantDetuning(
        RampWaveform(t_fall, Omega_max, 0.0), delta_f, 0.0
    )

    # Create the sequence and run the simulation
    seq = Sequence(reg, AnalogDevice)
    seq.declare_channel("ising", "rydberg_global")
    seq.add(rise, "ising")
    seq.add(sweep, "ising")
    seq.add(fall, "ising")

    simul = QutipEmulator.from_sequence(seq, sampling_rate=0.02)
    sim_results_states = simul.run(progress_bar=True)

    # Calculate the correlation values and final states
    sim_final_state = sim_results_states.states[-1]
    correlation_value = get_full_corr_function(reg, sim_final_state, R_interatomic, N_SIDE)

    return t_tot, sim_results_states, correlation_value



def run_simulation(reg, R_interatomic, Omega_max, delta_0, delta_f, t_rise, t_fall, t_sweep_range, N_SIDE):
    """
    Run a loop over multiple sweep times and aggregate the results of quantum simulations.

    Parameters
    ----------
    reg : Register
        Quantum register containing qubit positions.
    R_interatomic : float
        Interatomic distance scaling factor.
    Omega_max : float
        Maximum Rabi frequency for the pulse.
    delta_0 : float
        Initial detuning of the pulse.
    delta_f : float
        Final detuning of the pulse.
    t_rise : float
        Rise time of the pulse.
    t_fall : float
        Fall time of the pulse.
    t_sweep_range : list or array-like
        Array of sweep times to simulate.
    N_SIDE : int
        Number of qubits per side of the square register, defining a grid of (N_SIDE x N_SIDE) qubits.
    Returns
    -------
    list of tuple
        A list where each tuple contains:
        - t_tot (float): Total duration of the sweep in microseconds.
        - sim_results_states (QutipEmulator.Result): Simulation results containing quantum states.
        - correlation_value (dict): Calculated correlation function values.
    """
    output_simulations = []

    for t_sweep in t_sweep_range:
        result = run_single_simulation(reg, R_interatomic, Omega_max, delta_0, delta_f, t_rise, t_fall, t_sweep, N_SIDE)
        output_simulations.append(result)

    return output_simulations



def prepare_and_show_first_figure_correlation_matrix(output_simulations, N_SIDE, ax=None): 
    """
    Prepare correlation matrix figure for for all different total time and display the first correlation matrix figure.

    Parameters
    ----------
    output_simulations : list of tuple
        A list of simulation results containing:
        - t_tot (float): Total duration of the sweep in microseconds.
        - sim_results_states (QutipEmulator.Result): Simulation results.
        - correlation_value (dict): Correlation function values.
    N_SIDE : int
        Number of qubits per side of the square register, defining a grid of (N_SIDE x N_SIDE) qubits.
    ax : matplotlib.axes.Axes, optional
        Axes on which to draw the plot. If None, a new plot is created.

    Returns
    -------
    None
    """
    
    for t_tot, _, correlation_value in output_simulations:
        IMAGE_MATRIX_LIST.append((prepare_correlation_matrix(correlation_value, N_SIDE), t_tot)) 
        
    create_figure_correlation_matrix(IMAGE_MATRIX_LIST[IMAGE_MATRIX_LIST_INDEX][0], IMAGE_MATRIX_LIST[IMAGE_MATRIX_LIST_INDEX][1], N_SIDE, ax = ax)



def prepare_and_show_Neel_structure_factor(reg, R_interatomic, output_simulations,  N_SIDE, ax=None):

    """
    Prepare and display the Néel structure factor as a function of time sweep.

    Parameters
    ----------
    reg : Register
        Quantum register containing qubit positions.
    R_interatomic : float
        Interatomic distance scaling factor.
    output_simulations : list of tuple
        A list of simulation results containing:
        - t_tot (float): Total duration of the sweep in microseconds.
        - sim_results_states (QutipEmulator.Result): Simulation results.
        - correlation_value (dict): Correlation function values.
    N_SIDE : int
        Number of qubits per side of the square register, defining a grid of (N_SIDE x N_SIDE) qubits.
    ax : matplotlib.axes.Axes, optional
        Axes on which to draw the plot. If None, a new plot is created.

    Returns
    -------
    None
    """
    neel_structure_factors_list = []
    t_tot_list = []

    
    for t_tot, sim_results_states, _ in output_simulations:
        neel_factor = get_neel_structure_factor(reg, R_interatomic, sim_results_states.states[-1], N_SIDE)
        neel_structure_factors_list.append(neel_factor)
        t_tot_list.append(t_tot)

  
    neel_structure_factors = np.array(neel_structure_factors_list)
    t_tot = np.array(t_tot_list)

    create_figure_Neel_structure_factor(t_tot, neel_structure_factors, ax = ax)

