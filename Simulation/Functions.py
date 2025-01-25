import numpy as np
import matplotlib.pyplot as plt
from pulser import Pulse, Sequence, Register
from pulser_simulation import QutipEmulator
from pulser.waveforms import RampWaveform
from pulser.devices import AnalogDevice
from Utilities import *
#from Main import main_window
from Config import *


def prepare_register(U):
    R_interatomic = AnalogDevice.rydberg_blockade_radius(U)
    reg = Register.square(N_SIDE, R_interatomic, prefix="q")
    return reg, R_interatomic

def run_simulation(reg, R_interatomic, Omega_max, delta_0, delta_f, t_rise, t_fall, t_sweep_range):
   
    results_correlation = []

    for t_sweep in t_sweep_range:
        t_tot = (t_sweep + t_rise + t_fall)*1e-3 # Convert to µs

        rise = Pulse.ConstantDetuning(
            RampWaveform(t_rise, 0.0, Omega_max), delta_0, 0.0
        )
        sweep = Pulse.ConstantAmplitude(
            Omega_max, RampWaveform(t_sweep, delta_0, delta_f), 0.0
        )
        fall = Pulse.ConstantDetuning(
            RampWaveform(t_fall, Omega_max, 0.0), delta_f, 0.0
        )

        seq = Sequence(reg, AnalogDevice)
        seq.declare_channel("ising", "rydberg_global")
        seq.add(rise, "ising")
        seq.add(sweep, "ising")
        seq.add(fall, "ising")
        simul = QutipEmulator.from_sequence(seq, sampling_rate=0.02)
        results = simul.run(progress_bar=True)
        final_state = results.states[-1]
        correlation_function = get_full_corr_function(reg, final_state, R_interatomic)

        results_correlation.append((t_tot, results, correlation_function))

    return results_correlation

def plot_correlation_vs_t(results_correlation, ax=None):
    """
    Plot the correlation function for varying time sweeps.

    Parameters
    ----------
    correlation_function : dict
        Dictionary of correlation function values.
    t_sweep : float
        Duration of the time sweep.
    t_rise : float
        Rise time of the pulse.
    t_fall : float
        Fall time of the pulse.
    N_side : int
        Number of sites per side in the lattice.

    Returns
    -------
    None
    """
    
    for t_tot, _, correlation_function in results_correlation:
        IMAGE_LIST.append((prepare_correlation_matrix(correlation_function, N_SIDE), t_tot))
        
    Create_figure_correlation_matrix(IMAGE_LIST[IMAGE_LIST_INDEX][0], IMAGE_LIST[IMAGE_LIST_INDEX][1], N_SIDE, ax = ax)
    """global ax_2
    ax_2 = ax"""


"""def show_next_correlation_plot(ax=None):
    global image_list_index 
    global ax_2
    
    if image_list_index < len(image_list):
        image_list_index += 1
        print(image_list_index)
        print(N_side)
        print(ax_2)
        Create_figure_correlation_matrix(image_list[image_list_index][0], image_list[image_list_index][1], N_side, ax = ax_2)"""


def plot_Neel_structure_factor_different_t(reg, R_interatomic, correlation_function, ax=None):

    """
    Plot the Néel structure factor as a function of time sweep.

    Parameters
    ----------
    results_final_state : list of tuple
        Each tuple contains a time sweep value and the corresponding simulation results.
    t_rise : float
        Rise time of the pulse.
    t_fall : float
        Fall time of the pulse.
    reg : object
        Quantum register containing qubit positions.
    R_interatomic : float
        Interatomic distance scaling factor.

    Returns
    -------
    None
    """
    neel_structure_factors = []
    t_tot_list = []

    for t_tot, results, _ in correlation_function:
        neel_factor = get_neel_structure_factor(reg, R_interatomic, results.states[-1])
        neel_structure_factors.append(neel_factor)
        t_tot_list.append(t_tot)

  
    neel_structure_factors = np.array(neel_structure_factors)
    t_tot = np.array(t_tot_list)

    Create_figure_Neel_structure_factor(t_tot, neel_structure_factors, ax = ax)

