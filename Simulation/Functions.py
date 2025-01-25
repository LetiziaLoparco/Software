import numpy as np
from pulser import Pulse, Sequence, Register
from pulser_simulation import QutipEmulator
from pulser.waveforms import RampWaveform
from pulser.devices import AnalogDevice
from Utilities import *
from Config import *


def prepare_register(U):
    R_interatomic = AnalogDevice.rydberg_blockade_radius(U)
    reg = Register.square(N_SIDE, R_interatomic, prefix="q")
    return reg, R_interatomic


def run_simulation(reg, R_interatomic, Omega_max, delta_0, delta_f, t_rise, t_fall, t_sweep_range):
   
    output_simulations = []

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

        sim_results_states = simul.run(progress_bar=True)
        sim_final_state = sim_results_states.states[-1]
        correlation_value = get_full_corr_function(reg, sim_final_state, R_interatomic) 

        output_simulations.append((t_tot, sim_results_states, correlation_value)) 

    return output_simulations


def prepare_and_show_first_figure_correlation_matrix(output_simulations, ax=None): 
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
    
    for t_tot, _, correlation_value in output_simulations:
        IMAGE_LIST.append((prepare_correlation_matrix(correlation_value), t_tot)) # lista di tuple
        
    create_figure_correlation_matrix(IMAGE_LIST[IMAGE_LIST_INDEX][0], IMAGE_LIST[IMAGE_LIST_INDEX][1], ax = ax)


def prepare_and_show_Neel_structure_factor(reg, R_interatomic, output_simulations, ax=None):

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
    neel_structure_factors_list = []
    t_tot_list = []

    # Commento su cosa fa
    for t_tot, sim_results_states, _ in output_simulations:
        neel_factor = get_neel_structure_factor(reg, R_interatomic, sim_results_states.states[-1])
        neel_structure_factors_list.append(neel_factor)
        t_tot_list.append(t_tot)

  
    neel_structure_factors = np.array(neel_structure_factors)
    t_tot = np.array(t_tot_list)

    Create_figure_Neel_structure_factor(t_tot, neel_structure_factors, ax = ax)

