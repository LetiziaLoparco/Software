import pytest
import numpy as np
from pulser import Register
from pulser.devices import AnalogDevice


from Functions import prepare_register 
from Functions import run_single_simulation  



def test_prepare_register():
    """
    Test the `prepare_register` function with a valid interaction strength.

    The test checks:
    - Correct calculation of R_interatomic.
    - Proper creation of the quantum register.
    """
    N_SIDE = 3 
    U = 10.0 #MHz
    C6 = 5*10**5 #Mhz*um**6 for Rb87 at n=70

    expected_R_interatomic = (C6 / U)**(1/6)
    expected_num_qubits = N_SIDE ** 2  

    # Call the function
    reg, R_interatomic = prepare_register(U, N_SIDE)

    # Assertions
    assert R_interatomic == pytest.approx(expected_R_interatomic, abs =1.5), "R_interatomic is not calculated correctly."
    assert len(reg.qubits) == expected_num_qubits, f"Expected {expected_num_qubits} qubits, but got {len(reg.qubits)}."

####################################################################################################


def test_run_single_simulation_output_properties():
    """
    Test the `run_single_simulation` function with fixed inputs .

    Ensures reproducible results and checks output properties.
    """

    # Define minimal test inputs
    R_interatomic = 6.0
    N_SIDE = 3
    reg = Register.square(N_SIDE, R_interatomic, prefix="q")  
    Omega_max = 1.8 *2 * np.pi 
    delta_0 = -6.0
    delta_f = 4.5
    t_rise = 250.0  
    t_fall = 250.0  
    t_sweep = 500.0  

    # Run the simulation
    t_tot, sim_results_states, correlation_value = run_single_simulation(
        reg, R_interatomic, Omega_max, delta_0, delta_f, t_rise, t_fall, t_sweep, N_SIDE
    )

    # Check total time
    expected_t_tot = (t_sweep + t_rise + t_fall) * 1e-3  # Convert to µs
    assert t_tot == pytest.approx(expected_t_tot, rel=1e-3), "Total time is incorrect."

    # Check that simulation results contain states
    assert hasattr(sim_results_states, "states"), "Simulation results must contain states."
    assert len(sim_results_states.states) > 0, "Simulation results should have at least one state."

    # Check correlation values
    assert isinstance(correlation_value, dict), "Correlation values should be a dictionary."
    assert len(correlation_value) > 0, "Correlation values should not be empty."



def test_run_single_simulation_staggered_profile():
    """
    Test if the final state has at least 60% probability of being in '101010101'
    when measured 1000 times.

    Since the sweep time (t_sweep) is long (1000.0 ns), the system has enough time to evolve 
    towards the target state. The expected final state, '101010101', represents a staggered 
    profile, which indicates strong antiferromagnetic correlations.

    A high probability of measuring this state confirms that the system successfully 
    forms an antiferromagnetic order, as desired.
    
    This test ensures that the quantum register exhibits the correct staggered ordering 
    after a sufficiently long evolution.
    """
    # Define simulation parameters
    U = 2.7*2*np.pi  
    N_SIDE = 3
    R_interatomic = AnalogDevice.rydberg_blockade_radius(U) 
    reg = Register.square(N_SIDE, R_interatomic, prefix="q") 
    Omega_max = 1.8 *2 * np.pi 
    delta_0 = -6.0 *2 * np.pi 
    delta_f = 4.5 *2 * np.pi 
    t_rise = 250.0  
    t_fall = 250.0  
    t_sweep = 1000.0  

    # Run the simulation
    _, sim_results_states, _ = run_single_simulation(
        reg, R_interatomic, Omega_max, delta_0, delta_f, t_rise, t_fall, t_sweep, N_SIDE
    )

    # Sample final state 1000 times
    num_samples = 1000
    counts = sim_results_states.sample_final_state(N_samples=num_samples)

    # Define the target state
    target_state = "101010101"

    # Count occurrences of the target state
    target_count = counts.get(target_state, 0)
    percentage = (target_count / num_samples) * 100  # Convert to percentage

    # Check if at least 60% of the samples are in the target state
    assert percentage >= 60, f"Test failed! Probability of '{target_state}' is only {percentage:.2f}%"



def test_run_single_simulation_ground_state():
    """
    Test if the final state has at least 60% probability of NOT being in '101010101'
    when measured 1000 times.

    Since the sweep time (t_sweep) is short (100.0 ns), the system does not have enough 
    time to fully evolve into the Rydberg state. Instead, the final state should remain 
    close to the ground state, meaning that the antiferromagnetic ordering ('101010101') 
    should appear infrequently.

    This test ensures that, with a small evolution time, the atoms remain mostly 
    in their initial ground state, validating the expected slow dynamical evolution 
    of the system.
    """
    # Define simulation parameters
    U = 2.7 * 2 * np.pi  
    N_SIDE = 3  
    R_interatomic = AnalogDevice.rydberg_blockade_radius(U)  
    reg = Register.square(N_SIDE, R_interatomic, prefix="q")  
    Omega_max = 1.8 * 2 * np.pi  
    delta_0 = -6.0 * 2 * np.pi  
    delta_f = 4.5 * 2 * np.pi  
    t_rise = 250.0  
    t_fall = 250.0  
    t_sweep = 100.0  

    # Run the simulation
    _, sim_results_states, _ = run_single_simulation(
        reg, R_interatomic, Omega_max, delta_0, delta_f, t_rise, t_fall, t_sweep, N_SIDE
    )

    # Sample final state 1000 times
    num_samples = 1000
    counts = sim_results_states.sample_final_state(N_samples=num_samples)

    # Define state (that should NOT appear frequently)
    rydberg_target_state = "101010101"  

    # Count occurrences
    rydberg_count = counts.get(rydberg_target_state, 0)
    rydberg_percentage = (rydberg_count / num_samples) * 100  # Probability of '101010101'
    non_rydberg_percentage = 100 - rydberg_percentage  # Probability of NOT being '101010101'

    # Ensure that at least 60% of the states are NOT in "101010101"
    assert non_rydberg_percentage >= 60, (
        f"Test failed! Only {non_rydberg_percentage:.2f}% of samples are NOT in '{rydberg_target_state}'"
    )


