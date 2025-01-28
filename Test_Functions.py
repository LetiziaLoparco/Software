import pytest
import numpy as np
import qutip
from pulser import Register
from pulser.devices import AnalogDevice


from Functions import prepare_register 
from Functions import run_single_simulation  
from Functions import get_corr_function


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
    reg, R_interatomic = prepare_register(U)

    # Assertions
    assert R_interatomic == pytest.approx(expected_R_interatomic, abs =1.5), "R_interatomic is not calculated correctly."
    assert len(reg.qubits) == expected_num_qubits, f"Expected {expected_num_qubits} qubits, but got {len(reg.qubits)}."

####################################################################################################


def test_run_single_simulation():
    """
    Test the `run_single_simulation` function with fixed inputs and a random seed.

    Ensures reproducible results and checks output properties.
    """
    # Set random seed for reproducibility
    np.random.seed(42)

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
        reg, R_interatomic, Omega_max, delta_0, delta_f, t_rise, t_fall, t_sweep
    )

    # Assertions: Check outputs
    expected_t_tot = (t_sweep + t_rise + t_fall) * 1e-3  # Convert to µs
    assert t_tot == pytest.approx(expected_t_tot, rel=1e-3), "Total time is incorrect."

    # Check that simulation results contain states
    assert hasattr(sim_results_states, "states"), "Simulation results must contain states."
    assert len(sim_results_states.states) > 0, "Simulation results should have at least one state."

    # Check correlation values
    assert isinstance(correlation_value, dict), "Correlation values should be a dictionary."
    assert len(correlation_value) > 0, "Correlation values should not be empty."


"""
def test_run_single_simulation_ground_state():
    
    #Test that for a low sweep time (e.g., t_sweep = 50 µs), 
    #most qubits remain in the ground state (approximately 70% or more).
    

    np.random.seed(42)

    N_SIDE = 3
    U = 2.7 * 2 * np.pi  # MHz
    reg, R_interatomic = prepare_register(U)   
    Omega_max = 2 * np.pi * 1.8  # MHz
    delta_0 = -6.0  
    delta_f = 4.5  
    t_rise = 250.0
    t_fall = 250.0  
    t_sweep = 520.0  

    # Run the simulation
    t_tot, sim_results_states, correlation_value = run_single_simulation(
        reg, R_interatomic, Omega_max, delta_0, delta_f, t_rise, t_fall, t_sweep
    )

    # Final state from the simulation
    final_state = sim_results_states.states[-1]

    # Multi-qubit operators for each qubit
    num_qubits = len(reg.qubits)
    operators = [
        qutip.tensor(
            [qutip.qeye(2) if j != i else qutip.basis(2, 0) * qutip.basis(2, 0).dag()
             for j in range(num_qubits)]
        )
        for i in range(num_qubits)
    ]

    # Calculate the expectation values for the ground state
    ground_state_expectations = [
        qutip.expect(operators[i], final_state) for i in range(num_qubits)
    ]

    # Count qubits in the ground state (expectation > 0.7)
    qubits_in_ground_state = sum(1 for exp in ground_state_expectations if exp > 0.7)

    # Assert that at least 70% of qubits are in the ground state
    assert qubits_in_ground_state >= 0.7 * num_qubits, (
        f"Expected at least 70% of qubits in the ground state, "
        f"but got {qubits_in_ground_state / num_qubits * 100:.1f}%"
    )

"""