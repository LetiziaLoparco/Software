import numpy as np
import qutip
import pytest

from pulser import Register


from Utilities import occupation  
from Utilities import get_corr_pairs
from Utilities import get_corr_function
from Utilities import prepare_correlation_matrix
from Utilities import get_neel_structure_factor


####################################################################################################

def test_occupation_valid_index():
    """
    Test the `occupation` function with valid indices within the range of the quantum register.

    The function is tested with indices that are valid for a given number of qubits. We expect
    the resulting operator to have the correct tensor structure with the occupation operator
    applied only to the specified site.
    """
    # Starting data
    j = 2  # Index of the site
    N_qubits = 5  # Total number of qubits in the register

    # Expected result
    up = qutip.basis(2, 0)
    expected_operator = qutip.tensor( [qutip.qeye(2), qutip.qeye(2), up * up.dag(), qutip.qeye(2), qutip.qeye(2)] )

    # Calculated result
    result = occupation(j, N_qubits)

    # Assert equality
    assert result == expected_operator, "Occupation operator does not match the expected structure."


def test_occupation_first_site():
    """
    Test the `occupation` function for the first site in the quantum register.

    The function is tested with `j=0` to verify that the occupation operator is applied
    to the first qubit, leaving the rest as identity operators.
    """
    # Starting data
    j = 0  # First site
    N_qubits = 3  # Total number of qubits

    # Expected result
    up = qutip.basis(2, 0)
    expected_operator = qutip.tensor([up * up.dag(), qutip.qeye(2), qutip.qeye(2)])

    # Calculated result
    result = occupation(j, N_qubits)

    # Assert equality
    assert result == expected_operator, "Occupation operator is incorrect for the first site."


def test_occupation_last_site():
    """
    Test the `occupation` function for the last site in the quantum register.

    The function is tested with `j=N_qubits-1` to verify that the occupation operator
    is applied to the last qubit, leaving the others as identity operators.
    """
    # Starting data
    j = 4  # Last site
    N_qubits = 5  # Total number of qubits

    # Expected result
    up = qutip.basis(2, 0)
    expected_operator = qutip.tensor(
        [qutip.qeye(2), qutip.qeye(2), qutip.qeye(2), qutip.qeye(2), up * up.dag()]
    )

    # Calculated result
    result = occupation(j, N_qubits)

    # Assert equality
    assert result == expected_operator, "Occupation operator is incorrect for the last site."


def test_occupation_invalid_index():
    """
    Test the `occupation` function with an invalid index (out of range).

    If the index `j` is greater than the number of qubits,
    the function is expected to raise an `IndexError`.
    """
    # Invalid index
    j = 6  # Out of range
    N_qubits = 5  # Total number of qubits

    # Assert that an IndexError is raised
    with pytest.raises(IndexError, match="Index out of range"):
        occupation(j, N_qubits)

def test_occupation_raises_negative_index():
    """
    Test that the `occupation` function raises an IndexError for a negative index.
    """
    j = -1  # Negative index
    N_qubits = 5

    # Assert that an IndexError is raised
    with pytest.raises(IndexError, match="Index out of range"):
        occupation(j, N_qubits)

def test_occupation_single_qubit():
    """
    Test the `occupation` function for a single-qubit register.

    The function is tested with `N_qubits=1` to verify that the occupation operator
    is created correctly for the single qubit.
    """
    # Starting data
    j = 0  # Only site in the register
    N_qubits = 1  # Single-qubit register

    # Expected result
    up = qutip.basis(2, 0)
    expected_operator = up * up.dag()

    # Calculated result
    result = occupation(j, N_qubits)

    # Assert equality
    assert result == expected_operator, "Occupation operator is incorrect for a single-qubit register."


####################################################################################################

def test_get_corr_pairs_valid():
    """
    Test the `get_corr_pairs` function with valid inputs and expected pairs.

    The function is tested with a simple register layout and displacement values (k, l).
    Expected pairs are calculated based on the distance condition.
    """
    # Starting data
    R_interatomic = 1.0  
    reg = Register.rectangle(2, 2, spacing=R_interatomic)  # Create a 2x2 lattice
    k, l = 1, 0  # Displacement vector

    # Expected pairs: [(1,0), (3,2)]
    expected_pairs = [[1,0], [3,2]]

    # Calculated result
    result_pairs = get_corr_pairs(k, l, reg, R_interatomic)

    # Assert that the result matches the expected pairs
    assert result_pairs == expected_pairs, f"Expected {expected_pairs}, got {result_pairs}"

def test_get_corr_pairs_zero_displacement_vector(): 
    """
    Test the `get_corr_pairs` function with valid inputs and expected pairs.

    The function is tested with a simple register layout and displacement values (k, l).
    Expected pairs are calculated based on the distance condition.
    """
    # Starting data
    R_interatomic = 1.0 
    reg = Register.rectangle(2, 2, spacing=R_interatomic)  
    k, l = 0, 0  

    # Expected pairs
    expected_pairs = [[0,0], [1,1], [2,2], [3,3]]

    # Calculated result
    result_pairs = get_corr_pairs(k, l, reg, R_interatomic)

    # Assert that the result matches the expected pairs
    assert result_pairs == expected_pairs, f"Expected {expected_pairs}, got {result_pairs}"


def test_get_corr_pairs_diagonal_pairs(): 
    """
    Test the `get_corr_pairs` function with diagonal displacement vector.

    The function is tested with a register and displacement values where no qubits
    should meet the distance requirement.
    """

    # Starting data
    R_interatomic = 2.0 
    N_SIDE = 3
    reg = Register.square(N_SIDE, spacing=R_interatomic)  
    k, l = 1, 1  

    # Expected result
    expected_pairs = [[4, 0], [5, 1], [7, 3], [8, 4]]

    # Calculated result
    result_pairs = get_corr_pairs(k, l, reg, R_interatomic)

    # Assert that the result 
    assert result_pairs == expected_pairs, f"Expected no pairs, got {result_pairs}"


def test_get_corr_pairs_single_qubit():
    """
    Test the `get_corr_pairs` function with a single-qubit register.

    The function is tested with a single qubit in the register, where no pairs
    can be formed regardless of the displacement.
    """
    # Starting data
    R_interatomic = 1.0
    reg = Register.from_coordinates([[0, 0]]) 
    k, l = 1, 0  

    # Expected result: no pairs
    expected_pairs = []

    # Calculated result
    result_pairs = get_corr_pairs(k, l, reg, R_interatomic)

    # Assert that the result is an empty list
    assert result_pairs == expected_pairs, f"Expected no pairs, got {result_pairs}"


def test_get_corr_pairs_multiple_displacements():
    """
    Test the `get_corr_pairs` function with various displacement values.

    The function is tested with a 3x3 register and multiple displacement vectors
    to ensure correct pair generation for each case.
    """
    # Starting data
    R_interatomic = 1.0
    N_SIDE = 3
    reg = Register.square(N_SIDE, spacing=R_interatomic)  
    displacements = [(1, 0), (0, 1), (1, 1)]  

    # Expected results for each displacement
    expected_results = {
        (1, 0): [[1, 0], [2, 1], [4, 3], [5, 4], [7, 6], [8, 7]],
        (0, 1): [[3, 0], [4, 1], [5, 2], [6, 3], [7, 4], [8, 5]],
        (1, 1): [[4, 0], [5, 1], [7, 3], [8, 4]],
    }

    for (k, l), expected_pairs in expected_results.items():
        # Calculated result for each displacement
        result_pairs = get_corr_pairs(k, l, reg, R_interatomic)
        assert result_pairs == expected_pairs, f"Displacement ({k}, {l}): Expected {expected_pairs}, got {result_pairs}"




####################################################################################################

def test_get_corr_function_entangled_states(): 
    """
    Test the `get_corr_function` for a quantum register initialized with an entangled state.

    This test validates the calculation of the correlation function, a key measure of 
    quantum correlations between qubits, for a register with an entangled state. The test 
    uses a predefined 9-qubit system where the entangled state is constructed as a 
    superposition of two alternating spin configurations: |010101010> and |101010101>.
    """
    # Constants
    R_interatomic = 1.0
    N_SIDE = 3 
    k, l = 1, 0  
    N = N_SIDE * N_SIDE
    # Create a square register of qubits
    reg = Register.square(N_SIDE, spacing=R_interatomic)


    # Define the basis states for 4 qubits
    basis_1 = qutip.tensor(*[qutip.basis(2, i % 2) for i in range(N)])  # |010101010>
    basis_2 = qutip.tensor(*[qutip.basis(2, (i + 1) % 2) for i in range(N)])  # |101010101>

    # Create the entangled state
    final_state = (1 / np.sqrt(2)) * (basis_1 + basis_2)


 
    expected_covariance = -0.25


    # Calculate the correlation function using the tested function
    result = get_corr_function(k, l, reg, R_interatomic, final_state, N_SIDE)

    # Print debugging information
    print(f"Result: {result}, Expected: {expected_covariance}")

    # Assert the result matches the expected value
    assert np.isclose(result, expected_covariance), (
        f"Expected {expected_covariance}, got {result}"
    )


def test_get_corr_function_uncorrelated_qubits():
    """
    Test the `get_corr_function` for a quantum register with an uncorrelated state.

    This test verifies that the correlation function correctly evaluates to zero 
    for a system where qubits are in an uncorrelated product state. The state 
    used is a simple alternating spin configuration: |010101010>, where each qubit 
    is either in the Rydberg state (|0> = |r>) or the ground state (|1>) based on its index.

    The test focuses on the displacement vector (k, l) = (1, 0), representing 
    nearest neighbors along the first axis. Since the qubits are in a product state, 
    no correlations exist, and the expected covariance is zero.

    This test ensures that the function accurately handles scenarios without 
    quantum correlations, validating its correctness for non-entangled states.
    """
    # Starting data
    R_interatomic = 1.0  
    N_SIDE = 3
    reg = Register.square(N_SIDE, spacing=R_interatomic)  
    k, l = 1, 0  

    # Create a simple final state: all qubits in the |r> state
    N = N_SIDE * N_SIDE
    final_state = qutip.tensor(*[qutip.basis(2, i % 2) for i in range(N)])  # |010101010>

    # Expected correlation value
    expected_covariance = 0

    # Calculated result
    result = get_corr_function(k, l, reg, R_interatomic, final_state, N_SIDE)
    
    
    # Assert equality
    assert np.isclose(result, expected_covariance), f"Expected {expected_covariance}, got {result}"

def test_get_corr_function_uncorrelated_qubits_with_interaction():
    """
    Test the `get_corr_function` for a quantum register with an uncorrelated state.

    This test verifies that the correlation function correctly evaluates to zero 
    for a system where qubits are in an uncorrelated product state. The state 
    used is a simple alternating spin configuration: |010101010>, where each qubit 
    is either in the Rydberg state (|0> = |r>) or the ground state (|1>) based on its index.

    The test focuses on the displacement vector (k, l) = (1, 1). Even though for 
    pairs such as (4, 0) and (8, 4) the joint expectation value ⟨ni nj⟩ = 1, 
    the covariance remains zero because the individual expectation values satisfy 
    ⟨ni⟩⟨nj⟩ = 1. As a result, the correlation function evaluates to zero, 
    confirming that the state exhibits no quantum correlations.

    This test ensures that the `get_corr_function` behaves correctly for uncorrelated 
    product states and displacement vectors that span multiple qubits.
    """

    # Starting data
    R_interatomic = 1.0  
    N_SIDE = 3
    reg = Register.square(N_SIDE, spacing=R_interatomic)  
    k, l = 1, 1  

    # Create a simple final state: all qubits in the |r> state
    N = N_SIDE * N_SIDE
    final_state = qutip.tensor(*[qutip.basis(2, i % 2) for i in range(N)])  # |010101010>

    # Expected correlation value
    expected_covariance = 0

    # Calculated result
    result = get_corr_function(k, l, reg, R_interatomic, final_state, N_SIDE)
    
    
    # Assert equality
    assert np.isclose(result, expected_covariance), f"Expected {expected_covariance}, got {result}"

def test_get_corr_function_self_correlation():
    """
    Test the `get_corr_function` for self-correlation in a uniform quantum state.

    This test evaluates the correlation function for the self-correlation case 
    (k, l) = (0, 0), where the displacement vector represents no spatial separation 
    (i.e., correlations of qubits with themselves).

    The register is initialized as a 3x3 lattice, and the final state consists 
    of all qubits in the Rydberg state (|r>). For such a uniform state, the expected 
    covariance for self-correlation should be zero.

    This test ensures that the function correctly calculates the self-correlation 
    and handles the trivial displacement vector case.
    """
    # Starting data
    R_interatomic = 1.0 
    N_SIDE = 3 
    reg = Register.square(N_SIDE, spacing=R_interatomic) 
    k, l = 0, 0  

    # Create a simple final state: all qubits in the |r> state
    N_qubits = N_SIDE * N_SIDE
    final_state = qutip.tensor([qutip.basis(2, 0)] * N_qubits)

   

    expected_covariance = 0
    result = get_corr_function(k, l, reg, R_interatomic, final_state, N_SIDE)
    
    
    # Assert equality
    assert np.isclose(result, expected_covariance), f"Expected {expected_covariance}, got {result}"



####################################################################################################



def test_prepare_correlation_matrix_valid():
    """
    Test the `prepare_correlation_matrix` function with valid inputs.

    The function is tested with a valid correlation function dictionary, ensuring that
    it correctly reshapes and normalizes the matrix.
    """
    # Starting data:
    N_SIDE = 3
    correlation_function = {
        (-2, -2): 1,
        (-2, -1): 2,
        (-2, 0): 3,
        (-2, 1): 4,
        (-2, 2): 5,
        (-1, -2): 6,
        (-1, -1): 7,
        (-1, 0): 8,
        (-1, 1): 9,
        (-1, 2): 10,
        (0, -2): 11,
        (0, -1): 12,
        (0, 0): 13,
        (0, 1): 14,
        (0, 2): 15,
        (1, -2): 16,
        (1, -1): 17,
        (1, 0): 18,
        (1, 1): 19,
        (1, 2): 20,
        (2, -2): 21,
        (2, -1): 22,
        (2, 0): 23,
        (2, 1): 24,
        (2, 2): 25,
    }

    # Expected reshaped and normalized matrix
    expected_matrix = np.array([
        [0.04, 0.08, 0.12, 0.16, 0.2],
        [0.24, 0.28, 0.32, 0.36, 0.4],
        [0.44, 0.48, 0.52, 0.56, 0.6],
        [0.64, 0.68, 0.72, 0.76, 0.8],
        [0.84, 0.88, 0.92, 0.96, 1.0]
    ])

    # Call the function
    result_matrix = prepare_correlation_matrix(correlation_function, N_SIDE)

    # Assertions
    assert result_matrix.shape == (2 * N_SIDE - 1, 2 * N_SIDE - 1), \
        f"Expected shape {(2 * N_SIDE - 1, 2 * N_SIDE - 1)}, got {result_matrix.shape}"
    assert np.allclose(result_matrix, expected_matrix), \
        f"Expected matrix:\n{expected_matrix}\nGot:\n{result_matrix}"


def test_prepare_correlation_matrix_empty():
    """
    Test the `prepare_correlation_matrix` function with an empty correlation function.

    The function is expected to return an empty matrix or raise a meaningful error.
    """
    # Starting data
    N_SIDE = 3
    correlation_function = {}

    # Assert that a ValueError is raised
    with pytest.raises(ValueError, match="Correlation function is empty."):
        prepare_correlation_matrix(correlation_function, N_SIDE)


def test_prepare_correlation_matrix_unbalanced():
    """
    Test the `prepare_correlation_matrix` function with unbalanced input.

    The function is tested with a dictionary that has fewer than 
    expected for the lattice size, ensuring proper handling of incorrect inputs.
    """
    # Starting data:
    N_SIDE = 3
    correlation_function = {
        (-1, 0): 3,
        (0, 0): 4,
        (1, 0): 5  # Insufficient entries for (2 * N_SIDE - 1)^2
    }

    # Assert that a ValueError is raised
    with pytest.raises(ValueError, match="Expected .* entries, but got .*"):
        prepare_correlation_matrix(correlation_function, N_SIDE)


def test_prepare_correlation_matrix_normalization():
    """
    Test the `prepare_correlation_matrix` function to ensure proper normalization.

    The function is tested with values that have a known maximum, verifying that 
    the matrix is normalized to the range [0, 1].
    """
    # Starting data
    N_SIDE = 3
    correlation_function = {
        (k, l): k * l for k in range(-2, 3) for l in range(-2, 3)
    }

    # Call the function
    result_matrix = prepare_correlation_matrix(correlation_function, N_SIDE)

    # Check normalization
    max_value = np.max(result_matrix)
    assert np.isclose(max_value, 1), f"Expected max value 1 after normalization, got {max_value}"

####################################################################################################

def test_neel_structure_factor_positive():
    """
    Test if `get_neel_structure_factor` always returns positive values.
    
    This test verifies that the Néel structure factor is non-negative for various 
    quantum states, including uniform, alternating, and entangled states, across 
    a predefined 3x3 quantum register.
    """
    # Constants
    R_interatomic = 1.0
    N_SIDE = 3
    reg = Register.square(N_SIDE, spacing=R_interatomic)

    # Test cases
    test_states = [
        # Uniform state: All |r>=|0>
        qutip.tensor(*[qutip.basis(2, 0) for _ in range(N_SIDE * N_SIDE)]),

        # Uniform state: All |1>
        qutip.tensor(*[qutip.basis(2, 1) for _ in range(N_SIDE * N_SIDE)]),

        # Alternating state: |010101...>
        qutip.tensor(*[qutip.basis(2, i % 2) for i in range(N_SIDE * N_SIDE)]),

        # Entangled state: |010101...> + |101010...>
        (1 / np.sqrt(2)) * (
            qutip.tensor(*[qutip.basis(2, i % 2) for i in range(N_SIDE * N_SIDE)]) +
            qutip.tensor(*[qutip.basis(2, (i + 1) % 2) for i in range(N_SIDE * N_SIDE)])
        ),
    ]

    # Test each state
    for idx, state in enumerate(test_states):
        result = get_neel_structure_factor(reg, R_interatomic, state, N_SIDE)
        
        assert result >= 0, f"Néel structure factor is negative for test case {idx + 1}."
