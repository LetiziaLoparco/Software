import numpy as np
import qutip

from pulser import Register
from Config import N_SIDE


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

    If the index `j` is greater than or equal to the number of qubits or less than 0,
    the function is expected to raise an `IndexError`.
    """
    # Invalid index
    j = 6  # Out of range
    N_qubits = 5  # Total number of qubits

    # Assert that an IndexError is raised
    try:
        occupation(j, N_qubits)
    except IndexError:
        pass
    else:
        assert False, "Expected IndexError for invalid index."


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
    R_interatomic = 1.0  # Scaling factor for interatomic distance
    reg = Register.rectangle(2, 2, spacing=R_interatomic)  # Create a 2x2 lattice
    k, l = 1, 0  # Displacement vector

    # Expected pairs: [(1,0), (3,2)]
    expected_pairs = [[1,0], [3,2]]

    # Calculated result
    result_pairs = get_corr_pairs(k, l, reg, R_interatomic)

    # Assert that the result matches the expected pairs
    assert sorted(result_pairs) == sorted(expected_pairs), f"Expected {expected_pairs}, got {result_pairs}"


def test_get_corr_pairs_no_pairs():
    """
    Test the `get_corr_pairs` function when no pairs satisfy the distance condition.

    The function is tested with a register and displacement values where no qubits
    should meet the distance requirement.
    """
    # Starting data
    R_interatomic = 2.0  # Large interatomic distance
    reg = Register.square(3, spacing=R_interatomic)  # Create a 3x3 lattice
    k, l = 10, 10  # Large displacement vector

    # Expected result: no pairs
    expected_pairs = []

    # Calculated result
    result_pairs = get_corr_pairs(k, l, reg, R_interatomic)

    # Assert that the result is an empty list
    assert result_pairs == expected_pairs, f"Expected no pairs, got {result_pairs}"


def test_get_corr_pairs_single_qubit():
    """
    Test the `get_corr_pairs` function with a single-qubit register.

    The function is tested with a single qubit in the register, where no pairs
    can be formed regardless of the displacement.
    """
    # Starting data
    R_interatomic = 1.0
    reg = Register.from_coordinates([[0, 0]])  # Single qubit at the origin
    k, l = 1, 0  # Displacement vector

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
    reg = Register.square(3, spacing=R_interatomic)  # Create a 3x3 lattice
    displacements = [(1, 0), (0, 1), (1, 1)]  # Different (k, l) values

    # Expected results for each displacement
    expected_results = {
        (1, 0): [[1, 0], [2, 1], [4, 3], [5, 4], [7, 6], [8, 7]],
        (0, 1): [[3, 0], [4, 1], [5, 2], [6, 3], [7, 4], [8, 5]],
        (1, 1): [[4, 0], [5, 1], [7, 3], [8, 4]],
    }

    for (k, l), expected_pairs in expected_results.items():
        # Calculated result for each displacement
        result_pairs = get_corr_pairs(k, l, reg, R_interatomic)
        assert sorted(result_pairs) == sorted(expected_pairs), f"Displacement ({k}, {l}): Expected {expected_pairs}, got {result_pairs}"




####################################################################################################




def test_get_corr_function_valid():
    """
    Test the `get_corr_function` for valid inputs and a simple final state.

    The function is tested with a predefined register and a uniform final state. The 
    correlation function is calculated for a valid (k, l) displacement.
    """
    # Starting data
    R_interatomic = 1.0  # Interatomic distance scaling factor
    reg = Register.square(N_SIDE, spacing=R_interatomic)  # Create a 3x3 lattice
    k, l = 1, 0  # Displacement vector

    # Create a simple final state: all qubits in the |0> state
    N_qubits = N_SIDE * N_SIDE
    final_state = qutip.tensor([qutip.basis(2, 0)] * N_qubits)

    # Expected correlation value: the theoretical calculation can be added if known
    corr_pairs = get_corr_pairs(k, l, reg, R_interatomic)
    operators = [occupation(j, N_qubits) for j in range(N_qubits)]
    expected_covariance = sum(
        qutip.expect(operators[qi] * operators[qj], final_state)
        - qutip.expect(operators[qi], final_state) * qutip.expect(operators[qj], final_state)
        for qi, qj in corr_pairs
    ) / len(corr_pairs)

    # Calculated result
    result = get_corr_function(k, l, reg, R_interatomic, final_state)

    # Assert equality
    assert np.isclose(result, expected_covariance), f"Expected {expected_covariance}, got {result}"


def test_get_corr_function_no_pairs():
    """
    Test the `get_corr_function` when no valid pairs exist for the given (k, l).

    The function is tested with a displacement vector that exceeds the range of 
    the register. The result should be zero or raise a controlled error.
    """
    # Starting data
    R_interatomic = 1.0
    reg = Register.square(N_SIDE, spacing=R_interatomic)  # Create a 3x3 lattice
    k, l = 10, 10  # Large displacement vector

    # Create a simple final state
    N_qubits = N_SIDE * N_SIDE
    final_state = qutip.tensor([qutip.basis(2, 0)] * N_qubits)

    # Expected result: No valid pairs -> covariance should be zero
    expected_result = 0

    # Calculated result
    result = get_corr_function(k, l, reg, R_interatomic, final_state)

    # Assert equality
    assert np.isclose(result, expected_result), f"Expected {expected_result}, got {result}"


def test_get_corr_function_negative_displacement():
    """
    Test the `get_corr_function` for negative displacement values.

    The function is tested with negative (k, l) values to ensure that the
    calculation correctly accounts for direction.
    """
    # Starting data
    R_interatomic = 1.0
    reg = Register.square(N_SIDE, spacing=R_interatomic)  # Create a 3x3 lattice
    k, l = -1, 0  # Negative displacement vector

    # Create a simple final state
    N_qubits = N_SIDE * N_SIDE
    final_state = qutip.tensor([qutip.basis(2, 0)] * N_qubits)

    # Expected correlation value
    corr_pairs = get_corr_pairs(k, l, reg, R_interatomic)
    operators = [occupation(j, N_qubits) for j in range(N_qubits)]
    expected_covariance = sum(
        qutip.expect(operators[qi] * operators[qj], final_state)
        - qutip.expect(operators[qi], final_state) * qutip.expect(operators[qj], final_state)
        for qi, qj in corr_pairs
    ) / len(corr_pairs)

    # Calculated result
    result = get_corr_function(k, l, reg, R_interatomic, final_state)

    # Assert equality
    assert np.isclose(result, expected_covariance), f"Expected {expected_covariance}, got {result}"


####################################################################################################



def test_prepare_correlation_matrix_valid():
    """
    Test the `prepare_correlation_matrix` function with valid inputs.

    The function is tested with a valid correlation function dictionary, ensuring that
    it correctly reshapes and normalizes the matrix.
    """
    # Starting data: A mock correlation function
    N_SIDE = 3  # Example for a 3x3 lattice
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
    result_matrix = prepare_correlation_matrix(correlation_function)

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
    correlation_function = {}

    # Call the function and handle empty input
    try:
        result_matrix = prepare_correlation_matrix(correlation_function)
        assert result_matrix.size == 0, "Expected an empty matrix for empty input."
    except ValueError:
        pass  # If the function raises an error for empty input, that's acceptable.


def test_prepare_correlation_matrix_unbalanced():
    """
    Test the `prepare_correlation_matrix` function with unbalanced input.

    The function is tested with a dictionary that has fewer or more entries than 
    expected for the lattice size, ensuring proper handling of incorrect inputs.
    """
    # Starting data: Missing or extra entries
    N_SIDE = 3
    correlation_function = {
        (-1, 0): 3,
        (0, 0): 4,
        (1, 0): 5  # Insufficient entries for (2 * N_SIDE - 1)^2
    }

    # Call the function and handle incorrect input size
    try:
        prepare_correlation_matrix(correlation_function)
        assert False, "Expected an error for unbalanced input."
    except ValueError:
        pass  # Properly raises an error


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
    result_matrix = prepare_correlation_matrix(correlation_function)

    # Check normalization
    max_value = np.max(result_matrix)
    assert np.isclose(max_value, 1), f"Expected max value 1 after normalization, got {max_value}"

####################################################################################################

def test_get_neel_structure_factor_valid():
    """
    Test `get_neel_structure_factor` with a valid quantum register and state.

    The function is tested with a simple 3x3 lattice and a predefined quantum state.
    """
    # Starting data
    R_interatomic = 1.0
    reg = Register.square(N_SIDE, spacing=R_interatomic)  # Create a 3x3 lattice
    N_qubits = N_SIDE * N_SIDE

    # Simple final state: all qubits in the |0> state
    state = qutip.tensor([qutip.basis(2, 0)] * N_qubits)

    # Mock get_corr_function to return a constant value (e.g., 1 for all (k, l))
    def mock_get_corr_function(k, l, reg, R_interatomic, state):
        return 1

    # Inject the mock function
    global get_corr_function
    original_get_corr_function = get_corr_function
    get_corr_function = mock_get_corr_function

    try:
        # Calculated result
        result = get_neel_structure_factor(reg, R_interatomic, state)

        # Expected Néel structure factor
        expected_result = sum(
            (-1) ** (np.abs(k) + np.abs(l)) for k in range(-N_SIDE + 1, N_SIDE) for l in range(-N_SIDE + 1, N_SIDE)
            if not (k == 0 and l == 0)
        ) * 4

        # Assert equality
        assert np.isclose(result, expected_result), f"Expected {expected_result}, got {result}"
    finally:
        # Restore the original function
        get_corr_function = original_get_corr_function


def test_get_neel_structure_factor_zero_state():
    """
    Test `get_neel_structure_factor` with a quantum state that leads to zero correlation.

    The function is tested with a state designed to produce zero correlation values.
    """
    # Starting data
    R_interatomic = 1.0
    reg = Register.square(N_SIDE, spacing=R_interatomic)  # Create a 3x3 lattice
    N_qubits = N_SIDE * N_SIDE

    # Simple final state: all qubits in the |1> state
    state = qutip.tensor([qutip.basis(2, 1)] * N_qubits)

    # Mock get_corr_function to return zero for all (k, l)
    def mock_get_corr_function(k, l, reg, R_interatomic, state):
        return 0

    # Inject the mock function
    global get_corr_function
    original_get_corr_function = get_corr_function
    get_corr_function = mock_get_corr_function

    try:
        # Calculated result
        result = get_neel_structure_factor(reg, R_interatomic, state)

        # Expected result: Zero Néel structure factor
        expected_result = 0

        # Assert equality
        assert np.isclose(result, expected_result), f"Expected {expected_result}, got {result}"
    finally:
        # Restore the original function
        get_corr_function = original_get_corr_function

