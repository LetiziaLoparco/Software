import numpy as np
import qutip
import matplotlib.pyplot as plt
from Config import N_SIDE

def occupation(j, N_qubits):
    """
    Create the occupation operator |r><r| for each site j of the register, where the Rydberg state |r> represents the excited state.

    Parameters
    ----------
    j : int
        Index of the site for which the occupation operator is created.
    N_qubits : int
        Number of qubits in the quantum register.

    Returns
    -------
    qutip.Qobj
        Occupation operator for the specified site.
    """
    up = qutip.basis(2, 0) # creates a basis vector in the 2-dimensional Hilbert space, state index=0 it means |r> in this case.
    prod = [qutip.qeye(2) for _ in range(N_qubits)] # creates a list of 2x2 identity matrix (for N_qubits).
    prod[j] = up * up.dag() #forms a projector operator for |r> at the state j.
    return qutip.tensor(prod) # computes the tensor product of the list of operators.



def get_corr_pairs(k, l, register, R_interatomic):
    """
    Determine all couples (i,j) for a given (k,l).

    Parameters
    ----------
    k : int
        Displacement in the first dimension.
    l : int
        Displacement in the second dimension.
    register : object
        Quantum register containing qubit positions.
    R_interatomic : float
        Interatomic distance scaling factor.

    Returns
    -------
    list of list
        Pairs of indices of qubits satisfying the condition.
    """
    corr_pairs = []
    for i, qi in enumerate(register.qubits):
        for j, qj in enumerate(register.qubits):
            r_ij = register.qubits[qi] - register.qubits[qj]
            distance = np.linalg.norm(r_ij - R_interatomic * np.array([k, l]))
            if distance < 1 :
                corr_pairs.append([i, j])


    return corr_pairs



def get_corr_function(k, l, reg, R_interatomic, final_state):
    """
    Calculate the covariance in order to calculate the correlation function for specific displacements.

    Parameters
    ----------
    k : int
        Displacement in the first dimension.
    l : int
        Displacement in the second dimension.
    reg : object
        Quantum register containing qubit positions.
    R_interatomic : float
        Interatomic distance scaling factor.
    final_state : qutip.Qobj
        Final quantum state of the system.

    Returns
    -------
    float
        Correlation function for specific displacements.
    """
    
    N_qubits = N_SIDE * N_SIDE
    corr_pairs = get_corr_pairs(k, l, reg, R_interatomic)

    # Handle the case where no pairs exist
    if not corr_pairs:  
        return 0

    operators = [occupation(j, N_qubits) for j in range(N_qubits)]
    covariance = 0
    for qi, qj in corr_pairs:
        covariance += qutip.expect(operators[qi] * operators[qj], final_state)
        covariance -= qutip.expect(operators[qi], final_state) * qutip.expect(operators[qj], final_state)

    return covariance / len(corr_pairs)



def get_full_corr_function(reg, final_state, R_interatomic):
    """
    Calculate the full correlation function for all displacement vectors.

    Parameters
    ----------
    reg : object
        Quantum register containing qubit positions.
    final_state : qutip.Qobj
        Final quantum state of the system.
    R_interatomic : float
        Interatomic distance scaling factor.

    Returns
    -------
    dict
        Full correlation function with displacement vectors as keys.
    """
    correlation_function = {}
    
    for k in range(-N_SIDE + 1, N_SIDE):
        for l in range(-N_SIDE + 1, N_SIDE):
            correlation_function[(k, l)] = get_corr_function(k, l, reg, R_interatomic, final_state)
    return correlation_function



def prepare_correlation_matrix(correlation_function):
    """
    Prepare the correlation matrix for visualization.

    Parameters
    ----------
    correlation_function : dict
        Dictionary of correlation function values.

    Returns
    -------
    np.ndarray
        Normalized correlation matrix for visualization.
    """
    A = 4 * np.reshape(
        list(correlation_function.values()), (2 * N_SIDE - 1, 2 * N_SIDE - 1)
    )
    return A / np.max(A)



def create_figure_correlation_matrix(A, t_tot, ax = None):
    """
    Plot the correlation matrix.

    Parameters
    ----------
    A : np.ndarray
        Correlation matrix.
    t_tot : float
        Total time of the pulse sequence.
    ax : matplotlib.axes.Axes, optional
        Axes object for embedding the plot. If None, a new plot is created.

    Returns
    -------
    None
    """

    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 5))  # Create a new figure and subplot if none is provided

    cax = ax.imshow(A, cmap="coolwarm", vmin=-0.6, vmax=0.6)
    ax.set_xticks(range(len(A)))
    ax.set_xticklabels([f"{x}" for x in range(-N_SIDE + 1, N_SIDE)])
    ax.set_xlabel(r"$\mathscr{k}$", fontsize=22)
    ax.set_yticks(range(len(A)))
    ax.set_yticklabels([f"{-y}" for y in range(-N_SIDE + 1, N_SIDE)])
    ax.set_ylabel(r"$\mathscr{l}$", rotation=0, fontsize=22, labelpad=10)

    # Add colorbar to the figure
    if ax is None:
        fig.colorbar(cax, ax=ax, fraction=0.047, pad=0.02)
    else:
        plt.colorbar(cax, ax=ax, fraction=0.047, pad=0.02)

    ax.set_title(f"Correlation function at $t_{{tot}}={t_tot:.1f} \\mu s$", fontsize=14)

    # Show the plot only if a new figure was created
    if ax is None:
        plt.show()



def get_neel_structure_factor(reg, R_interatomic, state):
    """
    Calculate the Néel structure factor for a quantum register.

    Parameters
    ----------
    reg : object
        Quantum register containing qubit positions.
    R_interatomic : float
        Interatomic distance scaling factor.
    state : qutip.Qobj
        Final quantum state of the system.

    Returns
    -------
    float
        Calculated Néel structure factor.
    """

    st_fac = 0
    for k in range(-N_SIDE + 1, N_SIDE):
        for l in range(-N_SIDE + 1, N_SIDE):
            kk = np.abs(k)
            ll = np.abs(l)
            if not (k == 0 and l == 0):
                st_fac += (
                     (-1) ** (kk + ll) * get_corr_function(k, l, reg, R_interatomic, state)
                )
    return 4 * st_fac



def create_figure_Neel_structure_factor(t_tot, neel_structure_factors, ax=None):
    """
    Plot the Néel structure factor as a function of total time.

    Parameters
    ----------
    t_tot : np.ndarray
        Total times for each sweep.
    neel_structure_factors : np.ndarray
        Néel structure factors corresponding to each time.
    ax : matplotlib.axes.Axes, optional
        Axes object for embedding the plot. If None, a new plot is created.

    Returns
    -------
    None
    """
    # Use the provided axes or create a new figure and axes
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 5))  # Create new figure and axes if ax is not provided

    # Plot the Néel structure factor
    ax.plot(t_tot, neel_structure_factors, "o--", label="Néel Structure Factor $S_{Neel}$")
    ax.set_xlabel(r"$t_{{tot}} \, \mu s$", fontsize=14)
    ax.set_ylabel(r"$S_{Neel}$", fontsize=14)
    ax.set_title(r"Néel Structure Factor vs. $t_{{tot}}$", fontsize=16)
    ax.grid(alpha=0.5)
    ax.legend()

    # If a new figure was created, show the plot
    if ax is None:
        plt.show()




