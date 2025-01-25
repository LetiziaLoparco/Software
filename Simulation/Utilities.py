import numpy as np
import qutip
import matplotlib.pyplot as plt
from Config import N_SIDE

def occupation(j, N):
    """
    Create the occupation operator for a single site in a quantum register.

    Parameters
    ----------
    j : int
        Index of the site for which the occupation operator is created.
    N : int
        Size of the quantum register.

    Returns
    -------
    qutip.Qobj
        Occupation operator for the specified site.
    """
    up = qutip.basis(2, 0)
    prod = [qutip.qeye(2) for _ in range(N)]
    prod[j] = up * up.dag()
    return qutip.tensor(prod)

def get_corr_pairs(k, l, register, R_interatomic):
    """
    Determine all pairs of qubits satisfying a specific distance condition.

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
            if distance < 1:
                corr_pairs.append([i, j])
    return corr_pairs

def get_corr_function(k, l, reg, R_interatomic, state):
    """
    Calculate the covariance of the correlation function for specific displacements.

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
    state : qutip.Qobj
        Quantum state of the system.

    Returns
    -------
    float
        Normalized covariance of the correlation function.
    """
    N_qubits = len(reg.qubits)
    corr_pairs = get_corr_pairs(k, l, reg, R_interatomic)

    operators = [occupation(j, N_qubits) for j in range(N_qubits)]
    covariance = 0
    for qi, qj in corr_pairs:
        covariance += qutip.expect(operators[qi] * operators[qj], state)
        covariance -= qutip.expect(operators[qi], state) * qutip.expect(operators[qj], state)
    return covariance / len(corr_pairs)

def get_full_corr_function(reg, state, R_interatomic):
    """
    Calculate the full correlation function for all displacement vectors.

    Parameters
    ----------
    reg : object
        Quantum register containing qubit positions.
    state : qutip.Qobj
        Quantum state of the system.
    R_interatomic : float
        Interatomic distance scaling factor.

    Returns
    -------
    dict
        Full correlation function with displacement vectors as keys.
    """
    N_qubits = len(reg.qubits)

    correlation_function = {}
    N_side = int(np.sqrt(N_qubits))
    for k in range(-N_side + 1, N_side):
        for l in range(-N_side + 1, N_side):
            correlation_function[(k, l)] = get_corr_function(k, l, reg, R_interatomic, state)
    return correlation_function

def prepare_correlation_matrix(correlation_function):
    """
    Prepare the correlation matrix for visualization.

    Parameters
    ----------
    correlation_function : dict
        Dictionary of correlation function values.
    N_side : int
        Number of sites per side in the lattice.

    Returns
    -------
    np.ndarray
        Normalized correlation matrix.
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
    N_side : int
        Number of sites per side in the lattice.
    ax : matplotlib.axes.Axes, optional
        Axes object for embedding the plot (default: None).

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
        Quantum state of the system.

    Returns
    -------
    float
        Néel structure factor.
    """
    N_qubits = len(reg.qubits)
    N_side = int(np.sqrt(N_qubits))

    st_fac = 0
    for k in range(-N_side + 1, N_side):
        for l in range(-N_side + 1, N_side):
            kk = np.abs(k)
            ll = np.abs(l)
            if not (k == 0 and l == 0):
                st_fac += (
                     (-1) ** (kk + ll) * get_corr_function(k, l, reg, R_interatomic, state)
                )
    return 4 * st_fac

def Create_figure_Neel_structure_factor(t_tot, neel_structure_factors, ax=None):
    """
    Plot the Néel structure factor as a function of total time.

    Parameters
    ----------
    t_tot : np.ndarray
        Total times for each sweep.
    neel_structure_factors : np.ndarray
        Néel structure factors corresponding to each time.
    ax : matplotlib.axes.Axes, optional
        Axes object for embedding the plot (default: None).

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




