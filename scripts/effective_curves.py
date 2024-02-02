import numpy as np
from scipy.optimize import curve_fit

def get_m_eff_log(C: np.ndarray) -> np.ndarray:
    """Effective mass curve

    log: log(C[t]/C[t+1])

    Args:
        C (np.ndarray): Correlator C(t)

    Returns:
        np.ndarray: M_eff(t) = log(C(t)/C(t+1))
    """
    T = C.shape[0] ## temporal extent
    return np.array([np.log(C[t]/C[t+1]) if C[t]/C[t+1] > 0 else 0.0 for t in range(T-1)])
####

def get_m_eff(C: np.ndarray, strategy: str) -> np.ndarray:
    """Effective mass curve from the correlator

    Args:
        C (np.ndarray): correlator C(t)
        strategy (str): computation strategy. Supported: ["log"]

    Raises:
        ValueError: if strategy is not in the list of supported types

    Returns:
        np.ndarray: array of effective mass values M_eff(t)
    """
    if strategy == "log":
        return get_m_eff_log(C)
    else:
        raise ValueError("Illegal strategy for calculation of effective mass: {strategy}".format(strategy=strategy))
    ####
####

def fit_eff_mass(m_eff: np.ndarray, dm_eff: np.ndarray) -> float:
    """Fit the effective mass (in the plateau) to a constant

    Args:
        m_eff (np.ndarray): array of values for M_eff(t) ONLY in the plateau
        dm_eff (np.ndarray): uncertainty on M_eff(t)

    Returns:
        float: best fit value of M_eff
    """
    T = m_eff.shape[0]
    t = [i for i in range(T)]
    ansatz = lambda x, m: x
    par, cov = curve_fit(ansatz, t, m_eff, sigma=dm_eff, p0=np.average(m_eff))
    return par[0]
####

## example code
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    T = 24
    t = np.array([i for i in range(T)])
    m = 0.15
    C = np.random.normal(1.0, 0.001, T)*np.exp(-m*t)

    m_eff = get_m_eff(C, strategy="log")[4:16]

    plt.plot(m_eff)
    plt.show()

    print(fit_eff_mass(m_eff=m_eff, dm_eff=0.01*m_eff))
####
