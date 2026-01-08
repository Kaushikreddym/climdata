import time

import numpy as np

from imputegap.tools import utils
from imputegap.wrapper.AlgoPython.MissNet.recoveryMissNet import MissNet

def miss_net(incomp_data, alpha, beta, L, n_cl, max_iteration, tol, random_init, tr_ratio=0.9, logs=True, verbose=True):
    """
    Perform imputation using the Multivariate Recurrent Neural Network (MRNN) algorithm.

    Parameters
    ----------
    incomp_data : numpy.ndarray
        The input matrix with contamination (missing values represented as NaNs).
    alpha : float, optional
        Trade-off parameter controlling the contribution of contextual matrix
        and time-series. If alpha = 0, network is ignored. (default 0.5)
    beta : float, optional
        Regularization parameter for sparsity. (default 0.1)
    L : int, optional
        Hidden dimension size. (default 10)
    n_cl : int, optional
        Number of clusters. (default 1)
    max_iteration : int, optional
        Maximum number of iterations for convergence. (default 20)
    tol : float, optional
        Tolerance for early stopping criteria.  (default 5)
    random_init : bool, optional
        Whether to use random initialization for latent variables. (default False)
    tr_ratio: float, optional
        Split ratio between training and testing sets (default is 0.9).
    logs : bool, optional
        Whether to log the execution time (default is True).
    verbose : bool, optional
        Whether to display the contamination information (default is True).

    Returns
    -------
    numpy.ndarray
        The imputed matrix with missing values recovered.

    Example
    -------
        >>> recov_data = miss_net(incomp_data, alpha=0.5, beta=0.1, L=10, n_cl=1, max_iteration=20, tol=5, random_init=False)
        >>> print(recov_data)

    References
    ----------
    Kohei Obata, Koki Kawabata, Yasuko Matsubara, and Yasushi Sakurai. 2024. Mining of Switching Sparse Networks for Missing Value Imputation in Multivariate Time Series. In Proceedings of the 30th ACM SIGKDD Conference on Knowledge Discovery and Data Mining (KDD '24). Association for Computing Machinery, New York, NY, USA, 2296–2306. https://doi.org/10.1145/3637528.3671760

    """

    recov = np.copy(incomp_data)
    m_mask = np.isnan(incomp_data)

    if verbose:
        print(f"(IMPUTATION) MISS NET\n\tMatrix Shape: {incomp_data.shape[0]}, {incomp_data.shape[1]}\n\talpha: {alpha}\n\tbeta: {beta}\n\tL: {L} \n\tn_cl: {n_cl}\n\tmax_iteration: {max_iteration}\n\ttr_ratio: {tr_ratio}\n")

    cont_data_matrix, mask_train, mask_test, mask_val, error = utils.dl_integration_transformation(incomp_data, tr_ratio=tr_ratio, inside_tr_cont_ratio=0.4, split_ts=1, split_val=0, nan_val=None, prevent_leak=False, block_selection=False, offset=0.05, seed=42, verbose=False)
    if error:
        return incomp_data

    start_time = time.time()  # Record start time

    missnet_model = MissNet(alpha=alpha, beta=beta, L=L, n_cl=n_cl)
    missnet_model.fit(cont_data_matrix, random_init=random_init, max_iter=max_iteration, tol=tol, verbose=verbose)  # Train the model
    recov_data = missnet_model.imputation()  # Get the imputed data

    end_time = time.time()

    recov[m_mask] = recov_data[m_mask]

    if logs and verbose:
        print(f"\n> logs: imputation miss_net - Execution Time: {(end_time - start_time):.4f} seconds\n")

    return recov
