import time

from imputegap.wrapper.AlgoPython.MPIN.runnerMPIN import recoverMPIN


def mpin(incomp_data=None, incre_mode="alone", window=2, k=10, lr=0.01, weight_decay=0.1, epochs=200, num_of_iteration=5, thre=0.25, base="SAGE", tr_ratio=0.9, logs=True, verbose=True):
    """
    Perform imputation using the Missing Value Imputation for Multi-attribute Sensor Data Streams via Message Propagation algorithm.

    Parameters
    ----------
    incomp_data : numpy.ndarray
        The input matrix with contamination (missing values represented as NaNs).
    incre_mode : str, optional
        The mode of incremental learning. Options are: 'alone',  'data', 'state', 'state+transfer', 'data+state', 'data+state+transfer' (default is "alone").
    window : int, optional
        The size of the sliding window for processing data streams (default is 2).
    k : int, optional
        The number of neighbors to consider during message propagation (default is 10).
    lr : float, optional
        The learning rate for optimizing the message propagation algorithm (default is 0.01).
    weight_decay : float, optional
        The weight decay (regularization) term to prevent overfitting during training (default is 0.1).
    epochs : int, optional
        The number of epochs to run the training process (default is 200).
    num_of_iteration : int, optional
        The number of iteration of the whole training (default is 5).
    thre : float, optional
        The threshold for considering a missing value as imputed (default is 0.25).
    base : str, optional
        The base model used for graph representation and message propagation. Common options include "SAGE" and "GCN" (default is "SAGE").
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
        >>> recov_data = mpin(incomp_data, incre_mode="alone", window=2, k=10, lr=0.01, weight_decay=0.1, epochs=200, thre=0.25, base="SAGE")
        >>> print(recov_data)

    References
    ----------
    Li, X., Li, H., Lu, H., Jensen, C.S., Pandey, V. & Markl, V. Missing Value Imputation for Multi-attribute Sensor Data Streams via Message Propagation (Extended Version). arXiv (2023). https://arxiv.org/abs/2311.07344
    https://github.com/XLI-2020/MPIN
    """
    start_time = time.time()  # Record start time

    recov_data = recoverMPIN(input=incomp_data, mode=incre_mode, window=window, k=k, lr=lr, weight_decay=weight_decay, epochs=epochs, num_of_iteration=num_of_iteration, thre=thre, base=base, out_channels=64, eval_ratio=0.05, state=True, verbose=verbose)

    end_time = time.time()
    if logs and verbose:
        print(f"\n> logs: imputation mpin - Execution Time: {(end_time - start_time):.4f} seconds\n")

    return recov_data
