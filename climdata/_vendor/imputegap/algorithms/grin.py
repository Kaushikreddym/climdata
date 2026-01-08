import time

from imputegap.wrapper.AlgoPython.GRIN.recoveryGRIN import recoveryGRIN


def grin(incomp_data, d_hidden=32, lr=0.001, batch_size=32, window=10, alpha=10.0, patience=4, epochs=20, workers=2, tr_ratio=0.9, logs=True, verbose=True):
    """
    Perform imputation using the Multivariate Recurrent Neural Network (MRNN) algorithm.

    Parameters
    ----------
    incomp_data : numpy.ndarray
        The input matrix with contamination (missing values represented as NaNs).

    d_hidden : int, optional, default=32
        The number of hidden units in the model's recurrent and graph layers.

    lr : float, optional, default=0.001
        Learning rate for the optimizer.

    batch_size : int, optional, default=32
        The number of samples per training batch.

    window : int, optional, default=10
        The size of the time window used for modeling temporal dependencies.

    alpha : float, optional, default=10.0
        The weight assigned to the adversarial loss term during training.

    patience : int, optional, default=4
        Number of epochs without improvement before early stopping is triggered.

    epochs : int, optional, default=20
        The maximum number of training epochs.

    workers : int, optional, default=2
        The number of worker processes for data loading.

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
        >>> recov_data = grin(incomp_data, d_hidden=32, lr=0.001, batch_size=32, window=10, alpha=10.0, patience=4, epochs=20, workers=2)
        >>> print(recov_data)

    References
    ----------
    A. Cini, I. Marisca, and C. Alippi, "Multivariate Time Series Imputation by Graph Neural Networks," CoRR, vol. abs/2108.00298, 2021
    https://github.com/Graph-Machine-Learning-Group/grin
    """
    start_time = time.time()  # Record start time

    recov_data = recoveryGRIN(input=incomp_data, d_hidden=d_hidden, lr=lr, batch_size=batch_size, window=window, alpha=alpha, patience=patience, epochs=epochs, workers=workers, tr_ratio=tr_ratio, verbose=verbose)

    end_time = time.time()
    if logs and verbose:
        print(f"\n> logs: imputation grin - Execution Time: {(end_time - start_time):.4f} seconds\n")

    return recov_data
