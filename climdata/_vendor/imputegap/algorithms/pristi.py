import time

from imputegap.wrapper.AlgoPython.Pristi_.runner_pristi import recov_pristi


def pristi(incomp_data, target_strategy="hybrid", unconditional=True, batch_size=-1, embedding=-1, num_workers=0, tr_ratio=0.9, seed=42, logs=True, verbose=True):
    """
    Perform imputation using the priSTI (Probabilistic Imputation via Sequential Targeted Imputation) algorithm.

    Parameters
    ----------
    incomp_data : numpy.ndarray
        The input matrix with contamination (missing values represented as NaNs).
    target_strategy : str, optional
        The strategy to use for targeting missing values. Options include: "hybrid", "random", "historical" (default is "hybrid").
    unconditional : bool, optional
        Whether to use an unconditional imputation model (default is True).
        If False, conditional imputation models are used, depending on available data patterns.
    seed : int, optional
        Random seed for reproducibility (default is 42).
    batch_size : int, optional
        Size of the batch to train the deep learning model (-1 means compute automatically based on the dataset shape).
    embedding : int, optional
        Size of the embedding used to train the deep learning model (-1 means compute automatically based on the dataset shape).
    num_workers: int, optional
         Number of worker for multiprocess (default is 0).
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
        >>> recov_data = recov_pristi(incomp_data=ts_input, target_strategy="hybrid", unconditional=True, seed=42)
        >>> print(recov_data)

    References
    ----------
    M. Liu, H. Huang, H. Feng, L. Sun, B. Du and Y. Fu, "PriSTI: A Conditional Diffusion Framework for Spatiotemporal Imputation," 2023 IEEE 39th International Conference on Data Engineering (ICDE), Anaheim, CA, USA, 2023, pp. 1927-1939, doi: 10.1109/ICDE55515.2023.00150.
    https://github.com/LMZZML/PriSTI
    """
    start_time = time.time()  # Record start time

    recov_data = recov_pristi(data=incomp_data, target_strategy=target_strategy, unconditional=unconditional, batch_size=batch_size, embedding=embedding, num_workers=num_workers, tr_ratio=tr_ratio, seed=seed, verbose=verbose)

    end_time = time.time()
    if logs and verbose:
        print(f"\n> logs: imputation priSTI - Execution Time: {(end_time - start_time):.4f} seconds\n")

    return recov_data
