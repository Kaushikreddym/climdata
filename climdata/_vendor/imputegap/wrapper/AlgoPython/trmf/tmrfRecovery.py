import numpy as np

from imputegap.wrapper.AlgoPython.trmf.trmf import trmf


def recoveryTRMF(data, lags=[], K=-1, lambda_f=1.0, lambda_x=1.0, lambda_w=1.0, eta=1.0, alpha=1000.0, max_iter=100):
    """Temporal Regularized Matrix Factorization : https://github.com/SemenovAlex/trmf

    Parameters
    ----------
    data : numpy.ndarray
        The input matrix with contamination (missing values represented as NaNs).
    lags : array-like, optional
        Set of lag indices to use in model.
    K : int, optional
        Length of latent embedding dimension
    lambda_f : float, optional
        Regularization parameter used for matrix F.
    lambda_x : float, optional
        Regularization parameter used for matrix X.
    lambda_w : float, optional
        Regularization parameter used for matrix W.
    alpha : float, optional
        Regularization parameter used for make the sum of lag coefficient close to 1.
        That helps to avoid big deviations when forecasting.
    eta : float, optional
        Regularization parameter used for X when undercovering autoregressive dependencies.
    max_iter : int, optional
        Number of iterations of updating matrices F, X and W.

    Returns
    -------
    numpy.ndarray
        The imputed matrix with missing values recovered.
    """


    if not lags:
        lags = list(range(1, 11))

    if K == -1:
        n = data.shape[0]
        K = n + 1

    print("(IMPUTATION) TRMF: Matrix Shape: (", data.shape[0], ", ", data.shape[1], ") for lags ", lags, ", K ", K,
          ", lambda_f ", lambda_f, " lambda_x", lambda_x, ", lambda_w ", lambda_w, ", eta ", eta, " alpha", alpha,
          ", and max_iter ", max_iter, ")...")

    incomp_data = np.copy(data)  # Copy data to avoid modifying original

    model = trmf(lags, K, lambda_f, lambda_x, lambda_w, alpha, eta, max_iter)
    model.fit(incomp_data)
    data_imputed = model.impute_missings()
    data_imputed = np.array(data_imputed)

    return data_imputed