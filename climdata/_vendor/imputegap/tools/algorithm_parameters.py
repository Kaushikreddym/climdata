import numpy as np
from skopt.space import Integer, Real
from ray import train, tune

# CDRec parameters
CDREC_RANK_RANGE = [i for i in range(1, 11)]  # This will generate a range from 1 to 10
CDREC_EPS_RANGE = np.logspace(-6, 0, num=10)  # log scale for eps
CDREC_ITERS_RANGE = [i * 100 for i in range(1, 11)]  # replace with actual range

# IIM parameters
IIM_LEARNING_NEIGHBOR_RANGE = [i for i in range(1, 100)]  # Test up to 100 learning neighbors

# MRNN parameters
MRNN_LEARNING_RATE_CHANGE = np.logspace(-6, 0, num=20)  # log scale for learning rate
MRNN_HIDDEN_DIM_RANGE = [i for i in range(0, 300, 25)]  # number of epochs
MRNN_SEQ_LEN_RANGE = [i for i in range(10)]  # sequence length
MRNN_NUM_ITER_RANGE = [i for i in range(0, 100, 5)]  # number of epochs
MRNN_KEEP_PROB_RANGE = np.logspace(-6, 0, num=10)  # dropout keep probability

# STMVL parameters
STMVL_WINDOW_SIZE_RANGE = [i for i in range(2, 100)]  # window size
STMVL_GAMMA_RANGE = np.logspace(-6, 0, num=10, endpoint=False)  # smoothing parameter gamma
STMVL_ALPHA_RANGE = [i for i in range(1, 10)]  # smoothing parameter alpha

# Define the search space for each algorithm separately
SEARCH_SPACES = {
    'cdrec': [Integer(1, 9, name='rank'), Real(1e-6, 1, "log-uniform", name='epsilon'), Integer(100, 1000, name='iteration')],
    'iim': [Integer(1, 100, name='learning_neighbors')],
    'mrnn': [Integer(10, 15, name='hidden_dim'), Real(1e-6, 1e-1, "log-uniform", name='learning_rate'), Integer(10, 95, name='iterations')],
    'stmvl': [Integer(2, 99, name='window_size'), Real(1e-6, 0.999999, "log-uniform", name='gamma'), Integer(1, 9, name='alpha')],
}

SEARCH_SPACES_PSO = {
    'cdrec': [(1, 9), (1e-6, 1), (100, 1000)],
    'iim': [(1, 100)],
    'mrnn': [(1, 15), (1e-6, 1e-1), (10, 95)],
    'stmvl': [(2, 99), (1e-6, 0.999999), (1, 9)]
}

# Define the parameter names for each algorithm
PARAM_NAMES = {
    'cdrec': ['rank', 'epsilon', 'iteration'],
    'iim': ['learning_neighbors'],
    'mrnn': ['hidden_dim', 'learning_rate', 'iterations'],
    'stmvl': ['window_size', 'gamma', 'alpha']
}


CDREC_PARAMS = {'rank': CDREC_RANK_RANGE, 'epsilon': CDREC_EPS_RANGE, 'iteration': CDREC_ITERS_RANGE}
IIM_PARAMS = {'learning_neighbors': IIM_LEARNING_NEIGHBOR_RANGE}
MRNN_PARAMS = {'learning_rate': MRNN_LEARNING_RATE_CHANGE, 'hidden_dim': MRNN_HIDDEN_DIM_RANGE, 'iterations': MRNN_NUM_ITER_RANGE}
STMVL_PARAMS = {'window_size': STMVL_WINDOW_SIZE_RANGE, 'gamma': STMVL_GAMMA_RANGE, 'alpha': STMVL_ALPHA_RANGE}

# Create a dictionary to hold all parameter dictionaries for each algorithm
ALL_ALGO_PARAMS = {'cdrec': CDREC_PARAMS, 'iim': IIM_PARAMS, 'mrnn': MRNN_PARAMS, 'stmvl': STMVL_PARAMS}


RAYTUNE_PARAMS = {
    'cdrec': {
        "rank": tune.grid_search([i for i in range(2, 16, 1)]),
        "eps": tune.loguniform(1e-6, 1),
        "iters": tune.grid_search([i * 50 for i in range(1, 4)])
    },

    "iim": {
        "learning_neighbors": tune.grid_search([i for i in range(1, 20)])  # Up to 100 learning neighbors
    },

    "mrnn": {
        "hidden_dim":  tune.grid_search([i for i in range(10, 100, 20)]),  # Hidden dimension
        "learning_rate": tune.grid_search([0.001, 0.01, 0.1, 1]),  # Log scale for learning rate
        "num_iter": tune.grid_search([i for i in range(5, 10, 5)]),  # Number of epochs
        "seq_len": 7  # tune.grid_search([i for i in range(5, 7, 2)]),  # Sequence length
        #"keep_prob": tune.loguniform(1e-6, 1)  # Dropout keep probability
    },

    "stmvl": {
        "window_size": tune.grid_search([i for i in range(10, 100, 10)]),  # Window size
        #"gamma": tune.loguniform(1e-6, 1),  # Smoothing parameter gamma
        "gamma": tune.grid_search([1e-6, 1e-3, 0.1]),
        "alpha": tune.grid_search([i for i in range(1, 6, 2)]),  # Alpha values
    },

    "iterative_svd": {
        "rank": tune.grid_search([i for i in range(2, 16, 1)])  # Testing rank from 2 to 18
    },

    "grouse": {
        "max_rank": tune.grid_search([i for i in range(2, 16, 1)])  # Testing rank from 2 to 18
    },

    "rosl": {
        "rank": tune.grid_search([i for i in range(2, 15, 2)]),  # Testing rank from 2 to 18
        "regularization": tune.grid_search([0.1, 0.2, 0.4, 0.6, 0.8]) # Regularization parameter
    },

    "soft_impute": {
        "max_rank": tune.grid_search([i for i in range(2, 16, 1)])   # Testing max_rank from 5 to 15
    },

    "spirit": {
        "k": tune.grid_search([i for i in range(1, 10, 2)]),  # Number of components
        "w": tune.grid_search([i for i in range(1, 10, 2)]),  # Window size
        "lvalue": tune.grid_search([0.1, 0.5, 1, 3, 5])  # Eigenvalue scaling
    },

    "svt": {
        "tau": tune.grid_search([0.01, 0.1, 0.2, 0.5, 1, 5]),  # Singular value thresholding parameter
        "delta": tune.grid_search([0.001, 0.01, 0.1, 1]),  # Step size for SVT
        "max_iter": tune.grid_search([i * 10 for i in range(5, 20, 5)])  # Max iterations (50 to 150)
    },

    # --- Newly Added Pattern-Based Algorithms ---

    "dynammo": {
        "h": tune.grid_search([i for i in range(3, 10)]),  # Pattern length, range 3 to 9
        "max_iteration": tune.grid_search([i for i in range(3, 15, 3)]),  # Iteration range 3 to 12
        "approximation": True
    },

    "tkcm": {
        "rank": tune.grid_search([i for i in range(2, 16, 1)])   # Testing rank from 2 to 8
    },

    # --- Newly Added Deep Learning-Based Algorithms ---

    "brits": {
        "model": tune.grid_search(["brits_i_univ"]),  # Support both univariate & multivariate models
        "epoch": tune.grid_search([i for i in range(5, 20, 5)]),  # Epochs from 5 to 15
        "batch_size": tune.grid_search([8, 16, 32]),  # Test different batch sizes
        "nbr_features": 1,  # tune.grid_search([1, 2, 5]),  # Number of features
        "hidden_layers": tune.grid_search([32, 64, 128]),
        "num_workers": 0
    },

    "deep_mvi": {
        "max_epoch": tune.grid_search([10, 50, 100]),  # Testing from 500 to 1500 epochs
        "patience": tune.grid_search([2, 5, 10]),  # Number of early stopping patience
        "lr": tune.grid_search([0.001, 0.1])  # learning rate
    },

    "mpin": {
        "incre_mode": tune.choice(["alone", "data", "state", "state+transfer", "data+state", "data+state+transfer"]),  # Different incremental modes
        "window": 1,  # Window size variations
        "k": tune.grid_search([5, 10]),  # Number of neighbors
        "learning_rate": tune.grid_search([0.01, 0.1]),  # Learning rate range
        "weight_decay": tune.grid_search([0.1, 0.5]),  # Weight decay regularization
        "epochs": 200,  # Number of epochs
        "num_of_iteration": 5,  # Number of epochs
        "threshold":  tune.grid_search([0.25, 0.75]),  # Threshold range
        "base": tune.choice(["SAGE", "GAT", "GCN"])  # Model architectures
    },

    "pristi": {
        "target_strategy": "block",  # Different strategies
        "unconditional": tune.choice([True, False]),  # Use unconditional or not
        "batch_size": tune.grid_search([-1, 8]),
        "embedding": tune.grid_search([-1, 8]),
        "num_workers": 0,
        "seed": 42
    },

    "knn_impute": {
        "k": tune.grid_search([1, 3, 5, 7, 10]),
        "weights": tune.choice(["uniform", "distance"])
    },


    "knn": {
        "k": tune.grid_search([1, 12, 1]),
        "weights": tune.choice(["uniform", "distance"])
    },

    "interpolation": {
        "method": tune.choice(["nearest", "spline", "polynomial", "linear"]),
        "poly_order": tune.grid_search([2, 10, 1])
    },

    "trmf": {
        "lags": tune.grid_search([[], [1, 2, 3], [1, 5, 10]]),  # Different lag configurations
        "K": tune.choice([-1, 5, 10, 20]),  # Latent dimensions
        "lambda_f": tune.grid_search([0.1, 1.0, 10.0]),  # Regularization parameter for factors
        "lambda_x": tune.grid_search([0.1, 1.0, 10.0]),  # Regularization parameter for observations
        "lambda_w": tune.grid_search([0.1, 1.0, 10.0]),  # Regularization parameter for weights
        "eta": tune.grid_search([0.1, 1.0, 5.0]),  # Learning rate-like parameter
        "alpha": tune.grid_search([100.0, 500.0]),  # Temporal regularization strength
        "max_iter": tune.choice([100])  # Maximum number of iterations
    },

    "mice": {
        "max_iter": tune.grid_search([2, 3]),
        "tol": tune.grid_search([0.001, 0.1]),
        "initial_strategy": tune.choice(["mean", "median", "most_frequent", "constant"]),
        "seed": 42
    },

    "miss_forest": {
        "n_estimators": tune.grid_search([2, 10, 15]),
        "max_iter": tune.grid_search([2, 5, 10]),
        "max_features": tune.choice(["auto", "sqrt", "log2"]),
        "seed": 42
    },

    "xgboost": {
        "n_estimators": tune.grid_search([1, 2, 7, 10, 15]),
        "seed": 42
    },

    "miss_net": {
        "alpha": tune.grid_search([0.1, 0.5, 1.0]),  # Example search space for alpha
        "beta": tune.grid_search([0.01, 0.1, 0.5]),  # Example search space for beta
        "L": tune.grid_search([5, 10, 20]),  # Example search space for L (hidden dimension)
        "n_cl": 1,  # Number of clusters
        "max_iter": tune.grid_search([10]),  # Max iterations
        "tol": tune.grid_search([1, 5, 10]),  # Tolerance values
        "random_init": False
    },

    "grin": {
        "d_hidden": tune.grid_search([16, 64]),  # Example search space for d_hidden
        "lr": tune.grid_search([0.001, 0.01]),  # Example search space for learning rate
        "batch_size": -1,  # Example search space for batch size
        "window": 1,  # Example search space for window size
        "alpha": tune.grid_search([10, 20]),  # Example search space for alpha
        "patience": 10,  # Patience for early stopping
        "epochs": 100,  # Max training epochs
        "workers": 4 # Number of workers
    },

    "gain": {
        "batch_size": tune.grid_search([1, 10, 32]),
        "hint_rate": tune.grid_search([0.01, 0.1, 0.5, 0.75, 0.9]),
        "alpha": tune.grid_search([1, 5, 10]),
        "epoch": 100
    },

    "bay_otide": {
        "K_trend": tune.grid_search([15, 30]),  # Trend factor search space
        "K_season": tune.grid_search([1]),  # Seasonal factor search space
        "n_season": tune.grid_search([3, 7, 10]),  # Seasonal components per factor
        "K_bias": tune.grid_search([0, 1]),  # Bias factor inclusion
        "time_scale": tune.grid_search([1, 2]),  # Time scaling factor
        "a0": tune.grid_search([0.1, 0.6, 1]),  # Prior hyperparameter a0
        "b0": tune.grid_search([1, 2.5, 5]),  # Prior hyperparameter b0
        "v": tune.grid_search([0.1, 0.5, 1.0]),  # Variance parameter
        "tr_ratio": tune.grid_search([0.8])  # tr parameter
    },

    "bit_graph": {
            "node_number": -1,  # Trend factor search space
            "kernel_set": tune.grid_search([[1,2,3,4], [2,3,6,7]]),  # Seasonal factor search space
            "dropout": tune.grid_search([0.05, 0.3]),  # Seasonal components per factor
            "subgraph_size": tune.grid_search([2, 5]),  # Bias factor inclusion
            "node_dim": 3,  # Time scaling factor
            "seq_len": -1,  # Prior hyperparameter a0
            "lr": 0.001,  # Prior hyperparameter b0
            "batch_size": tune.grid_search([8, 32]),   # Prior hyperparameter b0
            "epoch": 50,  # Variance parameter
            "num_workers": 0,  # Variance parameter
            "seed": 42 # Variance parameter
        },

    "hkmf_t": {
        "tags": [],
        "data_names": [],
        "epoch": tune.grid_search([2, 10])
    },

    "nuwats": {
        "seq_length": tune.grid_search([-1, 10]),
        "patch_size": -1,
        "batch_size": tune.grid_search([-1, 32]),
        "pred_length": -1,
        "label_length":-1,
        "enc_in": tune.grid_search([-1, 16]),
        "dec_in": tune.grid_search([-1, 16]),
        "c_out": tune.grid_search([-1, 16]),
        "gpt_layers":tune.grid_search([2, 6, 16]),
        "num_workers":0,
        "seed":42
    },
    "gpt4ts": {
        "seq_length": tune.grid_search([-1, 10]),
        "patch_size": -1,
        "batch_size": tune.grid_search([-1, 32]),
        "pred_length": -1,
        "label_length":-1,
        "enc_in": tune.grid_search([-1, 16]),
        "dec_in": tune.grid_search([-1, 16]),
        "c_out": tune.grid_search([-1, 16]),
        "gpt_layers":tune.grid_search([2, 6, 16]),
        "num_workers":0,
        "seed":42
    }

    # your ray-tune limitation #contributing
    #
    #'your_algo_name': {
    #    "param_1": tune.grid_search([i for i in range(2, 16, 1)]),
    #    "param_2": tune.loguniform(1e-6, 1),
    #    "param_3": tune.grid_search([i * 50 for i in range(1, 4)])
    #},
}
