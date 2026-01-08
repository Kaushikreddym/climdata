import ctypes
import os
import toml
import importlib.resources
import numpy as __numpy_import
import platform


def load_parameters(query: str = "default", algorithm: str = "cdrec", dataset: str = "chlorine", optimizer: str = "b", path=None, verbose=False):
    """
    Load default or optimal parameters for algorithms from a TOML file.

    Parameters
    ----------
    query : str, optional
        'default' or 'optimal' to load default or optimal parameters (default is "default").
    algorithm : str, optional
        Algorithm to load parameters for (default is "cdrec").
    dataset : str, optional
        Name of the dataset (default is "chlorine").
    optimizer : str, optional
        Optimizer type for optimal parameters (default is "b").
    path : str, optional
        Custom file path for the TOML file (default is None).
    verbose : bool, optional
        Whether to display the contamination information (default is False).

    Returns
    -------
    tuple
        A tuple containing the loaded parameters for the given algorithm.
    """
    if query == "default":
        if path is None:
            filepath = importlib.resources.files('imputegap.env').joinpath("./default_values.toml")
            if not filepath.is_file():
                filepath = "./env/default_values.toml"
        else:
            filepath = path
            if not os.path.exists(filepath):
                filepath = "./env/default_values.toml"

    elif query == "optimal":
        if path is None:
            filename = "./optimal_parameters_" + str(optimizer) + "_" + str(dataset) + "_" + str(algorithm) + ".toml"
            filepath = importlib.resources.files('imputegap.params').joinpath(filename)
            if not filepath.is_file():
                filepath = "./params/optimal_parameters_" + str(optimizer) + "_" + str(dataset) + "_" + str(algorithm) + ".toml"
        else:
            filepath = path
            if not os.path.exists(filepath):
                filepath = "./params/optimal_parameters_" + str(optimizer) + "_" + str(dataset) + "_" + str(algorithm) + ".toml"

    else:
        filepath = None
        print("Query not found for this function ('optimal' or 'default')")

    if not os.path.exists(filepath):
        filepath = "./params/optimal_parameters_" + str(optimizer) + "_" + str(dataset) + "_" + str(algorithm) + ".toml"
        if not os.path.exists(filepath):
            filepath = filepath[1:]

    with open(filepath, "r") as _:
        config = toml.load(filepath)

    if verbose:
        print("\n(SYS) Inner files loaded : ", filepath, "\n")

    if algorithm == "cdrec":
        truncation_rank = int(config[algorithm]['rank'])
        epsilon = float(config[algorithm]['epsilon'])
        iterations = int(config[algorithm]['iteration'])
        return (truncation_rank, epsilon, iterations)
    elif algorithm == "stmvl":
        window_size = int(config[algorithm]['window_size'])
        gamma = float(config[algorithm]['gamma'])
        alpha = int(config[algorithm]['alpha'])
        return (window_size, gamma, alpha)
    elif algorithm == "iim":
        learning_neighbors = int(config[algorithm]['learning_neighbors'])
        if query == "default":
            algo_code = config[algorithm]['algorithm_code']
            return (learning_neighbors, algo_code)
        else:
            return (learning_neighbors,)
    elif algorithm == "mrnn":
        hidden_dim = int(config[algorithm]['hidden_dim'])
        learning_rate = float(config[algorithm]['learning_rate'])
        iterations = int(config[algorithm]['iterations'])
        if query == "default":
            sequence_length = int(config[algorithm]['sequence_length'])
            return (hidden_dim, learning_rate, iterations, sequence_length)
        else:
            return (hidden_dim, learning_rate, iterations)
    elif algorithm == "iterative_svd":
        truncation_rank = int(config[algorithm]['rank'])
        return (truncation_rank)
    elif algorithm == "grouse":
        max_rank = int(config[algorithm]['max_rank'])
        return (max_rank)
    elif algorithm == "dynammo":
        h = int(config[algorithm]['h'])
        max_iteration = int(config[algorithm]['max_iteration'])
        approximation = bool(config[algorithm]['approximation'])
        return (h, max_iteration, approximation)
    elif algorithm == "rosl":
        rank = int(config[algorithm]['rank'])
        regularization = float(config[algorithm]['regularization'])
        return (rank, regularization)
    elif algorithm == "soft_impute":
        max_rank = int(config[algorithm]['max_rank'])
        return (max_rank)
    elif algorithm == "spirit":
        k = int(config[algorithm]['k'])
        w = int(config[algorithm]['w'])
        lvalue = float(config[algorithm]['lvalue'])
        return (k, w, lvalue)
    elif algorithm == "svt":
        tau = float(config[algorithm]['tau'])
        return (tau)
    elif algorithm == "tkcm":
        rank = int(config[algorithm]['rank'])
        return (rank)
    elif algorithm == "deep_mvi":
        max_epoch = int(config[algorithm]['max_epoch'])
        patience = int(config[algorithm]['patience'])
        lr = float(config[algorithm]['lr'])
        return (max_epoch, patience, lr)
    elif algorithm == "brits":
        model = str(config[algorithm]['model'])
        epoch = int(config[algorithm]['epoch'])
        batch_size = int(config[algorithm]['batch_size'])
        nbr_features = int(config[algorithm]['nbr_features'])
        hidden_layers = int(config[algorithm]['hidden_layers'])
        num_workers = int(config[algorithm]['num_workers'])
        return (model, epoch, batch_size, nbr_features, hidden_layers, num_workers)
    elif algorithm == "mpin":
        incre_mode = str(config[algorithm]['incre_mode'])
        window = int(config[algorithm]['window'])
        k = int(config[algorithm]['k'])
        learning_rate = float(config[algorithm]['learning_rate'])
        weight_decay = float(config[algorithm]['weight_decay'])
        epochs = int(config[algorithm]['epochs'])
        num_of_iteration = int(config[algorithm]['num_of_iteration'])
        threshold = float(config[algorithm]['threshold'])
        base = str(config[algorithm]['base'])
        return (incre_mode, window, k, learning_rate, weight_decay, epochs, num_of_iteration, threshold, base)
    elif algorithm == "pristi":
        target_strategy = str(config[algorithm]['target_strategy'])
        unconditional = bool(config[algorithm]['unconditional'])
        batch_size = int(config[algorithm]['batch_size'])
        embedding = int(config[algorithm]['embedding'])
        num_workers = int(config[algorithm]['num_workers'])
        seed = int(config[algorithm]['seed'])
        return (target_strategy, unconditional, batch_size, embedding, num_workers, seed)
    elif algorithm == "knn" or algorithm == "knn_impute":
        k = int(config[algorithm]['k'])
        weights = str(config[algorithm]['weights'])
        return (k, weights)
    elif algorithm == "interpolation":
        method = str(config[algorithm]['method'])
        poly_order = int(config[algorithm]['poly_order'])
        return (method, poly_order)
    elif algorithm == "trmf":
        lags = list(config[algorithm]['lags'])
        K = int(config[algorithm]['K'])
        lambda_f = float(config[algorithm]['lambda_f'])
        lambda_x = float(config[algorithm]['lambda_x'])
        lambda_w = float(config[algorithm]['lambda_w'])
        eta = float(config[algorithm]['eta'])
        alpha = float(config[algorithm]['alpha'])
        max_iter = int(config[algorithm]['max_iter'])
        return (lags, K, lambda_f, lambda_x, lambda_w, eta, alpha, max_iter)
    elif algorithm == "mice":
        max_iter = int(config[algorithm]['max_iter'])
        tol = float(config[algorithm]['tol'])
        initial_strategy = str(config[algorithm]['initial_strategy'])
        seed = int(config[algorithm]['seed'])
        return (max_iter, tol, initial_strategy, seed)
    elif algorithm == "miss_forest":
        n_estimators = int(config[algorithm]['n_estimators'])
        max_iter = int(config[algorithm]['max_iter'])
        max_features = str(config[algorithm]['max_features'])
        seed = int(config[algorithm]['seed'])
        return (n_estimators, max_iter, max_features, seed)
    elif algorithm == "xgboost":
        n_estimators = int(config[algorithm]['n_estimators'])
        seed = int(config[algorithm]['seed'])
        return (n_estimators, seed)
    elif algorithm == "miss_net":
        alpha = float(config[algorithm]['alpha'])
        beta = float(config[algorithm]['beta'])
        L = int(config[algorithm]['L'])
        n_cl = int(config[algorithm]['n_cl'])
        max_iter = int(config[algorithm]['max_iter'])
        tol = float(config[algorithm]['tol'])
        random_init = bool(config[algorithm]['random_init'])
        return (alpha, beta, L, n_cl, max_iter, tol, random_init)
    elif algorithm == "gain":
        batch_size = int(config[algorithm]['batch_size'])
        hint_rate = float(config[algorithm]['hint_rate'])
        alpha = int(config[algorithm]['alpha'])
        epoch = int(config[algorithm]['epoch'])
        return (batch_size, hint_rate, alpha, epoch)
    elif algorithm == "grin":
        d_hidden = int(config[algorithm]['d_hidden'])
        lr = float(config[algorithm]['lr'])
        batch_size = int(config[algorithm]['batch_size'])
        window = int(config[algorithm]['window'])
        alpha = int(config[algorithm]['alpha'])
        patience = int(config[algorithm]['patience'])
        epochs = int(config[algorithm]['epochs'])
        workers = int(config[algorithm]['workers'])
        return (d_hidden, lr, batch_size, window, alpha, patience, epochs, workers)
    elif algorithm == "bay_otide":
        K_trend = int(config[algorithm]['K_trend'])
        K_season = int(config[algorithm]['K_season'])
        n_season = int(config[algorithm]['n_season'])
        K_bias = int(config[algorithm]['K_bias'])
        time_scale = int(config[algorithm]['time_scale'])
        a0 = float(config[algorithm]['a0'])
        b0 = float(config[algorithm]['b0'])
        v = float(config[algorithm]['v'])
        num_workers = int(config[algorithm]['num_workers'])
        tr_ratio = float(config[algorithm]['tr_ratio'])
        return (K_trend, K_season, n_season, K_bias, time_scale, a0, b0, v, num_workers, tr_ratio)
    elif algorithm == "hkmf_t":
        tags = config[algorithm]['tags']
        data_names = config[algorithm]['data_names']
        epoch = int(config[algorithm]['epoch'])
        return (tags, data_names, epoch)
    elif algorithm == "nuwats":
        seq_length = int(config[algorithm]['seq_length'])
        patch_size = int(config[algorithm]['patch_size'])
        batch_size = int(config[algorithm]['batch_size'])
        pred_length = int(config[algorithm]['pred_length'])
        label_length = int(config[algorithm]['label_length'])
        enc_in = int(config[algorithm]['enc_in'])
        dec_in = int(config[algorithm]['dec_in'])
        c_out = int(config[algorithm]['c_out'])
        gpt_layers = int(config[algorithm]['gpt_layers'])
        num_workers = int(config[algorithm]['num_workers'])
        seed = int(config[algorithm]['seed'])
        return (seq_length, patch_size, batch_size, pred_length, label_length, enc_in, dec_in, c_out, gpt_layers, num_workers, seed)
    elif algorithm == "gpt4ts":
        seq_length = int(config[algorithm]['seq_length'])
        patch_size = int(config[algorithm]['patch_size'])
        batch_size = int(config[algorithm]['batch_size'])
        pred_length = int(config[algorithm]['pred_length'])
        label_length = int(config[algorithm]['label_length'])
        enc_in = int(config[algorithm]['enc_in'])
        dec_in = int(config[algorithm]['dec_in'])
        c_out = int(config[algorithm]['c_out'])
        gpt_layers = int(config[algorithm]['gpt_layers'])
        num_workers = int(config[algorithm]['num_workers'])
        seed = int(config[algorithm]['seed'])
        return (seq_length, patch_size, batch_size, pred_length, label_length, enc_in, dec_in, c_out, gpt_layers, num_workers, seed)
    elif algorithm == "bit_graph":
        node_number = int(config[algorithm]['node_number'])
        kernel_set = config[algorithm]['kernel_set']
        dropout = float(config[algorithm]['dropout'])
        subgraph_size = int(config[algorithm]['subgraph_size'])
        node_dim = int(config[algorithm]['node_dim'])
        seq_len = int(config[algorithm]['seq_len'])
        lr = float(config[algorithm]['lr'])
        batch_size = int(config[algorithm]['batch_size'])
        epoch = int(config[algorithm]['epoch'])
        num_workers = int(config[algorithm]['num_workers'])
        seed = int(config[algorithm]['seed'])
        return (node_number, kernel_set, dropout, subgraph_size, node_dim, seq_len, lr, batch_size, epoch, num_workers, seed)
    elif algorithm == "greedy":
        n_calls = int(config[algorithm]['n_calls'])
        metrics = config[algorithm]['metrics']
        return (n_calls, [metrics])
    elif algorithm.lower() in ["bayesian", "bo", "bayesopt"]:
        n_calls = int(config['bayesian']['n_calls'])
        n_random_starts = int(config['bayesian']['n_random_starts'])
        acq_func = str(config['bayesian']['acq_func'])
        metrics = config['bayesian']['metrics']
        return (n_calls, n_random_starts, acq_func, [metrics])
    elif algorithm.lower() in ['pso', "particle_swarm"]:
        n_particles = int(config['pso']['n_particles'])
        c1 = float(config['pso']['c1'])
        c2 = float(config['pso']['c2'])
        w = float(config['pso']['w'])
        iterations = int(config['pso']['iterations'])
        n_processes = int(config['pso']['n_processes'])
        metrics = config['pso']['metrics']
        return (n_particles, c1, c2, w, iterations, n_processes, [metrics])
    elif algorithm.lower() in  ['sh', "successive_halving"]:
        num_configs = int(config['sh']['num_configs'])
        num_iterations = int(config['sh']['num_iterations'])
        reduction_factor = int(config['sh']['reduction_factor'])
        metrics = config['sh']['metrics']
        return (num_configs, num_iterations, reduction_factor, [metrics])
    elif algorithm.lower() in ['ray_tune', "ray"]:
        metrics = config['ray_tune']['metrics']
        n_calls = int(config['ray_tune']['n_calls'])
        max_concurrent_trials = int(config['ray_tune']['max_concurrent_trials'])
        return ([metrics], n_calls, max_concurrent_trials)
    elif algorithm == "forecaster-naive":
        strategy = str(config[algorithm]['strategy'])
        window_length = int(config[algorithm]['window_length'])
        sp = int(config[algorithm]['sp'])
        return {"strategy": strategy, "window_length": window_length, "sp": sp}
    elif algorithm == "forecaster-exp-smoothing":
        trend = str(config[algorithm]['trend'])
        seasonal = str(config[algorithm]['seasonal'])
        sp = int(config[algorithm]['sp'])
        return {"trend": trend, "seasonal": seasonal, "sp": sp}
    elif algorithm == "forecaster-prophet":
        seasonality_mode = str(config[algorithm]['seasonality_mode'])
        n_changepoints = int(config[algorithm]['n_changepoints'])
        return {"seasonality_mode": seasonality_mode, "n_changepoints": n_changepoints}
    elif algorithm == "forecaster-nbeats":
        input_chunk_length = int(config[algorithm]['input_chunk_length'])
        output_chunk_length = int(config[algorithm]['output_chunk_length'])
        num_blocks = int(config[algorithm]['num_blocks'])
        layer_widths = int(config[algorithm]['layer_widths'])
        random_state = int(config[algorithm]['random_state'])
        n_epochs = int(config[algorithm]['n_epochs'])
        pl_trainer_kwargs = str(config[algorithm]['pl_trainer_kwargs'])
        if pl_trainer_kwargs == "cpu":
            drive = {"accelerator": pl_trainer_kwargs}
        else:
            drive = {"accelerator": pl_trainer_kwargs, "devices": [0]}
        return {"input_chunk_length": input_chunk_length, "output_chunk_length": output_chunk_length, "num_blocks": num_blocks,
                "layer_widths": layer_widths, "random_state": random_state, "n_epochs": n_epochs, "pl_trainer_kwargs": drive}
    elif algorithm == "forecaster-xgboost":
        lags = int(config[algorithm]['lags'])
        return {"lags": lags}
    elif algorithm == "forecaster-lightgbm":
        lags = int(config[algorithm]['lags'])
        verbose = int(config[algorithm]['verbose'])
        return {"lags": lags, "verbose": verbose}
    elif algorithm == "forecaster-lstm":
        input_chunk_length = int(config[algorithm]['input_chunk_length'])
        model = str(config[algorithm]['model'])
        random_state = int(config[algorithm]['random_state'])
        n_epochs = int(config[algorithm]['n_epochs'])
        pl_trainer_kwargs = str(config[algorithm]['pl_trainer_kwargs'])
        if pl_trainer_kwargs == "cpu":
            drive = {"accelerator": pl_trainer_kwargs}
        else:
            drive = {"accelerator": pl_trainer_kwargs, "devices": [0]}
        return {"input_chunk_length": input_chunk_length, "model": model, "random_state": random_state, "n_epochs": n_epochs, "pl_trainer_kwargs": drive}
    elif algorithm == "forecaster-deepar":
        input_chunk_length = int(config[algorithm]['input_chunk_length'])
        model = str(config[algorithm]['model'])
        random_state = int(config[algorithm]['random_state'])
        n_epochs = int(config[algorithm]['n_epochs'])
        pl_trainer_kwargs = str(config[algorithm]['pl_trainer_kwargs'])
        if pl_trainer_kwargs == "cpu":
            drive = {"accelerator": pl_trainer_kwargs}
        else:
            drive = {"accelerator": pl_trainer_kwargs, "devices": [0]}
        return {"input_chunk_length": input_chunk_length, "model": model, "random_state": random_state, "n_epochs": n_epochs, "pl_trainer_kwargs": drive}
    elif algorithm == "forecaster-transformer":
        input_chunk_length = int(config[algorithm]['input_chunk_length'])
        output_chunk_length = int(config[algorithm]['output_chunk_length'])
        random_state = int(config[algorithm]['random_state'])
        n_epochs = int(config[algorithm]['n_epochs'])
        pl_trainer_kwargs = str(config[algorithm]['pl_trainer_kwargs'])
        if pl_trainer_kwargs == "cpu":
            drive = {"accelerator": pl_trainer_kwargs}
        else:
            drive = {"accelerator": pl_trainer_kwargs, "devices": [0]}
        return {"input_chunk_length": input_chunk_length, "output_chunk_length": output_chunk_length, "random_state": random_state, "n_epochs": n_epochs, "pl_trainer_kwargs": drive}

    elif algorithm == "forecaster-hw-add":
        sp = int(config[algorithm]['sp'])
        trend = str(config[algorithm]['trend'])
        seasonal = str(config[algorithm]['seasonal'])
        return {"sp": sp, "trend": trend, "seasonal": seasonal}
    elif algorithm == "forecaster-arima":
        sp = int(config[algorithm]['sp'])
        suppress_warnings = bool(config[algorithm]['suppress_warnings'])
        start_p = int(config[algorithm]['start_p'])
        start_q = int(config[algorithm]['start_q'])
        max_p = int(config[algorithm]['max_p'])
        max_q = int(config[algorithm]['max_q'])
        start_P = int(config[algorithm]['start_P'])
        seasonal = int(config[algorithm]['seasonal'])
        d = int(config[algorithm]['d'])
        D = int(config[algorithm]['D'])
        return {"sp": sp, "suppress_warnings": suppress_warnings, "start_p": start_p, "start_q": start_q,
                "max_p": max_p, "max_q": max_q, "start_P": start_P, "seasonal": seasonal, "d": d, "D": D}
    elif algorithm == "forecaster-sf-arima":
        sp = int(config[algorithm]['sp'])
        start_p = int(config[algorithm]['start_p'])
        start_q = int(config[algorithm]['start_q'])
        max_p = int(config[algorithm]['max_p'])
        max_q = int(config[algorithm]['max_q'])
        start_P = int(config[algorithm]['start_P'])
        seasonal = int(config[algorithm]['seasonal'])
        d = int(config[algorithm]['d'])
        D = int(config[algorithm]['D'])
        return {"sp": sp, "start_p": start_p, "start_q": start_q,
                "max_p": max_p, "max_q": max_q, "start_P": start_P, "seasonal": seasonal, "d": d, "D": D}
    elif algorithm == "forecaster-bats":
        sp = int(config[algorithm]['sp'])
        use_trend = bool(config[algorithm]['use_trend'])
        use_box_cox = bool(config[algorithm]['use_box_cox'])
        return {"sp": sp, "use_trend": use_trend, "use_box_cox": use_box_cox}
    elif algorithm == "forecaster-ets":
        sp = int(config[algorithm]['sp'])
        auto = bool(config[algorithm]['auto'])
        return {"sp": sp, "auto": auto}
    elif algorithm == "forecaster-croston":
        smoothing = float(config[algorithm]['smoothing'])
        return {"smoothing": smoothing}
    elif algorithm == "forecaster-unobs":
        level = bool(config[algorithm]['level'])
        trend = bool(config[algorithm]['trend'])
        sp = int(config[algorithm]['sp'])
        return {"level": level, "trend": trend, "seasonal": sp}
    elif algorithm == "forecaster-theta":
        sp = int(config[algorithm]['sp'])
        deseasonalize = bool(config[algorithm]['deseasonalize'])
        return {"sp": sp, "deseasonalize": deseasonalize}
    elif algorithm == "forecaster-rnn":
        input_size = int(config[algorithm]['input_size'])
        inference_input_size = int(config[algorithm]['inference_input_size'])
        return {"input_size": input_size, "inference_input_size": inference_input_size}
    elif algorithm == "colors":
        colors = config[algorithm]['plot']
        return colors
    elif algorithm == "other":
        return config

    # Your own default parameters #contributing
    #
    #elif algorithm == "your_algo_name":
    #    param_1 = int(config[algorithm]['param_1'])
    #    param_2 = config[algorithm]['param_2']
    #    param_3 = float(config[algorithm]['param_3'])
    #    return (param_1, param_2, param_3)

    else:
        print("(SYS) Default/Optimal config not found for this algorithm")
        return None


def config_impute_algorithm(incomp_data, algorithm, verbose=True):
    """
    Configure and execute algorithm for selected imputation imputer and pattern.

    Parameters
    ----------
    incomp_data : TimeSeries
        TimeSeries object containing dataset.

    algorithm : str
        Name of algorithm

    verbose : bool, optional
        Whether to display the contamination information (default is False).

    Returns
    -------
    BaseImputer
        Configured imputer instance with optimal parameters.
    """

    from imputegap.recovery.imputation import Imputation
    from imputegap.recovery.manager import TimeSeries

    alg_low = algorithm.lower()
    alg = alg_low.replace('_', '').replace('-', '')

    # 1st generation
    if alg == "cdrec":
        imputer = Imputation.MatrixCompletion.CDRec(incomp_data)
    elif alg == "stmvl":
        imputer = Imputation.PatternSearch.STMVL(incomp_data)
    elif alg == "iim":
        imputer = Imputation.MachineLearning.IIM(incomp_data)
    elif alg == "mrnn":
        imputer = Imputation.DeepLearning.MRNN(incomp_data)

    # 2nd generation
    elif alg == "iterativesvd" or alg == "itersvd":
        imputer = Imputation.MatrixCompletion.IterativeSVD(incomp_data)
    elif alg == "grouse":
        imputer = Imputation.MatrixCompletion.GROUSE(incomp_data)
    elif alg == "dynammo":
        imputer = Imputation.PatternSearch.DynaMMo(incomp_data)
    elif alg == "rosl":
        imputer = Imputation.MatrixCompletion.ROSL(incomp_data)
    elif alg == "softimpute" or alg == "softimp":
        imputer = Imputation.MatrixCompletion.SoftImpute(incomp_data)
    elif alg == "spirit":
        imputer = Imputation.MatrixCompletion.SPIRIT(incomp_data)
    elif alg == "svt":
        imputer = Imputation.MatrixCompletion.SVT(incomp_data)
    elif alg == "tkcm":
        imputer = Imputation.PatternSearch.TKCM(incomp_data)
    elif alg == "deepmvi":
        imputer = Imputation.DeepLearning.DeepMVI(incomp_data)
    elif alg == "brits":
        imputer = Imputation.DeepLearning.BRITS(incomp_data)
    elif alg == "mpin":
        imputer = Imputation.DeepLearning.MPIN(incomp_data)
    elif alg == "pristi":
        imputer = Imputation.DeepLearning.PRISTI(incomp_data)

    # 3rd generation
    elif alg == "knn" or alg == "knnimpute":
        imputer = Imputation.Statistics.KNNImpute(incomp_data)
    elif alg == "interpolation":
        imputer = Imputation.Statistics.Interpolation(incomp_data)
    elif alg == "meanseries" or alg == "meanimputebyseries":
        imputer = Imputation.Statistics.MeanImputeBySeries(incomp_data)
    elif alg == "minimpute":
        imputer = Imputation.Statistics.MinImpute(incomp_data)
    elif alg == "zeroimpute":
        imputer = Imputation.Statistics.ZeroImpute(incomp_data)
    elif alg == "trmf":
        imputer = Imputation.MatrixCompletion.TRMF(incomp_data)
    elif alg == "mice":
        imputer = Imputation.MachineLearning.MICE(incomp_data)
    elif alg == "missforest":
        imputer = Imputation.MachineLearning.MissForest(incomp_data)
    elif alg == "xgboost":
        imputer = Imputation.MachineLearning.XGBOOST(incomp_data)
    elif alg == "missnet":
        imputer = Imputation.DeepLearning.MissNet(incomp_data)
    elif alg == "gain":
        imputer = Imputation.DeepLearning.GAIN(incomp_data)
    elif alg == "grin":
        imputer = Imputation.DeepLearning.GRIN(incomp_data)
    elif alg == "bayotide":
        imputer = Imputation.DeepLearning.BayOTIDE(incomp_data)
    elif alg == "hkmft":
        imputer = Imputation.DeepLearning.HKMF_T(incomp_data)
    elif alg == "bitgraph":
        imputer = Imputation.DeepLearning.BitGraph(incomp_data)
    elif alg == "meanimpute":
        imputer = Imputation.Statistics.MeanImpute(incomp_data)

    # 4th generation
    elif alg == "nuwats":
        imputer = Imputation.LLMs.NuwaTS(incomp_data)
    elif alg == "gpt4ts":
        imputer = Imputation.LLMs.GPT4TS(incomp_data)

    # your own implementation #contributing
    #
    #elif alg == "your_algo_name":
    #    imputer = Imputation.MyFamily.NewAlg(incomp_data)

    else:
        raise ValueError(f"(IMP) Algorithm '{algorithm}' not recognized, please choose your algorithm from this list:\n\t{TimeSeries().algorithms}")
        imputer = None

    if imputer is not None:
        imputer.verbose = verbose

    return imputer


def save_optimization(optimal_params, algorithm="cdrec", dataset="", optimizer="b", file_name=None):
    """
    Save the optimization parameters to a TOML file for later use without recomputing.

    Parameters
    ----------
    optimal_params : dict
        Dictionary of the optimal parameters.
    algorithm : str, optional
        The name of the imputation algorithm (default is 'cdrec').
    dataset : str, optional
        The name of the dataset (default is an empty string).
    optimizer : str, optional
        The name of the optimizer used (default is 'b').
    file_name : str, optional
        The name of the TOML file to save the results (default is None).

    Returns
    -------
    None
    """
    if file_name is None:
        file_name = "./imputegap_assets/params/optimal_parameters_" + str(optimizer) + "_" + str(dataset) + "_" + str(algorithm) + ".toml"
    else:
        file_name += "optimal_parameters_" + str(optimizer) + "_" + str(dataset) + "_" + str(algorithm) + ".toml"

    dir_name = os.path.dirname(file_name)
    if dir_name and not os.path.exists(dir_name):
        os.makedirs(dir_name)

    if algorithm == "cdrec" or algorithm == "CDRec":
        params_to_save = {
            "rank": int(optimal_params[0]),
            "eps": optimal_params[1],
            "iters": int(optimal_params[2])
    }
    elif algorithm == "mrnn" or algorithm == "MRNN":
        params_to_save = { "hidden_dim": int(optimal_params[0]),
            "learning_rate": optimal_params[1],
            "num_iter": int(optimal_params[2]),
            "seq_len": 7  # Default value
        }
    elif algorithm == "stmvl" or algorithm == "STMVL" or algorithm == "ST-MVL":
        params_to_save = {
            "window_size": int(optimal_params[0]),
            "gamma": optimal_params[1],
            "alpha": int(optimal_params[2])
        }
    elif algorithm == "iim" or algorithm == "IIM":
        params_to_save = {
            "learning_neighbors": int(optimal_params[0])
        }

    elif algorithm == "iterative_svd" or algorithm == "IterativeSVD":
        params_to_save = {
            "rank": int(optimal_params[0])
        }
    elif algorithm == "grouse" or algorithm == "GROUSE":
        params_to_save= {
            "max_rank": int(optimal_params[0])
        }
    elif algorithm == "rosl" or algorithm == "ROSL":
        params_to_save = {
            "rank": int(optimal_params[0]),
            "regularization": optimal_params[1]
        }
    elif algorithm == "soft_impute" or algorithm == "SoftImpute":
        params_to_save = {
            "max_rank": int(optimal_params[0])
        }
    elif algorithm == "spirit" or algorithm == "SPIRIT":
        params_to_save = {
            "k": int(optimal_params[0]),
            "w": int(optimal_params[1]),
            "lvalue": optimal_params[2]
        }
    elif algorithm == "svt" or algorithm == "SVT":
        params_to_save = {
            "tau": optimal_params[0],
            "delta": optimal_params[1],
            "max_iter": int(optimal_params[2])
        }
    elif algorithm == "dynammo" or algorithm == "DynamoMO":
        params_to_save = {
            "h": int(optimal_params[0]),
            "max_iteration": int(optimal_params[1]),
            "approximation": bool(optimal_params[2])
        }
    elif algorithm == "tkcm":
        params_to_save = {
            "rank": int(optimal_params[0])
        }
    elif algorithm == "brits":
        params_to_save = {
            "model": optimal_params[0],
            "epoch": int(optimal_params[1]),
            "batch_size": int(optimal_params[2]),
            "hidden_layers": int(optimal_params[3]),
            "num_workers": int(optimal_params[4])
        }
    elif algorithm == "deep_mvi":
        params_to_save = {
            "max_epoch": int(optimal_params[0]),
            "patience": int(optimal_params[1]),
            "lr": float(optimal_params[2])
        }
    elif algorithm == "mpin":
        params_to_save = {
            "incre_mode": optimal_params[0],
            "window": int(optimal_params[1]),
            "k": int(optimal_params[2]),
            "learning_rate": optimal_params[3],
            "weight_decay": optimal_params[4],
            "epochs": int(optimal_params[5]),
            "num_of_iteration": int(optimal_params[6]),
            "threshold": optimal_params[7],
            "base": optimal_params[8]
        }
    elif algorithm == "pristi":
        params_to_save = {
            "target_strategy": optimal_params[0],
            "unconditional": bool(optimal_params[1]),
            "batch_size": bool(optimal_params[2]),
            "embedding": bool(optimal_params[3]),
            "num_workers": bool(optimal_params[4]),
            "seed": 42,  # Default seed
        }
    elif algorithm == "knn" or algorithm == "knn_impute":
        params_to_save = {
            "k": int(optimal_params[0]),
            "weights": str(optimal_params[1])
        }
    elif algorithm == "interpolation":
        params_to_save = {
            "method": str(optimal_params[0]),
            "poly_order": int(optimal_params[1])
        }
    elif algorithm == "mice":
        params_to_save = {
            "max_iter": int(optimal_params[0]),
            "tol": float(optimal_params[1]),
            "initial_strategy": str(optimal_params[2]),
            "seed": 42
        }
    elif algorithm == "miss_forest":
        params_to_save = {
            "n_estimators": int(optimal_params[0]),
            "max_iter": int(optimal_params[1]),
            "max_features": str(optimal_params[2]),
            "seed": 42
        }
    elif algorithm == "xgboost":
        params_to_save = {
            "n_estimators": int(optimal_params[0]),
            "seed": 42
        }
    elif algorithm == "miss_net":
        params_to_save = {
            "alpha": float(optimal_params[0]),
            "beta": float(optimal_params[1]),
            "L": int(optimal_params[2]),
            "n_cl": int(optimal_params[3]),
            "max_iter": int(optimal_params[4]),
            "tol": float(optimal_params[5]),
            "random_init": bool(optimal_params[6])
        }
    elif algorithm == "gain":
        params_to_save = {
            "batch_size": int(optimal_params[0]),
            "hint_rate": float(optimal_params[1]),
            "alpha": int(optimal_params[2]),
            "epoch": int(optimal_params[3])
        }
    elif algorithm == "grin":
        params_to_save = {
            "d_hidden": int(optimal_params[0]),
            "lr": float(optimal_params[1]),
            "batch_size": int(optimal_params[2]),
            "window": int(optimal_params[3]),
            "alpha": int(optimal_params[4]),
            "patience": int(optimal_params[5]),
            "epochs": int(optimal_params[6]),
            "workers": int(optimal_params[7])
        }
    elif algorithm == "bay_otide":
        params_to_save = {
            "K_trend": int(optimal_params[0]),
            "K_season": int(optimal_params[1]),
            "n_season": int(optimal_params[2]),
            "K_bias": int(optimal_params[3]),
            "time_scale": int(optimal_params[4]),
            "a0": float(optimal_params[5]),
            "b0": float(optimal_params[6]),
            "v": float(optimal_params[7]),
            "tr_ratio": float(optimal_params[8])
        }
    elif algorithm == "hkmf_t":
        params_to_save = {
            "tags": optimal_params[0],
            "data_names": optimal_params[1],
            "epoch": int(optimal_params[2]),
        }
    elif algorithm == "bit_graph":
        params_to_save = {
            "node_number": int(optimal_params[0]),
            "kernel_set": optimal_params[1],
            "dropout": float(optimal_params[2]),
            "subgraph_size": int(optimal_params[3]),
            "node_dim": int(optimal_params[4]),
            "seq_len": int(optimal_params[5]),
            "lr": float(optimal_params[6]),
            "batch_size": float(optimal_params[7]),
            "epoch": int(optimal_params[8]),
            "num_workers": int(optimal_params[9]),
            "seed": int(optimal_params[10]),
        }
    elif algorithm == "nuwats" or algorithm == "NUWATS":
        params_to_save = {
            "seq_length": int(optimal_params[0]),
            "patch_size": optimal_params[1],
            "batch_size": float(optimal_params[2]),
            "pred_length": int(optimal_params[3]),
            "label_length": int(optimal_params[4]),
            "enc_in": int(optimal_params[5]),
            "dec_in": float(optimal_params[6]),
            "c_out": float(optimal_params[7]),
            "gpt_layers": int(optimal_params[8]),
            "num_workers": int(optimal_params[9]),
            "seed": int(optimal_params[10]),
        }
    elif algorithm == "gpt4ts" or algorithm == "GPT4TS":
        params_to_save = {
            "seq_length": int(optimal_params[0]),
            "patch_size": optimal_params[1],
            "batch_size": float(optimal_params[2]),
            "pred_length": int(optimal_params[3]),
            "label_length": int(optimal_params[4]),
            "enc_in": int(optimal_params[5]),
            "dec_in": float(optimal_params[6]),
            "c_out": float(optimal_params[7]),
            "gpt_layers": int(optimal_params[8]),
            "num_workers": int(optimal_params[9]),
            "seed": int(optimal_params[10]),
        }


    # Your own optimal save parameters #contributing
    #
    #elif algorithm == "your_algo_name":
    #    params_to_save = {
    #        "param_1": int(optimal_params[0]),
    #    "param_2": optimal_params[1],
    #    "param_3": float(optimal_params[2]),
    #}


    else:
        print(f"\n\t\t(SYS) Algorithm {algorithm} is not recognized.")
        return

    try:
        with open(file_name, 'w') as file:
            toml.dump(params_to_save, file)
        print(f"\n(SYS) Optimization parameters successfully saved to {file_name}")
    except Exception as e:
        print(f"\n(SYS) An error occurred while saving the file: {e}")


def check_family(family, algorithm):
    # Normalize input
    norm_input = algorithm.lower().replace("_", "").replace("-", "")

    for full_name in list_of_algorithms_with_families():
        if full_name.startswith("DeepLearning."):
            suffix = full_name.split(".", 1)[1]
            norm_suffix = suffix.lower().replace("_", "").replace("-", "")

            if norm_input == norm_suffix:
                return True
    return False


def config_contamination(ts, pattern, dataset_rate=0.4, series_rate=0.4, block_size=10, offset=0.1, seed=True, limit=1, shift=0.05, std_dev=0.5, explainer=False, probabilities=None, verbose=True):
    """
    Configure and execute contamination for selected imputation algorithm and pattern.

    Parameters
    ----------
    rate : float
        Mean parameter for contamination missing percentage rate.
    ts_test : TimeSeries
        A TimeSeries object containing dataset.
    pattern : str
        Type of contamination pattern (e.g., "mcar", "mp", "blackout", "disjoint", "overlap", "gaussian").
    block_size_mcar : int
        Size of blocks removed in MCAR

    Returns
    -------
    TimeSeries
        TimeSeries object containing contaminated data.
    """

    from imputegap.recovery.manager import TimeSeries

    pattern_low = pattern.lower()
    ptn = pattern_low.replace('_', '').replace('-', '')

    if ptn == "mcar" or ptn == "missing_completely_at_random":
        incomp_data = ts.Contamination.mcar(input_data=ts.data, rate_dataset=dataset_rate, rate_series=series_rate, block_size=block_size, offset=offset, seed=seed, explainer=explainer, verbose=verbose)
    elif ptn == "mp" or ptn == "missingpercentage" or ptn == "aligned":
        incomp_data = ts.Contamination.aligned(input_data=ts.data, rate_dataset=dataset_rate, rate_series=series_rate, offset=offset, explainer=explainer, verbose=verbose)
    elif ptn == "ps" or ptn == "percentageshift" or ptn == "scattered" or ptn == "scatter":
        incomp_data = ts.Contamination.scattered(input_data=ts.data, rate_dataset=dataset_rate, rate_series=series_rate, offset=offset, seed=seed, explainer=explainer, verbose=verbose)
    elif ptn == "disjoint":
        incomp_data = ts.Contamination.disjoint(input_data=ts.data, rate_series=dataset_rate, limit=1, offset=offset, verbose=verbose)
    elif ptn == "overlap":
        incomp_data = ts.Contamination.overlap(input_data=ts.data, rate_series=dataset_rate, limit=limit, shift=shift, offset=offset, verbose=verbose)
    elif ptn == "gaussian":
        incomp_data = ts.Contamination.gaussian(input_data=ts.data, rate_dataset=dataset_rate, rate_series=series_rate, std_dev=std_dev, offset=offset, seed=seed, explainer=explainer, verbose=verbose)
    elif ptn == "distribution" or pattern == "dist":
        incomp_data = ts.Contamination.distribution(input_data=ts.data, rate_dataset=dataset_rate, rate_series=series_rate, probabilities_list=probabilities, offset=offset, seed=seed, explainer=explainer, verbose=verbose)
    elif ptn == "blackout":
        incomp_data = ts.Contamination.blackout(input_data=ts.data, series_rate=dataset_rate, offset=offset, verbose=verbose)
    else:
        raise ValueError(f"\n(CONT) Pattern '{pattern}' not recognized, please choose your algorithm on this list :\n\t{TimeSeries().patterns}\n")
        incomp_data = None

    return incomp_data


def config_forecaster(model, params):
        """
        Configure and execute forecaster model for downstream analytics

        Parameters
        ----------
        model : str
            name of the forcaster model
        params : list of params
            List of paramaters for a forcaster model

        Returns
        -------
        Forecaster object (SKTIME/DART)
            Forecaster object for downstream analytics
        """

        from imputegap.recovery.manager import TimeSeries

        model_low = model.lower()
        mdl = model_low.replace('_', '').replace('-', '')

        if mdl == "prophet":
            from sktime.forecasting.fbprophet import Prophet
            forecaster = Prophet(**params)
        elif mdl == "expsmoothing":
            from sktime.forecasting.exp_smoothing import ExponentialSmoothing
            forecaster = ExponentialSmoothing(**params)
        elif mdl == "nbeats":
            from darts.models import NBEATSModel
            forecaster = NBEATSModel(**params)
        elif mdl == "xgboost":
            from darts.models.forecasting.xgboost import XGBModel
            forecaster = XGBModel(**params)
        elif mdl == "lightgbm":
            from darts.models.forecasting.lgbm import LightGBMModel
            forecaster = LightGBMModel(**params)
        elif mdl == "lstm":
            from darts.models.forecasting.rnn_model import RNNModel
            forecaster = RNNModel(**params)
        elif mdl == "deepar":
            from darts.models.forecasting.rnn_model import RNNModel
            forecaster = RNNModel(**params)
        elif mdl == "transformer":
            from darts.models.forecasting.transformer_model import TransformerModel
            forecaster = TransformerModel(**params)
        elif mdl == "hwadd":
            from sktime.forecasting.exp_smoothing import ExponentialSmoothing
            forecaster = ExponentialSmoothing(**params)
        elif mdl == "arima":
            from sktime.forecasting.arima import AutoARIMA
            forecaster = AutoARIMA(**params)
        elif mdl == "sf-arima":
            from sktime.forecasting.statsforecast import StatsForecastAutoARIMA
            forecaster = StatsForecastAutoARIMA(**params)
            forecaster.set_config(warnings='off')
        elif mdl == "bats":
            from sktime.forecasting.bats import BATS
            forecaster = BATS(**params)
        elif mdl == "ets":
            from sktime.forecasting.ets import AutoETS
            forecaster = AutoETS(**params)
        elif mdl == "croston":
            from sktime.forecasting.croston import Croston
            forecaster = Croston(**params)
        elif mdl == "theta":
            from sktime.forecasting.theta import ThetaForecaster
            forecaster = ThetaForecaster(**params)
        elif mdl == "unobs":
            from sktime.forecasting.structural import UnobservedComponents
            forecaster = UnobservedComponents(**params)
        elif mdl == "naive":
            from sktime.forecasting.naive import NaiveForecaster
            forecaster = NaiveForecaster(**params)
        else:
            raise ValueError(f"\n(DOWN) Forecasting model '{model}' not recognized, please choose your algorithm on this list :\n\t{TimeSeries().forecasting_models}\n")
            forecaster = None

        return forecaster



def __marshal_as_numpy_column(__ctype_container, __py_sizen, __py_sizem):
    """
    Marshal a ctypes container as a numpy column-major array.

    Parameters
    ----------
    __ctype_container : ctypes.Array
        The input ctypes container (flattened matrix).
    __py_sizen : int
        The number of rows in the numpy array.
    __py_sizem : int
        The number of columns in the numpy array.

    Returns
    -------
    numpy.ndarray
        A numpy array reshaped to the original matrix dimensions (row-major order).
    """
    __numpy_marshal = __numpy_import.array(__ctype_container).reshape(__py_sizem, __py_sizen).T;

    return __numpy_marshal;


def __marshal_as_native_column(__py_matrix):
    """
    Marshal a numpy array as a ctypes flat container for passing to native code.

    Parameters
    ----------
    __py_matrix : numpy.ndarray
        The input numpy matrix (2D array).

    Returns
    -------
    ctypes.Array
        A ctypes array containing the flattened matrix (in column-major order).
    """
    __py_input_flat = __numpy_import.ndarray.flatten(__py_matrix.T);
    __ctype_marshal = __numpy_import.ctypeslib.as_ctypes(__py_input_flat);

    return __ctype_marshal;


def display_title(title="Master Thesis", aut="Quentin Nater", lib="ImputeGAP", university="University Fribourg"):
    """
    Display the title and author information.

    Parameters
    ----------
    title : str, optional
        The title of the thesis (default is "Master Thesis").
    aut : str, optional
        The author's name (default is "Quentin Nater").
    lib : str, optional
        The library or project name (default is "ImputeGAP").
    university : str, optional
        The university or institution (default is "University Fribourg").

    Returns
    -------
    None
    """

    print("=" * 100)
    print(f"{title} : {aut}")
    print("=" * 100)
    print(f"    {lib} - {university}")
    print("=" * 100)


def search_path(set_name="test"):
    """
    Find the accurate path for loading test files.

    Parameters
    ----------
    set_name : str, optional
        Name of the dataset (default is "test").

    Returns
    -------
    str
        The correct file path for the dataset.
    """

    if set_name in list_of_datasets():
        return set_name + ".txt"
    else:
        filepath = "../imputegap/datasets/" + set_name

        if not os.path.exists(filepath):
            filepath = filepath[1:]
        return filepath


def get_missing_ratio(incomp_data):
    """
    Check whether the proportion of missing values in the contaminated data is acceptable
    for training a deep learning model.

    Parameters
    ----------
    incomp_data : TimeSeries (numpy array)
            TimeSeries object containing dataset.

    Returns
    -------
    bool
        True if the missing data ratio is less than or equal to 40%, False otherwise.
    """
    import numpy as np

    miss_m = incomp_data
    total_values = miss_m.size
    missing_values = np.isnan(miss_m).sum()
    missing_ratio = missing_values / total_values

    return missing_ratio


def verification_limitation(percentage, low_limit=0.01, high_limit=1.0):
    """
    Format and verify that the percentage given by the user is within acceptable bounds.

    Parameters
    ----------
    percentage : float
        The percentage value to be checked and potentially adjusted.
    low_limit : float, optional
        The lower limit of the acceptable percentage range (default is 0.01).
    high_limit : float, optional
        The upper limit of the acceptable percentage range (default is 1.0).

    Returns
    -------
    float
        Adjusted percentage based on the limits.

    Raises
    ------
    ValueError
        If the percentage is outside the accepted limits.

    Notes
    -----
    - If the percentage is between 1 and 100, it will be divided by 100 to convert it to a decimal format.
    - If the percentage is outside the low and high limits, the function will print a warning and return the original value.
    """
    if low_limit <= percentage <= high_limit:
        return percentage  # No modification needed

    elif 1 <= percentage <= 100:
        print(f"The percentage {percentage} is between 1 and 100. Dividing by 100 to convert to a decimal.")
        return percentage / 100

    else:
        raise ValueError(f"The percentage {percentage} is out of the acceptable range.")


def dl_integration_transformation(input_matrix, tr_ratio=0.8, inside_tr_cont_ratio=0.2, split_ts=1, split_val=0, nan_val=-99999, prevent_leak=True, offset=0.05, block_selection=True, seed=42, verbose=False):
    """
        Prepares contaminated data and corresponding masks for deep learning-based imputation training,
        validation, and testing.

        This function simulates missingness in a controlled way, optionally prevents information leakage,
        and produces masks for training, testing, and validation using different contamination strategies.

        Parameters:
        ----------
        input_matrix : np.ndarray
            The complete input time series data matrix of shape [T, N] (time steps × variables).

        tr_ratio : float, default=0.8
            The fraction of data to reserve for training when constructing the test contamination mask.

        inside_tr_cont_ratio : float, default=0.2
            The proportion of values to randomly drop inside the training data for internal contamination.

        split_ts : float, default=1
            Proportion of the total contaminated data assigned to the test set.

        split_val : float, default=0
            Proportion of the total contaminated data assigned to the validation set.

        nan_val : float, default=-99999
            Value used to represent missing entries in the masked matrix.
            nan_val=-1 can be used to set mean values

        prevent_leak : bool, default=True
            Replace the value of NaN with a high number to prevent leakage.

        offset : float, default=0.05
            Minimum temporal offset in the begining of the series

        block_selection : bool, default=True
            Whether to simulate missing values in contiguous blocks (True) or randomly (False).

        seed : int, default=42
            Seed for NumPy random number generation to ensure reproducibility.

        verbose : bool, default=False
            Whether to print logging/debug information during execution.

        Returns:
        -------
        cont_data_matrix : np.ndarray
            The input matrix with synthetic missing values introduced.

        mask_train : np.ndarray
            Boolean mask of shape [T, N] indicating the training contamination locations (True = observed, False = missing).

        mask_test : np.ndarray
            Boolean mask of shape [T, N] indicating the test contamination locations.

        mask_valid : np.ndarray
            Boolean mask of shape [T, N] indicating the validation contamination locations.

        error : bool
            Tag which is triggered if the operation is impossible.
    """

    cont_data_matrix = input_matrix.copy()
    original_missing_ratio = get_missing_ratio(cont_data_matrix)

    cont_data_matrix, new_mask, error = prepare_testing_set(incomp_m=cont_data_matrix, original_missing_ratio=original_missing_ratio, block_selection=block_selection, tr_ratio=tr_ratio, verbose=verbose)

    if prevent_leak:
        if nan_val == -1:
            import numpy as np
            nan_val = np.nanmean(input_matrix)
            print(f"\nNaN replacement Mean Value : {nan_val}\n")
        cont_data_matrix = prevent_leakage(cont_data_matrix, new_mask, nan_val, verbose)

    mask_test, mask_valid, nbr_nans = split_mask_bwt_test_valid(cont_data_matrix, test_rate=split_ts, valid_rate=split_val, nan_val=nan_val, verbose=verbose, seed=seed)
    mask_train = generate_random_mask(gt=cont_data_matrix, mask_test=mask_test, mask_valid=mask_valid, droprate=inside_tr_cont_ratio, offset=offset, verbose=verbose, seed=seed)

    return cont_data_matrix, mask_train, mask_test, mask_valid, error


def prepare_fixed_testing_set(incomp_m, tr_ratio=0.8, offset=0.05, block_selection=True, verbose=True):
    """
    Introduces additional missing values (NaNs) into a data matrix to match a specified training ratio.

    This function modifies a copy of the input matrix `incomp_m` by introducing NaNs
    such that the proportion of observed (non-NaN) values matches the desired `tr_ratio`.
    It returns the modified matrix and the corresponding missing data mask.

    Parameters
    ----------
    incomp_m : np.ndarray
       A 2D NumPy array with potential pre-existing NaNs representing missing values.

    tr_ratio : float
       Desired ratio of observed (non-NaN) values in the output matrix. Must be in the range (0, 1).

    offset : float
        Protected zone in the begining of the series

    block_selection : bool
        Select the missing values by blocks or randomly (True, is by block)

    verbose : bool
        Whether to print debug info.

    Returns
    -------
    data_matrix_cont : np.ndarray
       The modified matrix with additional NaNs introduced to match the specified training ratio.

    new_mask : np.ndarray
       A boolean mask of the same shape as `data_matrix_cont` where True indicates missing (NaN) entries.

    Raises
    ------
    AssertionError:
       If the final observed and missing ratios deviate from the target by more than 1%.

    Notes
    -----
        - The function assumes that the input contains some non-NaN entries.
        - NaNs are added in row-major order from the list of available (non-NaN) positions.
    """

    import numpy as np

    data_matrix_cont = incomp_m.copy()

    target_ratio = 1 - tr_ratio
    total_values = data_matrix_cont.size
    target_n_nan = int(target_ratio * total_values)

    # 2) Current number of NaNs
    current_n_nan = np.isnan(data_matrix_cont).sum()
    n_new_nans = target_n_nan - current_n_nan

    available_mask = ~np.isnan(data_matrix_cont)

    offset_vals = int(offset * data_matrix_cont.shape[1])
    for row in range(data_matrix_cont.shape[0]):
        available_mask[row, :offset_vals] = False  # protect leftmost `offset` columns in each row

    available_indices = np.argwhere(available_mask)

    # 3) Pick indices to contaminate
    if n_new_nans > 0:
        if block_selection :
            chosen_indices = available_indices[:n_new_nans]
        else:
            np.random.seed(42)
            chosen_indices = available_indices[np.random.choice(len(available_indices), n_new_nans, replace=False)]

        for i, j in chosen_indices:
            data_matrix_cont[i, j] = np.nan

    # 4) check ratio
    n_total = data_matrix_cont.size
    n_nan = np.isnan(data_matrix_cont).sum()
    n_not_nan = n_total - n_nan

    # Compute actual ratios
    missing_ratio = n_nan / n_total
    observed_ratio = n_not_nan / n_total

    # Check if they match expectations (within a small tolerance)
    assert abs(missing_ratio - target_ratio) < 0.01, f"Missing ratio {missing_ratio} is not {target_ratio}"
    assert abs(observed_ratio - tr_ratio) < 0.01, f"Missing ratio {observed_ratio} is not {tr_ratio}"

    # Create the new mask
    new_mask = np.isnan(data_matrix_cont)
    new_m = data_matrix_cont.copy()

    if verbose:
        print(f"(DL): TEST-SET > Test set fixed to {int(round(target_ratio*100))}% of the dataset, for {target_n_nan} values, add test values: {n_new_nans}")

    return new_m, new_mask

def split_mask_bwt_test_valid(data_matrix, test_rate=0.8, valid_rate=0.2, nan_val=None, verbose=False, seed=42):
    """
    Dispatch NaN positions in data_matrix to test and validation masks only.

    Parameters
    ----------
    data_matrix : numpy.ndarray
        Input matrix containing NaNs to be split.

    test_rate : float
        Proportion of NaNs to assign to the test set (default is 0.8).

    valid_rate : float
        Proportion of NaNs to assign to the validation set (default is 0.2).
        test_rate + valid_rate must equal 1.0.

    verbose : bool
        Whether to print debug info.

    seed : int, optional
        Random seed for reproducibility.

    Returns
    -------
    tuple
        test_mask : numpy.ndarray
            Binary mask indicating positions of NaNs in the test set.

        valid_mask : numpy.ndarray
            Binary mask indicating positions of NaNs in the validation set.

        n_nan : int
            Total number of NaN values found in the input matrix.
    """
    import numpy as np

    assert np.isclose(test_rate + valid_rate, 1.0), "test_rate and valid_rate must sum to 1.0"

    if seed is not None:
        np.random.seed(seed)

    if nan_val is None:
        nan_mask = np.isnan(data_matrix)
    else:
        nan_mask = data_matrix == nan_val

    nan_indices = np.argwhere(nan_mask)
    np.random.shuffle(nan_indices)

    n_nan = len(nan_indices)
    n_test = int(n_nan * test_rate)
    n_valid = n_nan - n_test

    if verbose:
        print(f"\n(DL): MASKS > creating mask (testing, validation): Total NaNs = {n_nan}")
        print(f"(DL): TEST-MASK > creating mask: Assigned to test = {n_test}")
        print(f"(DL): VALID-MASK > creating mask: Assigned to valid = {n_valid}")

    test_idx = nan_indices[:n_test]
    valid_idx = nan_indices[n_test:]

    mask_test = np.zeros_like(data_matrix, dtype=np.uint8)
    mask_valid = np.zeros_like(data_matrix, dtype=np.uint8)

    mask_test[tuple(test_idx.T)] = 1
    mask_valid[tuple(valid_idx.T)] = 1

    if verbose:
        print(f"(DL): TEST-MASK > Test mask NaNs: {mask_test.sum()}")
        print(f"(DL): VALID-MASK > Valid mask NaNs: {mask_valid.sum()}\n")

    return mask_test, mask_valid, n_nan


def generate_random_mask(gt, mask_test, mask_valid, droprate=0.2, offset=None, verbose=False, seed=42):
    """
    Generate a random training mask over the non-NaN entries of gt, excluding positions
    already present in the test and validation masks.

    Parameters
    ----------
    gt : numpy.ndarray
        Ground truth data (no NaNs).
    mask_test : numpy.ndarray
        Binary mask indicating test positions.
    mask_valid : numpy.ndarray
        Binary mask indicating validation positions.
    droprate : float
        Proportion of eligible entries to include in the training mask.
    offset : float
        Protect of not the offset of the dataset
    verbose : bool
        Whether to print debug info.
    seed : int, optional
        Random seed for reproducibility.

    Returns
    -------
    numpy.ndarray
        Binary mask indicating training positions.
    """
    import numpy as np

    assert gt.shape == mask_test.shape == mask_valid.shape, "All input matrices must have the same shape"

    if seed is not None:
        np.random.seed(seed)


    # Valid positions: non-NaN and not in test/valid masks
    num_offset = 0
    mask_offset = np.zeros_like(gt, dtype=np.uint8)
    if offset is not None:
        if offset > droprate:
            offset = droprate
        mask_offset[:, :int(offset * gt.shape[1])] = 1
        num_offset = np.sum(mask_offset)


    occupied_mask = (mask_test + mask_valid + mask_offset).astype(bool)
    eligible_mask = (~np.isnan(gt)) & (~occupied_mask)

    eligible_indices = np.argwhere(eligible_mask)

    n_train = int(len(eligible_indices) * droprate) + int(num_offset*droprate)

    np.random.shuffle(eligible_indices)
    selected_indices = eligible_indices[:n_train]

    mask_train = np.zeros_like(gt, dtype=np.uint8)
    mask_train[tuple(selected_indices.T)] = 1

    if verbose:
        print(f"(DL): TRAIN-MASK > eligible entries: {len(eligible_indices)}")
        print(f"(DL): TRAIN-MASK > selected training entries: {n_train}\n")

    # Sanity check: no overlap between training and test masks
    overlap = np.logical_and(mask_train, mask_test).sum()
    assert overlap == 0, f"Overlap detected between training and test masks: {overlap} entries."

    # Sanity check: no overlap between training and test masks
    overlap = np.logical_and(mask_train, mask_valid).sum()
    assert overlap == 0, f"Overlap detected between training and test masks: {overlap} entries."

    if verbose:
        print(f"(DL): TRAIN-MASK > Train mask NaNs: {mask_train.sum()}\n")

    return mask_train

def prevent_leakage(matrix, mask, replacement=0, verbose=True):
    """
        Replaces missing values in a matrix to prevent data leakage during evaluation.

        This function replaces all entries in `matrix` that are marked as missing in `mask`
        with a specified `replacement` value (default is 0). It then checks to ensure that
        there are no remaining NaNs in the matrix and that at least one replacement occurred.

        Parameters
        ----------
        matrix : np.ndarray
            A NumPy array potentially containing missing values (NaNs).

        mask : np.ndarray
            A boolean mask of the same shape as `matrix`, where True indicates positions
            to be replaced (typically where original values were NaN).

        replacement : float or int, optional
            The value to use in place of missing entries. Defaults to 0.

        verbose : bool
            Whether to print debug info.

        Returns
        -------
        matrix : np.ndarray
            The matrix with missing entries replaced by the specified value.

        Raises
        ------
        AssertionError:
            If any NaNs remain in the matrix after replacement, or if no replacements were made.

        Notes
        -----
            - This function is typically used before evaluation to ensure the model does not
              access ground truth values where data was originally missing.
    """

    import numpy as np

    matrix[mask] = replacement

    assert not np.isnan(matrix).any(), "matrix still contains NaNs"
    assert (matrix == replacement).any(), "matrix does not contain any zeros"

    if verbose:
        print(f"\n(DL) Reset all testing matrix values to {replacement} to prevent data leakage.")

    return matrix

def prepare_testing_set(incomp_m, original_missing_ratio, block_selection=True, tr_ratio=0.8, verbose=True):
    import numpy as np

    error = False
    mask_original_nan = np.isnan(incomp_m)

    if verbose:
        print(f"\n(DL) TEST-SET : testing ratio to reach = {1-tr_ratio:.2%}")
        print(f"\n(DL) TEST-SET : original missing ratio = {original_missing_ratio:.2%}")
        print(f"(DL) TEST-SET : original missing numbers = {np.sum(mask_original_nan)}")

    if original_missing_ratio > 1-tr_ratio:
        print(f"\n(ERROR) The proportion of original missing values is too high and will corrupt the training set.\n\tPlease consider reducing the percentage contamination pattern [{original_missing_ratio:.2%}] or decreasing the training ratio [{tr_ratio:.2%}].\n")
        return incomp_m, mask_original_nan, True

    if abs((1-tr_ratio) - original_missing_ratio) > 0.01:
        new_m, new_mask = prepare_fixed_testing_set(incomp_m, tr_ratio, block_selection=block_selection, verbose=verbose)

        if verbose:
            print(f"(DL) TEST-SET : building of the test set to reach a fix ratio of {1 - tr_ratio:.2%}...")
            final_ratio = get_missing_ratio(new_m)
            print(f"(DL) TEST-SET : final artificially missing ratio for test set = {final_ratio:.2%}")
            print(f"(DL) TEST-SET : final number of rows with NaN values = {np.sum(np.isnan(new_m).any(axis=1))}")
            print(f"(DL) TEST-SET : final artificially missing numbers = {np.sum(new_mask)}\n")

    else:
        new_m = incomp_m
        new_mask = mask_original_nan.copy()

    return new_m, new_mask, error



def compute_rank_check(M, rank, verbose=True):
    """
    Validates and adjusts the rank used in matrix operations based on the number of time series.

    Parameters
    ----------
    M : int
        Number of series
    rank : int
        The desired rank (e.g., for matrix factorization or low-rank approximation).
    verbose : bool
        Print the error or not

    Returns
    -------
    rank: int
        A valid rank value, adjusted to avoid exceeding the number of available series.
    """
    if rank >= M-1:
        if verbose:
            print(f"ERROR: Rank choosen to high for the number of series: {rank} >= {M}.\n\tRank reduced to 2.")
        return 2
    else:
        return rank


def compute_seq_length(M):
    """
    Compute a sequence length based on the input length `M` using heuristic rules.

    Parameters
    ----------
    M : int
        Number of series in the dataset.

    Returns
    -------
    seq_length: int
        A derived sequence length appropriate for processing or windowing.
    """

    seq_length = 1
    if M > 5000:
        seq_length = 3000
    elif M > 3000:
        seq_length = 1400
    elif M > 2000:
        seq_length = 1000
    elif M > 1000:
        seq_length = 600
    elif M > 300:
        seq_length = 100
    elif M > 30:
        seq_length = 16
    else:
        if M % 5 == 0:
            seq_length = M // 5
        elif M % 6 == 0:
            seq_length = M // 6
        elif M % 2 == 0:
            seq_length = M // 2 - 2
            if seq_length < 1:
                seq_length = 1
        elif M % 3 == 0:
            seq_length = M // 3

    return seq_length


def compute_batch_size(data, min_size=4, max_size=16, divisor=2, verbose=True):
    """
    Compute an appropriate batch size based on the input data shape.

    The batch size is computed as `min(M // 2, max_size)`, where M is the number of samples.
    If this computed batch size is less than `min_size`, it is set to `min_size` instead.

    Parameters
    ----------
    data : np.ndarray or torch.Tensor
        Input 2D data of shape (M, N), where M is the number of samples.
    min_size : int, optional
        Minimum allowed batch size. Default is 4.
    max_size : int, optional
        Maximum allowed batch size. Default is 16.
    divisor : int, optional
        Divisor on the shape of the dataset. Default is 2.
    verbose : bool, optional
        If True, prints the computed batch size. Default is True.

    Returns
    -------
    int
        Computed batch size.
    """
    M, N = data.shape

    batch_size = min(M // divisor, max_size)

    if batch_size < min_size:
        batch_size = min_size

    if batch_size % 2 != 0:
        batch_size = batch_size + 1
        if batch_size > max_size:
            batch_size = batch_size -2

    if batch_size < 1:
        batch_size = 1

    if verbose:
        print(f"(Batch-Size) Computed batch size: {batch_size}\n")

    return batch_size



def load_share_lib(name="lib_cdrec", lib=True, verbose=True):
    """
    Load the shared library based on the operating system.

    Parameters
    ----------
    name : str, optional
        The name of the shared library (default is "lib_cdrec").
    lib : bool, optional
        If True, the function loads the library from the default 'imputegap' path; if False, it loads from a local path (default is True).
    verbose : bool, optional
        Whether to display the contamination information (default is True).

    Returns
    -------
    ctypes.CDLL
        The loaded shared library object.
    """
    system = platform.system()
    if system == "Windows":
        ext = ".so"
    elif system == "Darwin":
        ext = ".dylib"  # macOS uses .dylib for dynamic libraries
    else:
        ext = ".so"

    if lib:
        lib_path = importlib.resources.files('imputegap.algorithms.lib').joinpath("./" + str(name) + ext)
    else:
        local_path_lin = './algorithms/lib/' + name + ext

        if not os.path.exists(local_path_lin):
            local_path_lin = './imputegap/algorithms/lib/' + name + ext

        lib_path = os.path.join(local_path_lin)

    if verbose:
        print("\n(SYS) Wrapper files loaded for C++ : ", lib_path, "\n")

    return ctypes.CDLL(lib_path)



def list_of_algorithms():
    return sorted([
        "CDRec",
        "IterativeSVD",
        "GROUSE",
        "ROSL",
        "SPIRIT",
        "SoftImpute",
        "SVT",
        "TRMF",
        "STMVL",
        "DynaMMo",
        "TKCM",
        "IIM",
        "XGBOOST",
        "MICE",
        "MissForest",
        "KNNImpute",
        "Interpolation",
        "MinImpute",
        "MeanImpute",
        "ZeroImpute",
        "MeanImputeBySeries",
        "MRNN",
        "BRITS",
        "DeepMVI",
        "MPIN",
        "PRISTI",
        "MissNet",
        "GAIN",
        "GRIN",
        "BayOTIDE",
        "HKMF_T",
        "BitGraph",
        "NuwaTS",
        "GPT4TS"
    ])

def list_of_patterns():
    return sorted([
        "aligned",
        "disjoint",
        "overlap",
        "scattered",
        "mcar",
        "gaussian",
        "distribution"
    ])

def list_of_datasets(txt=False):

    list = sorted([
        "airq",
        "bafu",
        "chlorine",
        "climate",
        "drift",
        "eeg-alcohol",
        "eeg-reading",
        "electricity",
        "fmri-stoptask",
        "forecast-economy",
        "meteo",
        "motion",
        "soccer",
        "solar-plant",
        "sport-activity",
        "stock-exchange",
        "temperature",
        "traffic"
    ])

    if txt:
        list = [dataset + ".txt" for dataset in list]

    return list



def list_of_optimizers():
    return sorted([
        "ray_tune",
        "bayesian",
        "particle_swarm",
        "successive_halving",
        "greedy"
    ])

def list_of_downstreams():
    return sorted(list_of_downstreams_sktime() + list_of_downstreams_darts())


def list_of_downstreams_sktime():
    return sorted([
        "prophet",
        "exp-smoothing",
        "hw-add",
        "arima",
        "sf-arima",
        "bats",
        "ets",
        "croston",
        "theta",
        "unobs",
        "naive"
    ])

def list_of_downstreams_darts():
    return sorted([
        "nbeats",
        "xgboost",
        "lightgbm",
        "lstm",
        "deepar",
        "transformer"
    ])

def list_of_extractors():
    return sorted([
        "pycatch",
        "tsfel",
        "tsfresh"
    ])

def list_of_families():
    return sorted(["DeepLearning", "MatrixCompletion", "PatternSearch", "MachineLearning", "Statistics", "LLMs"])

def list_of_metrics():
    return ["RMSE", "MAE", "MI", "CORRELATION", "RUNTIME", "RUNTIME_LOG"]

def list_of_algorithms_with_families():
    return sorted([
        "MatrixCompletion.CDRec",
        "MatrixCompletion.IterativeSVD",
        "MatrixCompletion.GROUSE",
        "MatrixCompletion.ROSL",
        "MatrixCompletion.SPIRIT",
        "MatrixCompletion.SoftImpute",
        "MatrixCompletion.SVT",
        "MatrixCompletion.TRMF",
        "PatternSearch.STMVL",
        "PatternSearch.DynaMMo",
        "PatternSearch.TKCM",
        "MachineLearning.IIM",
        "MachineLearning.XGBOOST",
        "MachineLearning.MICE",
        "MachineLearning.MissForest",
        "Statistics.KNNImpute",
        "Statistics.Interpolation",
        "Statistics.MinImpute",
        "Statistics.MeanImpute",
        "Statistics.ZeroImpute",
        "Statistics.MeanImputeBySeries",
        "DeepLearning.MRNN",
        "DeepLearning.BRITS",
        "DeepLearning.DeepMVI",
        "DeepLearning.MPIN",
        "DeepLearning.PRISTI",
        "DeepLearning.MissNet",
        "DeepLearning.GAIN",
        "DeepLearning.GRIN",
        "DeepLearning.BayOTIDE",
        "DeepLearning.HKMF_T",
        "DeepLearning.BitGraph",
        "LLMs.NuwaTS",
        "LLMs.GPT4TS"
    ])

def list_of_normalizers():
    return ["z_score", "min_max"]

