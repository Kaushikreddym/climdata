//
// Created by zakhar on 01/11/18.
//

#pragma once

#include <cstdlib>
#include <cstdint>
#include <tuple>

#ifdef MLPACK_ACTIVE
#include <mlpack/core.hpp>
#else
#include <armadillo>
#endif


arma::mat
marshal_as_arma(double *matrixNative, size_t dimN, size_t dimM);

void
marshal_as_native(const arma::mat &matrixArma, double *container);

void
marshal_as_failed(double *container, size_t dimN, size_t dimM);

void
verifyRecovery(arma::mat &mat);

extern "C"
{
void
cd_decomposition(
        double *matrixNative, size_t dimN, size_t dimM,
        double *loadContainer, double *relContainer,
        size_t truncation
);

void
cdrec_imputation_simple(
        double *matrixNative, size_t dimN, size_t dimM,
        size_t truncation
);

void
cdrec_imputation_parametrized(
        double *matrixNative, size_t dimN, size_t dimM,
        size_t truncation, double epsilon, size_t iters
);

void
stmvl_imputation_simple(
        double *matrixNative, size_t dimN, size_t dimM
);

void
stmvl_imputation_parametrized(
        double *matrixNative, size_t dimN, size_t dimM,
        size_t window_size, double gamma, double alpha
);

}