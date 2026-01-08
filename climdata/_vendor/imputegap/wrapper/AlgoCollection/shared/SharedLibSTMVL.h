//
// Created by nater quentin 24/09/24
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
stmvl_imputation_parametrized(
        double *matrixNative, size_t dimN, size_t dimM,
        size_t window_size, double gamma, double alpha
);


}