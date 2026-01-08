//
// Created by quentin nater on 12/01/25.
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

bool
converged(const arma::mat &X_old, const arma::mat &X_new, const std::vector<arma::uvec> &indices, double threshold);


extern "C"
{
void recoveryIterativeSVD(double *matrixNative, size_t dimN, size_t dimM, size_t rank);
}