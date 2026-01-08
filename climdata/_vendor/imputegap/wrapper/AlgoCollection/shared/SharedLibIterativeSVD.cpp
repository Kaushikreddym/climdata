//
// Created by quentin nater on 12/01/25.
//

#include <cassert>
#include <chrono>

#include "SharedLibIterativeSVD.h"

#include "../Algebra/RSVD.h"
#include "../Algorithms/SVT.h"
#include "../Algorithms/IterativeSVD.h"


// Non-exposed functions
arma::mat
marshal_as_arma(double *matrixNative, size_t dimN, size_t dimM)
{
    return arma::mat(matrixNative, dimN, dimM, false, true);
}

void
marshal_as_native(const arma::mat &matrixArma, double *container)
{
    std::copy_n(matrixArma.memptr(), matrixArma.n_rows * matrixArma.n_cols, container);
}

void
marshal_as_failed(double *container, size_t dimN, size_t dimM)
{
    std::fill_n(container, dimN * dimM, std::numeric_limits<double>::quiet_NaN());
}

void
verifyRecovery(arma::mat &mat)
{
    for (uint64_t j = 0; j < mat.n_cols; ++j)
    {
        for (uint64_t i = 0; i < mat.n_rows; ++i)
        {
            if (std::isnan(mat.at(i, j)))
            {
                mat.at(i, j) = std::sqrt(std::numeric_limits<double>::max() / 100000.0);
            }
        }
    }
}

bool
converged(const arma::mat &X_old, const arma::mat &X_new, const std::vector<arma::uvec> &indices, double threshold)
{
    double delta = 0.0;
    double old_norm = 0.0;
    for (uint64_t i = 0; i < X_old.n_cols; ++i)
    {
        for (uint64_t j : indices[i])
        {
            old_norm += X_old.at(j, i) * X_old.at(j, i);
            double diff = X_old.at(j, i) - X_new.at(j, i);
            delta += diff * diff;
        }
    }

    return old_norm > arma::datum::eps && (delta / old_norm) < threshold;
}

// Exposed functions

extern "C"
{

void
recoveryIterativeSVD(double *matrixNative, size_t dimN, size_t dimM, size_t rank)
{
    arma::mat input = marshal_as_arma(matrixNative, dimN, dimM);

    // Local
    int64_t result;

    Algorithms::IterativeSVD iterativeSVD;

    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;

    // Recovery
    begin = std::chrono::steady_clock::now();
    iterativeSVD.recoveryIterativeSVD(input, rank);
    end = std::chrono::steady_clock::now();

    result = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();

    verifyRecovery(input);

    (void)result;
}

}