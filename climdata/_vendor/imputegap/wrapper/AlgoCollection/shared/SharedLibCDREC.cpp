//
// Created by nater quentin on 24/09/24
//

#include <cassert>
#include <chrono>

#include "SharedLibCDREC.h"

// cd as a factorization algo
#include "../Algebra/CentroidDecomposition.h"
#include "../Algebra/Auxiliary.h"

// imputation algorithms
#include "../Algorithms/CDMissingValueRecovery.h"


#ifdef MLPACK_ACTIVE
#include "../Algorithms/NMFMissingValueRecovery.h"
#endif

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

// Exposed functions

extern "C"
{

// cd matrix factorization

void
cd_decomposition(
        double *matrixNative, size_t dimN, size_t dimM,
        double *loadContainer, double *relContainer,
        size_t truncation
)
{
    arma::mat input = marshal_as_arma(matrixNative, dimN, dimM);

    arma::mat Load;
    arma::mat Rel;

    std::tie(Load, Rel) = Algorithms::CentroidDecomposition::PerformCentroidDecomposition(input, truncation);

    // [!] explicitly copy the output
    marshal_as_native(Load, loadContainer);
    marshal_as_native(Rel, relContainer);
}

void
cdrec_imputation_parametrized(
        double *matrixNative, size_t dimN, size_t dimM,
        size_t truncation, double epsilon, size_t iters
)
{
    arma::mat input = marshal_as_arma(matrixNative, dimN, dimM);

    // Local
    int64_t result;
    Algorithms::CDMissingValueRecovery rmv(input, iters, epsilon);
    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;

    // Recovery
    rmv.setReduction(truncation);
    rmv.disableCaching = false;
    rmv.useNormalization = false;

    begin = std::chrono::steady_clock::now();
    rmv.autoDetectMissingBlocks();
    rmv.performRecovery(truncation == 0);
    end = std::chrono::steady_clock::now();

    result = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();

    verifyRecovery(input);

    (void)result;

    // [!] input already modified
}
}