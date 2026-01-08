//
// Created by zakhar on 01/11/18.
//

#include <cassert>
#include <chrono>

#include "SharedLibFunctions.h"

// cd as a factorization algo
#include "../Algebra/CentroidDecomposition.h"

// imputation algorithms
#include "../Algorithms/CDMissingValueRecovery.h"
#include "../Algorithms/TKCM.h"
#include "../Algorithms/ST_MVL.h"
#include "../Algorithms/SPIRIT.h"
#include "../Algorithms/GROUSE.h"
#include "../Algorithms/DynaMMo.h"
#include "../Algorithms/SVT.h"
#include "../Algorithms/ROSL.h"
#include "../Algorithms/IterativeSVD.h"
#include "../Algorithms/SoftImpute.h"

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
cdrec_imputation_simple(
        double *matrixNative, size_t dimN, size_t dimM,
        size_t truncation
)
{
    arma::mat input = marshal_as_arma(matrixNative, dimN, dimM);
    
    // Local
    int64_t result;
    Algorithms::CDMissingValueRecovery rmv(input);
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

void
stmvl_imputation_simple(
        double *matrixNative, size_t dimN, size_t dimM
)
{
    arma::mat input = marshal_as_arma(matrixNative, dimN, dimM);
    
    // Local
    int64_t result;

    Algorithms::ST_MVL stmvl(input, 2.0, 0.85, 7);

    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;

    // Recovery

    begin = std::chrono::steady_clock::now();
    stmvl.Run();
    end = std::chrono::steady_clock::now();

    result = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
    
    verifyRecovery(input);
    
    (void)result;
    
    // [!] input already modified
}

void
stmvl_imputation_parametrized(
        double *matrixNative, size_t dimN, size_t dimM,
        size_t window_size, double gamma, double alpha
)
{
    arma::mat input = marshal_as_arma(matrixNative, dimN, dimM);
    
    // Local
    int64_t result;

    Algorithms::ST_MVL stmvl(input, alpha, gamma, window_size);

    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;

    // Recovery

    begin = std::chrono::steady_clock::now();
    stmvl.Run();
    end = std::chrono::steady_clock::now();

    result = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
    
    verifyRecovery(input);
    
    (void)result;
    
    // [!] input already modified
}

}