// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Dario Mangoni, Radu Serban
// =============================================================================

#include <bitset>
#include <mpi.h>

#include "chrono_mumps/ChMumpsEngine.h"

namespace chrono {

ChMumpsEngine::ChMumpsEngine() {

    //std::cout << "Initializing MPI for Mumps Engine" << std::endl;

    //// Initialize MPI
    //int argc;
    //char** argv;
    //int myid;
    //std::cout << "Calling MPI_Init" << std::endl;

    //int MPI_Init_ret = MPI_Init(&argc, &argv);
    //std::cout << "Called MPI_Init" << std::endl;

    //if (MPI_Init_ret!=MPI_SUCCESS) {
    //    std::cerr << "MPI_Init failed with error code: " << MPI_Init_ret << std::endl;
    //}
    //int MPI_Comm_rank_ret = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    //if (MPI_Comm_rank_ret!=MPI_SUCCESS) {
    //    std::cerr << "MPI_Comm_rank failed with error code: " << MPI_Comm_rank_ret << std::endl;
    //}
    //if (myid == 0) {
    //    std::cout << "MPI_Init successful. Current ID: " << myid << std::endl;
    //}

    /* Initialize a MUMPS instance. Use MPI_COMM_WORLD */
    mumps_id.job = INIT;
    mumps_id.par = 1;
    mumps_id.sym = UNSYMMETRIC;
    mumps_id.comm_fortran = USE_COMM_WORLD;

    dmumps_c(&mumps_id);

    /* Output messages */
    mumps_id.ICNTL(1) = 6;  // Error
    mumps_id.ICNTL(2) = 0;  // Diagnostic and warnings
    mumps_id.ICNTL(3) = 0;  // Global information
    mumps_id.ICNTL(4) = 1;  // Error, warning and diagnostic control

    /* Matrix control */
    mumps_id.ICNTL(5) = 0;   // COO Matrix format selection
    mumps_id.ICNTL(18) = 0;  // Matrix centralized on the host
}

ChMumpsEngine::~ChMumpsEngine() {
    mumps_id.job = END;
    dmumps_c(&mumps_id);  // Terminate instance

    //std::cout << "Finalize MPI instance" << std::endl;
    //// Terminate MPI instance
    //int MPI_Finalize_ret = MPI_Finalize();
    //if (MPI_Finalize_ret!=MPI_SUCCESS) {
    //    std::cerr << "MPI_Finalize failed with error code: " << MPI_Finalize_ret << std::endl;
    //}
    //std::cout << "MPI instance finalized" << std::endl;

}

void ChMumpsEngine::SetProblem(const ChSparseMatrix& Z, ChVectorRef rhs) {
    SetMatrix(Z);
    SetRhsVector(rhs);
}

void ChMumpsEngine::SetMatrix(const ChSparseMatrix& Z) {
    // Convert to COO representation (1-indexed)
    int dim = (int)Z.rows();
    int nnz = (int)Z.nonZeros();

    m_irn.resize(nnz);
    m_jcn.resize(nnz);
    m_a.resize(nnz);

    int i = 0;
    for (int k = 0; k < Z.outerSize(); ++k) {
        for (ChSparseMatrix::InnerIterator it(Z, k); it; ++it) {
            m_irn[i] = (int)it.row() + 1;
            m_jcn[i] = (int)it.col() + 1;
            m_a[i] = it.value();
            i++;
        }
    }

    mumps_id.n = dim;
    mumps_id.nz = nnz;
    mumps_id.irn = m_irn.data();
    mumps_id.jcn = m_jcn.data();
    mumps_id.a = m_a.data();
}

void ChMumpsEngine::SetMatrixSymmetry(mumps_SYM mat_type) {
    mumps_id.sym = mat_type;
}

void ChMumpsEngine::SetRhsVector(ChVectorRef b) {
    mumps_id.rhs = b.data();
}

void ChMumpsEngine::SetRhsVector(double* b) {
    mumps_id.rhs = b;
}

void ChMumpsEngine::EnableNullPivotDetection(bool val, double threshold) {
    mumps_id.ICNTL(24) = val;      ///< activates null pivot detection
    mumps_id.ICNTL(25) = 0;        ///< tries to compute one of the many solutions of AX = B
    mumps_id.CNTL(5) = 1e20;       ///< fixation value
    mumps_id.CNTL(3) = threshold;  ///< pivot threshold
}


int ChMumpsEngine::MumpsCall(mumps_JOB job_call) {
    /* Call the MUMPS package. */
    mumps_id.job = job_call;
    dmumps_c(&mumps_id);
    return mumps_id.INFOG(1);
}

void ChMumpsEngine::PrintINFOG() {
    if (mumps_id.INFOG(1) > 0) {
        std::cout << "WARN: INFOG(1)=";
        typedef std::bitset<sizeof(int)> IntBits;
        if (IntBits(mumps_id.INFOG(1)).test(0))
            std::cout << "+1: Row or column index out of range. Faulty entries: " << mumps_id.INFOG(2) << std::endl;
        if (IntBits(mumps_id.INFOG(1)).test(1))
            std::cout << "+2: Solution max-norm = 0" << std::endl;
        if (IntBits(mumps_id.INFOG(1)).test(3))
            std::cout << "+8: More than " << mumps_id.ICNTL(10) << " iterative refinements are required." << std::endl;
        return;
    }

    if (mumps_id.INFOG(1) == 0) {
        std::cout << "INFOG(1)=0: Mumps is successful!" << std::endl;
        return;
    }

    std::cout << "ERR: INFOG(1)=";
    switch (mumps_id.INFOG(1)) {
        case (0):
            std::cout << "0: Mumps is successful!" << std::endl;
            break;
        case (-1):
            std::cout << "-1: Error on processor " << mumps_id.INFOG(2) << std::endl;
            break;
        case (-2):
            std::cout << "-2: Number of nonzeros out of range NZ=" << mumps_id.INFOG(2) << std::endl;
            break;
        case (-3):
            std::cout << "-3: Mumps called with wrong JOB. JOB=" << mumps_id.INFOG(2) << std::endl;
            break;
        case (-4):
            std::cout << "-4: Error in user-provided permutation array PERM_IN at position: " << mumps_id.INFOG(2) << std::endl;
            break;
        case (-5):
            std::cout << "-5: Problem of real workspace allocation of size " << mumps_id.INFOG(2) << " during analysis" << std::endl;
            break;
        case (-6):
            std::cout << "-6: Matrix is singular in structure. Matrix rank: " << mumps_id.INFOG(2) << std::endl;
            break;
        case (-7):
            std::cout << "-7: Problem of integer workspace allocation of size " <<  mumps_id.INFOG(2) << " during analysis" << std::endl;
            break;
        case (-10):
            std::cout << "-10: Matrix is numerically singular." << std::endl;
            break;
        case (-16):
            std::cout << "-16: N is out of range. N=" << mumps_id.INFOG(2) << std::endl;
            break;
        case (-21):
            std::cout << "-21: PAR=1 not allowed because only one processor is available." << std::endl;
            break;
        case (-22):
            std::cout << "-22: Array pointers have problems. INFOG(2)=" << mumps_id.INFOG(2) << std::endl;
            break;
        default:
            std::cout << mumps_id.INFOG(1) << ": See the user guide. INFOG(2)=" << mumps_id.INFOG(2) << std::endl;
    }
}

}  // namespace chrono
