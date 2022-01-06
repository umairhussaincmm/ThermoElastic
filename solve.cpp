//
// Created by ubuntu on 1/6/22.
//
#include "ThermoElastic.h"

void ThermoElastic::solve() {
    SolverControl            solver_control(1000, 1e-12);
    SolverCG<Vector<double>> cg(solver_control);
    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);
    cg.solve(system_matrix, solution, system_rhs, preconditioner);
}
