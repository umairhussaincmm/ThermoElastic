//
// Created by ubuntu on 1/6/22.
//

#ifndef THERMO_ELASTIC_THERMOELASTIC_H
#define THERMO_ELASTIC_THERMOELASTIC_H

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <fstream>
#include <iostream>

using namespace dealii;

class ThermoElastic {
public:
  ThermoElastic();
  void run();
private:
  void setup_system();
  void assemble_system();
  void solve();
  void refine_grid();
  void output_results(const unsigned int cycle) const;
  Triangulation<2> triangulation;
  DoFHandler<2>    dof_handler;
  FESystem<2> fe;
  AffineConstraints<double> constraints;
  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;
  Vector<double> solution, old_solution;
  Vector<double> system_rhs;
  double       time;
  const double final_time, time_step;
  const double theta;
  const double alpha, k, Tsf, rho, ce, kappa;
};


#endif //THERMO_ELASTIC_THERMOELASTIC_H
