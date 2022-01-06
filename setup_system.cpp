//
// Created by ubuntu on 1/6/22.
//
#include "ThermoElastic.h"

void ThermoElastic::setup_system() {
    const Point<2,double> p1(0,0);
    const Point<2,double> p2(12,3);
    GridGenerator::hyper_rectangle(triangulation, p1, p2);
    triangulation.refine_global(4); //refining

    for (auto &cell : triangulation.active_cell_iterators())
        for (unsigned int face_n = 0;
             face_n < GeometryInfo<2>::faces_per_cell;
             ++face_n)
            if (cell->at_boundary(face_n)){
                if (std::fabs(cell->face(face_n)->center()[0] - 0) < 1e-10 )
                    cell->face(face_n)->set_boundary_id(55);

                else if (std::fabs(cell->face(face_n)->center()[0] - 12) < 1e-10 )
                    cell->face(face_n)->set_boundary_id(55);
            }

    dof_handler.distribute_dofs(fe);
    solution.reinit(dof_handler.n_dofs());
    old_solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp);
    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit(sparsity_pattern);
}
