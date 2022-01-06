//
// Created by ubuntu on 1/6/22.
//
#include "ThermoElastic.h"

void right_hand_side(const std::vector<Point<2>> &points,
                     std::vector<Tensor<1, 2>> &  values)
{
    AssertDimension(values.size(), points.size());
    Point<2> point_1, point_2;
    point_1(0) = 0.5;
    point_2(0) = -0.5;
    for (unsigned int point_n = 0; point_n < points.size(); ++point_n)
    {        if (((points[point_n] - point_1).norm_square() < 0.2 * 0.2) ||
                 ((points[point_n] - point_2).norm_square() < 0.2 * 0.2))
            values[point_n][0] = 1.0;
        else
            values[point_n][0] = 0.0;
        if (points[point_n].norm_square() < 0.2 * 0.2)
            values[point_n][1] = 1.0;
        else
            values[point_n][1] = 0.0;
    }
}

void ThermoElastic::assemble_system() {
    FEValuesExtractors::Scalar u(0);
    FEValuesExtractors::Scalar v(1);
    FEValuesExtractors::Scalar T(2);

    QGauss<2> quadrature_formula(fe.degree + 1);
    FEValues<2> fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);
    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    std::vector<double> lambda_values(n_q_points);
    std::vector<double> mu_values(n_q_points);
    Functions::ConstantFunction<2> lambda(1.), mu(1.);
    std::vector<Tensor<1, 2>> rhs_values(n_q_points);

    //Old time step iteration
    std::vector<Tensor<1, 2>> old_sol_grad_u(n_q_points);
    std::vector<double> old_sol_val_u(n_q_points);
    std::vector<Tensor<1, 2>> old_sol_grad_v(n_q_points);
    std::vector<double> old_sol_val_v(n_q_points);
    std::vector<Tensor<1, 2>> old_sol_grad_T(n_q_points);
    std::vector<double> old_sol_val_T(n_q_points);

    system_matrix = 0;
    system_rhs = 0;

    for (const auto &cell : dof_handler.active_cell_iterators()){
        cell_matrix = 0;
        cell_rhs    = 0;
        fe_values.reinit(cell);
        lambda.value_list(fe_values.get_quadrature_points(), lambda_values);
        mu.value_list(fe_values.get_quadrature_points(), mu_values);
        right_hand_side(fe_values.get_quadrature_points(), rhs_values);

        fe_values[u].get_function_values(old_solution, old_sol_val_u);
        fe_values[u].get_function_gradients(old_solution, old_sol_grad_u);
        fe_values[v].get_function_values(old_solution, old_sol_val_v);
        fe_values[v].get_function_gradients(old_solution, old_sol_grad_v);
        fe_values[T].get_function_values(old_solution, old_sol_val_T);
        fe_values[T].get_function_gradients(old_solution, old_sol_grad_T);

        for (unsigned int q = 0; q < n_q_points; ++q){
            for (unsigned int i = 0; i < dofs_per_cell; ++i){
                for(unsigned int j = 0; j < dofs_per_cell; ++j){
                    double psiU = fe_values[u].value(i,q);
                    double psiV = fe_values[v].value(i,q);
                    double psiT = fe_values[T].value(i,q);
                    double grad_U_x = fe_values[u].gradient(j,q)[0];
                    double grad_U_y = fe_values[u].gradient(j,q)[1];
                    double grad_V_x = fe_values[v].gradient(j,q)[0];
                    double grad_V_y = fe_values[v].gradient(j,q)[1];
                    double grad_T_x = fe_values[T].gradient(j,q)[0];
                    double grad_T_y = fe_values[T].gradient(j,q)[1];
                    double grad_psiU_x = fe_values[u].gradient(i,q)[0];
                    double grad_psiU_y = fe_values[u].gradient(i,q)[1];
                    double grad_psiV_x = fe_values[v].gradient(i,q)[0];
                    double grad_psiV_y = fe_values[v].gradient(i,q)[1];
                    double t1 = mu_values[q]*(grad_U_x*grad_psiU_x+grad_U_y*grad_psiU_y + grad_V_x*grad_psiV_x+grad_V_y*grad_psiV_y);
                    double t2 = mu_values[q]*(grad_U_x*grad_psiU_x+grad_U_y*grad_psiV_x + grad_V_x*grad_psiU_y+grad_V_y*grad_psiV_y);
                    double t3 = lambda_values[q]*(grad_U_x*grad_psiU_x+grad_U_x*grad_psiV_y + grad_V_y*grad_psiU_x+grad_V_y*grad_psiV_y);
                    double gamm = alpha*(3*lambda_values[q]+2*mu_values[q]);
                    double t4 = (grad_T_x*psiU+grad_T_y*psiV)*gamm;
                    double Eelast = t1+t2+t3+t4;
                    double K = fe_values[T].gradient(i,q)*fe_values[T].gradient(j,q);
                    double M = fe_values[T].value(i,q)*fe_values[T].value(j,q)/kappa;
                    double eta = gamm*Tsf/k;
                    double P = eta*fe_values[T].value(i,q)*fe_values[u].gradient(j,q)[0];
                    double Q = eta*fe_values[T].value(i,q)*fe_values[v].gradient(j,q)[1];
                    double Eheat = time_step*theta*K + M + P + Q;
                    cell_matrix(i, j) += (Eelast + Eheat)*fe_values.JxW(q);
                }
                //RHS
                double EelastR = fe_values[u].value(i,q)*rhs_values[q][0]+fe_values[v].value(i,q)*rhs_values[q][1];
                double Kn = fe_values[T].gradient(i,q)*old_sol_grad_T[q];
                double Mn = fe_values[T].value(i,q)*old_sol_val_T[q]/kappa;
                double gamm = alpha*(3*lambda_values[q]+2*mu_values[q]);
                double eta = gamm*Tsf/k;
                double Pn = fe_values[T].value(i,q)*old_sol_grad_u[q][0]*eta;
                double Qn = fe_values[T].value(i,q)*old_sol_grad_v[q][1]*eta;
                double EheatR = -time_step*(1-theta)*Kn + Mn + Pn + Qn;
                cell_rhs(i) += (EelastR+EheatR)*fe_values.JxW(q);
            }
        }
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
                system_matrix.add(local_dof_indices[i],
                                    local_dof_indices[j],
                                    cell_matrix(i, j));
            system_rhs(local_dof_indices[i]) += cell_rhs(i);
        }
    }
    //Applying zero BC
    std::map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             55,
                                             Functions::ZeroFunction<2>(3),
                                             boundary_values);
    MatrixTools::apply_boundary_values(boundary_values,
                                       system_matrix,
                                       solution,
                                       system_rhs);
}
