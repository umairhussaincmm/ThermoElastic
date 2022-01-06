//
// Created by ubuntu on 1/6/22.
//
#include "ThermoElastic.h"

void InitialValues::vector_value(const Point<2> &p,
                                 Vector<double> & values) const
{
    double u_initial = 0.;
    double v_initial = 0.;
    double T_initial = 0.;

    values(0)= u_initial;
    values(1)= v_initial;
    values(2)= T_initial;
}
