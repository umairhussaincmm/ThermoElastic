//
// Created by ubuntu on 1/6/22.
//

#include "ThermoElastic.h"

ThermoElastic::ThermoElastic()
        : dof_handler(triangulation)
        , fe(FE_Q<2>(1), 3)
        , time(0.0)
        , final_time(1)
        , time_step(.001)
        , timestep_number(0)
        , theta(0.5)
        , alpha(0.005)
        , k(10)
        , Tsf(0)
        , rho(1000)
        , ce(.02)
        , kappa(k/(rho*ce))
{}

