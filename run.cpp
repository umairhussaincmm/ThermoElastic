//
// Created by ubuntu on 1/6/22.
//
#include "ThermoElastic.h"

void ThermoElastic::run() {
    setup_system();
    InitialValues initial_value;
    VectorTools::interpolate(dof_handler,
                             initial_value,
                             old_solution);
    solution = old_solution;
    output_results();

    //Time steps begin here:
    unsigned int timestep_number = 1;
    for (; time <= final_time; time += time_step, ++timestep_number){
        assemble_system();
        solve();
        output_results();

        //Saving the solution
        old_solution = solution;
    }
}
