//
// Created by ubuntu on 1/6/22.
//
#include "ThermoElastic.h"

void ThermoElastic::output_results(const unsigned int timestep_number) const{
    DataOut<2> data_out;
    data_out.attach_dof_handler(dof_handler);
    std::vector<std::string> solution_names;

    solution_names.emplace_back("x_displacement");
    solution_names.emplace_back("y_displacement");
    solution_names.emplace_back("Temperature");
    //std::cout<< "Output Number: " << timestep_number << std::endl;
    data_out.add_data_vector(solution, solution_names);
    data_out.build_patches();
    std::ofstream output("results/trial2/solution-" + Utilities::int_to_string(timestep_number, 3) + ".vtk");
    data_out.write_vtk(output);
}
