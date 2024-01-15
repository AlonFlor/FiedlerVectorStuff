
#include <iostream>
#include <fstream>
#include <cmath>

#include "falling_constrained.hpp"
#include "falling.hpp"
#include "side_motion_constrained.hpp"


int main(int argc, char *argv[])
{
    int num_times = 0;
    std::string constitutive_model = argv[1];
    std::string scenario_name = argv[2];
    int number_of_time_steps = std::stoi(argv[3]);
    for(int i=4; i<argc; ++i){
        if(i%2==0){
            num_times=std::stoi(argv[i]);
         }else{
           printf("Simulating %s for %d times.\n", argv[i], num_times);

            //run the program and get the timings
            std::vector<double> timings;
            double average_time = 0.;
            if(!scenario_name.compare(0, scenario_name.size(), "falling_constrained")){
                for(int j=0; j<num_times; ++j){
                    double time_elapsed_one_run = run_one_sim_falling_constrained(argv[i], j, constitutive_model, number_of_time_steps);
                    timings.push_back(time_elapsed_one_run);
                    average_time += time_elapsed_one_run;
                }
            }else if(!scenario_name.compare(0, scenario_name.size(), "falling")){
                for(int j=0; j<num_times; ++j){
                    double time_elapsed_one_run = run_one_sim_falling(argv[i], j, constitutive_model, number_of_time_steps);
                    timings.push_back(time_elapsed_one_run);
                    average_time += time_elapsed_one_run;
                }
            }else if(!scenario_name.compare(0, scenario_name.size(), "side_motion_constrained")){
                for(int j=0; j<num_times; ++j){
                    double time_elapsed_one_run = run_one_sim_side_motion_constrained(argv[i], j, constitutive_model, number_of_time_steps);
                    timings.push_back(time_elapsed_one_run);
                    average_time += time_elapsed_one_run;
                }
            }
            average_time /= (double)num_times;

            //get standard deviation of timings
            double std_dev_time = 0.;
            for(int j=0; j<num_times; ++j){
                double diff = timings[j] - average_time;
                std_dev_time += diff*diff;
            }
            std_dev_time /= (double)num_times;
            std_dev_time = sqrt(std_dev_time);

            //print timings
            char file_path[70];
            std::string file_name = argv[i];
            file_name = file_name.substr(0,file_name.size()-4) + "_" + scenario_name + "_" + constitutive_model;
            sprintf(file_path, "timings/%s_timings.txt", file_name.data());
            FILE* fp = fopen(file_path, "w");
            for(int j=0; j<num_times; ++j){
                fprintf(fp, "%lf\n", timings[j]);
            }
            fprintf(fp, "\n\nAverage: %lf\nStandard deviation: %lf\n", average_time, std_dev_time);
            fclose(fp);

        }      
    }

    return 0;
}
