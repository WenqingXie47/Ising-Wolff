#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <cmath>


#include "monteCarlo.h"
#include "isingState.h"
#include "util.h"

int dim=2;
int length = 16;
double beta=0.5;
double h=0;
double J=1;

int n_samples=1000;
int n_iters_per_sample=100;
int n_thermalization_iters=10000;



double mean(std::vector<double> list);
double mean_abs(std::vector<double> list);
std::map<std::string, std::vector<double>> measure(std::vector<double>& beta_list);
void write_result(std::string fname, std::map<std::string, std::vector<double>> measurement);


int main()
{   
    std::vector<double> beta_list {};
    for (double beta = 0.1; beta< 1.0; beta+=0.02){
        beta_list.push_back(beta);
    }

    auto measurement = measure(beta_list);
    std::string fname {"../Data/data.csv"};
    write_result(fname, measurement);
    return 0;
}


std::map<std::string, std::vector<double>> measure(std::vector<double>& beta_list)
{
    std::vector<double> energy_list{};
    std::vector<double> magnetization_list{};
    std::vector<double> heat_capacity_list{};

    for (double beta : beta_list){
        ising::IsingState ising_state{dim, length, beta, h, J};
        ising::IsingWolff ising_mc(ising_state);
        auto history = ising_mc.sample(n_samples, n_iters_per_sample, n_thermalization_iters);

         

        energy_list.push_back(ising::mean(history["energy"])/ising_state.n_sites);
        magnetization_list.push_back(ising::mean_abs(history["magnetization"]));
        heat_capacity_list.push_back(ising::variance(history["energy"])*beta*beta/ising_state.n_sites);
    }
    std::map<std::string, std::vector<double>> measurement {
        {"beta", beta_list}, 
        {"energy", energy_list},
        {"magnetization", magnetization_list},
        {"heat_capacity", heat_capacity_list},
    };
    return measurement;
}


void write_result(std::string fname, std::map<std::string, std::vector<double>> measurement)
{
    std::ofstream myfile;
    myfile.open (fname);
    myfile << "beta, energy, magnetization, heat_capacity" << std::endl;

    for (int i=0; i<measurement["beta"].size(); i++){
        myfile << measurement["beta"][i] << ", " 
            << measurement["energy"][i] << ", " 
            << measurement["magnetization"][i] << ", " 
            << measurement["heat_capacity"][i] << ", "
            << std::endl;
    }
    myfile.close();
}

