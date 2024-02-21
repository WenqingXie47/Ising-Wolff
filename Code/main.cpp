#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <algorithm>

#include "monteCarlo.h"
#include "isingModel.h"
#include "util.h"

int dim=2;
int length = 16;
double beta=0.5;
double h=0.0;
double J=1.0;

int n_samples=5000;
int n_iters_per_sample=500;
int n_thermalization_iters=1000000;



std::vector<double> generate_beta_list();
std::map<std::string, std::vector<double>> measure(std::vector<double>& beta_list);
void write_result(std::string fname, std::map<std::string, std::vector<double>> measurement);


int main()
{   
    std::vector<double> beta_list = generate_beta_list();
    auto measurement = measure(beta_list);
    std::string fname {"../Data/data.csv"};
    write_result(fname, measurement);
    return 0;
}


std::vector<double> generate_beta_list()
{
    std::vector<double> beta_list {};
    for (double beta = 0.1; beta< 1.0; beta+=0.02){
        beta_list.push_back(beta);
    }
    for (double beta = 0.3; beta<= 0.5; beta+=0.005){
        beta_list.push_back(beta);
    }

    // sorted beta list
    std::sort(beta_list.begin(), beta_list.end());
    // remove repeated beta
    auto last = std::unique(beta_list.begin(), beta_list.end());
    beta_list.erase(last, beta_list.end());
    return beta_list;
}

std::map<std::string, std::vector<double>> measure(std::vector<double>& beta_list)
{
    std::vector<double> energy_list{};
    std::vector<double> magnetization_list{};
    std::vector<double> m2_list{};
    std::vector<double> heat_capacity_list{};

    ising::GridState grid{dim,length};

    for (double beta : beta_list){
        ising::IsingModel ising_model{grid, beta, J, h};
        ising::IsingWolff ising_mc(ising_model);
        auto history = ising_mc.sample(n_samples, n_iters_per_sample, n_thermalization_iters);

         
        energy_list.push_back(ising::mean(history["energy"])/ising_model.grid.n_sites);
        magnetization_list.push_back(ising::mean_abs(history["magnetization"]));
        m2_list.push_back(ising::mean_squared(history["magnetization"]));
        heat_capacity_list.push_back(ising::variance(history["energy"])*beta*beta/ising_model.grid.n_sites);
    }
    std::map<std::string, std::vector<double>> measurement {
        {"beta", beta_list}, 
        {"magnetization", magnetization_list},
        {"m2", m2_list},
        {"energy", energy_list},
        {"heat_capacity", heat_capacity_list},
    };
    return measurement;
}


void write_result(std::string fname, std::map<std::string, std::vector<double>> measurement)
{
    std::ofstream myfile;
    myfile.open (fname);
    myfile << "beta, magnetization, m2, energy, heat_capacity" << std::endl;

    for (int i=0; i<measurement["beta"].size(); i++){
        myfile << measurement["beta"][i] << ", " 
            << measurement["magnetization"][i] << ", " 
            << measurement["m2"][i] << ", " 
            << measurement["energy"][i] << ", " 
            << measurement["heat_capacity"][i] << ", "
            << std::endl;
    }
    myfile.close();
}

