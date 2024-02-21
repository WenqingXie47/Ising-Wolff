#include "monteCarlo.h"

#include <string>
#include <cmath>
#include <iostream>


namespace ising{

    IsingMonteCarlo::IsingMonteCarlo(IsingModel ising_model)
        :ising_model{ising_model}
    {
        generator = std::default_random_engine();
        int_distribution = std::uniform_int_distribution<int> (0,ising_model.grid.n_sites-1);
        real_distribution = std::uniform_real_distribution<double> (0.0,1.0);
    }


    void IsingMonteCarlo::thermalize(int n_iters)
    {
        this->ising_model.grid.initGrid();
        for (int i=0; i<n_iters; i++){
            this->update();
        }
    }

    std::map<std::string, std::vector<double>> IsingMonteCarlo::sample(
        int n_samples, int n_iters_per_sample, int n_thermalization_iters
    )
    {
        this->thermalize(n_thermalization_iters);

        std::map<std::string, std::vector<double>> history{
            {"energy", std::vector<double>{}},
            {"magnetization", std::vector<double>{}},
        };

        for (int i=0; i<n_samples; i++){
            history["energy"].push_back(this->ising_model.getEnergy());
            history["magnetization"].push_back(this->ising_model.getMagnetization());
            
            for (int j=0; j<n_iters_per_sample; j++){
                this->update();
            }
        }
        return history;
    }


    IsingMetropolis::IsingMetropolis(IsingModel ising_model)
        :IsingMonteCarlo(ising_model)
    {}

    void IsingMetropolis::update()
    {
        // randomly choose an index
        int index = this->int_distribution(this->generator);

        // delta energy = new energy - old energy
        double delta_energy = this->delta_energy(index);

        double random_real = this->real_distribution(this->generator);
        if (random_real < std::exp(-delta_energy*ising_model.beta)){
            this->ising_model.state[index] *= -1; 
        }

    }

    double IsingMetropolis::delta_energy(int flip_index)
    {
        int spin = this->ising_model.state[flip_index];

        // check the neighbours of the chosen spin
        int neighbour_sum = this->ising_model.sum_neighbour(flip_index);


        // delta energy = new energy - old energy
        double delta_energy = 2* (double) spin * (ising_model.h + (double) neighbour_sum * ising_model.J);
        return delta_energy;
    }


    IsingWolff::IsingWolff(IsingModel ising_model)
        :IsingMonteCarlo(ising_model)
    {
        this->p_add = 1- std::exp(-2*ising_model.beta*ising_model.J);
    }


    void IsingWolff::update()
    {
        // randomly choose an index
        int index = this->int_distribution(this->generator);
        flip_wolff(index);

    }


    void IsingWolff::flip_wolff(int index)
    {
        this->ising_model.state[index] *= -1;
        for (int neighbour_index: this->ising_model.get_neighbours()[index]){
            if (ising_model.state[neighbour_index]!= ising_model.state[index] &&
                this->real_distribution(this->generator) < this->p_add){
                flip_wolff(neighbour_index);
            }
        }
    }
} // end namespace 


