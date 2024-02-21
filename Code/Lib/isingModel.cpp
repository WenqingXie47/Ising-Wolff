#include "isingModel.h"
#include <iostream>

namespace ising
{

    SpinModel::SpinModel(GridState& grid, double beta, double J, double h)
        : grid{grid}, state{grid.state}, beta{beta}, h{h}, J{J}
    {
        this->grid.initGrid();
    }



    std::vector<std::vector<int>>& SpinModel::get_neighbours()
    {
        return this->grid.neighbours;
    }

    int SpinModel::sum_neighbour(int index)
    {
        int neighbour_sum = 0; 
        for (int neighbour_id : (this->get_neighbours())[index]){
            neighbour_sum += this->state[neighbour_id];
        }
        return neighbour_sum;
    }

    IsingModel::IsingModel(GridState& grid, double beta, double J, double h)
        : SpinModel(grid,beta,J,h)
    {}

    double IsingModel::getMagnetization()
    {
        double m = 0;
        for (int spin: this->state){
            m += (double) spin;
        }
        m /= (this->grid.n_sites);
        return m;
    }

    double IsingModel::getEnergy()
    {
        double energy_spin = 0;

        for (int spin_id=0; spin_id< (this->grid.n_sites); spin_id++){
            double spin = this->state[spin_id];
            double neighbour_sum = this->sum_neighbour(spin_id);
            energy_spin += (spin*neighbour_sum/2);
        }
        energy_spin *= -(this->J);
        double m = getMagnetization();
        m *= this->grid.n_sites;
        double energy_field = m * -(this->h);
        double energy = energy_spin + energy_field;
        return energy; 
    }

} // end namespace