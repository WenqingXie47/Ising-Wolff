#include "isingState.h"
#include <iostream>

namespace ising
{

    IsingState::IsingState(int dim, int length, double beta, double h, double J)
        : beta{beta}, h{h}, J{J}, GridState(dim,length)
    {}

    int IsingState::sum_neighbour(int index)
    {
        int neighbour_sum = 0; 
        for (int neighbour_id : (this->neighbours)[index]){
            neighbour_sum += this->state[neighbour_id];
        }
        return neighbour_sum;
    }

    double IsingState::getMagnetization()
    {
        double m = 0;
        for (int spin: this->state){
            m += (double) spin;
        }
        m /= (this->n_sites);
        return m;
    }

    double IsingState::getEnergy()
    {
        double energy_spin = 0;
        
        for (int spin_id=0; spin_id< (this->n_sites); spin_id++){
            double spin = this->state[spin_id];
            double neighbour_sum = this->sum_neighbour(spin_id);
            energy_spin += (spin*neighbour_sum/2);
        }
        energy_spin *= -(this->J);
        double energy_field = getMagnetization() * -(this->h);

        return energy_spin + energy_field; 
    }

} // end namespace