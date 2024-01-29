#ifndef monteCarlo_h
#define monteCarlo_h


#include "isingState.h"
#include <map>
#include <random>

namespace ising {

    class IsingMonteCarlo
    {
    public:
        IsingMonteCarlo(IsingState ising_state);
        virtual void update()=0;
        void thermalize(int n_iters=10000);
        std::map<std::string, std::vector<double>> sample(int n_samples=500, int n_iters_per_sample=20, int n_thermalization_iters=10000);


        // below are member attributes
        IsingState ising_state;
        std::default_random_engine generator;
        std::uniform_int_distribution<int> int_distribution;
        std::uniform_real_distribution<double> real_distribution;

    };

    class IsingMetropolis : public IsingMonteCarlo
    { 
    public:
        IsingMetropolis(IsingState ising_state);

        double delta_energy(int flip_index);
        void update();

    };

    class IsingWolff: public IsingMonteCarlo
    {
    public:
        IsingWolff(IsingState ising_state);
        void update();
        void flip_wolff(int pos);
    private:
        double p_add;
    };


}

#endif

