#ifndef monteCarlo_h
#define monteCarlo_h


#include "isingModel.h"
#include <map>
#include <random>

namespace ising {

    class IsingMonteCarlo
    {
    public:
        IsingMonteCarlo(IsingModel ising_model);
        virtual void update()=0;
        void thermalize(int n_iters=10000);
        std::map<std::string, std::vector<double>> sample(int n_samples=500, int n_iters_per_sample=20, int n_thermalization_iters=10000);


        // below are member attributes
        IsingModel ising_model;
        std::default_random_engine generator;
        std::uniform_int_distribution<int> int_distribution;
        std::uniform_real_distribution<double> real_distribution;

    };

    class IsingMetropolis : public IsingMonteCarlo
    { 
    public:
        IsingMetropolis(IsingModel ising_model);

        double delta_energy(int flip_index);
        void update();

    };

    class IsingWolff: public IsingMonteCarlo
    {
    public:
        IsingWolff(IsingModel ising_model);
        void update();
        void flip_wolff(int pos);
    private:
        double p_add;
    };


}

#endif

