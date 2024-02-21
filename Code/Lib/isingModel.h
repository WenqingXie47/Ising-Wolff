#ifndef isingState_h
#define isingState_h


#include <vector>
#include "gridState.h"


namespace ising
{
    class SpinModel
    {
    public:
        SpinModel(GridState& grid, double beta, double J=1, double h=0);
        std::vector<std::vector<int>>& get_neighbours();
        int sum_neighbour(int index);

        

        double beta;
        double h;
        double J;
        GridState grid;
        std::vector<int>& state;
    };


    class IsingModel : public SpinModel
    {
    public:
        IsingModel(GridState& grid, double beta, double J=1, double h=0);
        double getMagnetization();
        double getEnergy();
            
    }; // end of class IsingState

} // end namespace

#endif
