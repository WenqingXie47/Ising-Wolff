#ifndef isingState_h
#define isingState_h


#include <vector>
#include "gridState.h"


namespace ising
{

    class IsingState : public GridState 
    {
    public:
        IsingState(int dim, int length, double beta=1, double h=0, double J=1);
        double getMagnetization();
        double getEnergy();
        int sum_neighbour(int index);
            
        double beta;
        double h;
        double J;
    }; // end of class IsingState

} // end namespace

#endif
