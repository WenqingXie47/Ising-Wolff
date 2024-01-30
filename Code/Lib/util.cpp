#include "util.h"

#include <cmath>


namespace ising{

    double mean(std::vector<double> list)
    {
        double sum = 0;
        for (double n : list){
            sum += n;
        }
        return sum/list.size();
    }

    double mean_abs(std::vector<double> list)
    {
        double sum = 0;
        for (double n : list){
            sum += std::abs(n);
        }
        return sum/list.size();
    }

    double mean_squared(std::vector<double> list)
    {
        double sum = 0;
        for (double n : list){
            sum += n*n;
        }
        return sum/list.size();
    }

    double variance(std::vector<double> list)
    {
        double m_squared = mean_squared(list);
        double m = mean(list);
        return m_squared - m*m;
    }


} // end namespace