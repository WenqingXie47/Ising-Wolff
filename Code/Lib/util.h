#ifndef util_h
#define util_h


#include <vector>

namespace ising{

    double mean(std::vector<double> list);
    double mean_abs(std::vector<double> list);
    double mean_squared(std::vector<double> list);
    double variance(std::vector<double> list);

} // end namespace

#endif