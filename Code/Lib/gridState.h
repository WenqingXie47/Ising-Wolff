#ifndef gridState_h
#define gridState_h


#include <vector>



namespace ising
{

    class GridState 
    {
    public:
        // below are functions
        GridState(int dim, int length);
        void initGrid();
        void initNeighbours();
        int getFlattenIndex(std::vector<int> grid_index) const;
        std::vector<int> getGridIndex(int flatten_index) const;
        void get_adjacent_matrix() const;

        
        // below are attributes
        int dim;
        int length;
        int n_sites;
        int n_neighbours;
        std::vector<int> state;
        std::vector<std::vector<int>> neighbours;
    }; // end of  class GridState 

} // end namespace

#endif