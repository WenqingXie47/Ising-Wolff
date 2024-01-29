#include "gridState.h"
#include <cmath>
#include <random>
#include <iostream>


namespace ising{


    GridState::GridState(int dim, int length)
        : dim{dim}, length{length}
    {
        this->n_sites = pow(length,dim);
        this->n_neighbours = pow(2,dim);
        this->initGrid();
        this->initNeighbours();
    }

    void GridState::initGrid()
    {
        
        std::default_random_engine generator;
        std::bernoulli_distribution distribution(0.5);
        this->state = std::vector<int> (this->n_sites);
        for (int i=0; i<(this->n_sites); i++){
            if (distribution(generator)==true){
                this->state[i] = 1;
            }
            else {
                this->state[i] = -1;
            }
        } 
    }

    void GridState::initNeighbours()
    {
        this->neighbours = std::vector<std::vector<int>> (this->n_sites);
        
        for (int i=0; i<(this->n_sites); i++){
            this->neighbours[i] = std::vector<int> (this->n_neighbours);
            for (int d=0; d<(this->dim); d++){
                std::vector<int> grid_index = this->getGridIndex(i);
                std::vector<int> neighbour1_index =  grid_index;
                neighbour1_index[d] = ((grid_index[d]-1) + (this->length)) %(this->length);
                std::vector<int> neighbour2_index =  grid_index;
                neighbour2_index[d] = (grid_index[d]+1)%(this->length);
                this->neighbours[i][2*d] = this->getFlattenIndex(neighbour1_index);
                this->neighbours[i][2*d+1] = this->getFlattenIndex(neighbour2_index);
            }
        }
    }

    int GridState::getFlattenIndex(std::vector<int> grid_index) const
    {
        int index = 0;
        for (int d=0; d<(this->dim); d++){
            index += pow(this->length, d) * grid_index[d];
        }
        return index;
    }

    std::vector<int> GridState::getGridIndex(int flatten_index) const
    {
        std::vector<int> index = std::vector<int> (this->dim);
        for (int d=0; d<this->dim; d++) {
            index[d] = flatten_index % this->length;
            flatten_index = flatten_index / this->length;
        }
        return index;
    }

    void GridState::get_adjacent_matrix() const 
    {
        return;
    }


    // Below are getter and setter functions that no longer used

    // const std::vector<int>& GridState::getStateVector()
    // {
    //     return this->state;
    // }

    // void GridState::setStateVector(std::vector<int>& new_state)
    // {
    //     this->state = new_state;
    // }



} // end namespace