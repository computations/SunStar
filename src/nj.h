#pragma once

#include "tree.h"
#include <vector>
#include <string>

class nj_t{
    public:
        nj_t(const std::vector<float> &, const std::vector<std::string>&);
        tree_t neighbor_join(std::vector<float>);

        
    private:


        void compute_q();
        void compute_r();
        void recompute_dists();

        void find_pair();
        void join_pair();

        //these vectors need to stay in "lock step"
        //each entry needs to corrispond to a specific entry in another
        //these vectors are "flat"
        std::vector<node_t*> _tree;
        std::vector<float> _r_vec;

        //these vectors are "square"
        std::vector<float> _dists;
        std::vector<float> _q_vec;

        size_t _row_size;
        size_t _i, _j;
};
