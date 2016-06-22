#pragma once

#include "tree.h"
#include <vector>

class nj_t{
    public:
        nj_t(std::vector<float>);
        tree_t neighbor_join(std::vector<float>);

        void compute_q();
        void compute_r();

        void find_pair();
        void join_pair();
        
    private:
        std::vector<float> _dists;
        std::vector<float> _r_vec;
        std::vector<float> _q_vec;
        tree_t _tree;
        size_t _row_size;
        size_t _i, _j;
}
