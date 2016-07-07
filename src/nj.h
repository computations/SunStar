#pragma once

#include "tree.h"
#include <vector>
#include <string>

class nj_t{
    public:
        nj_t(const std::vector<float> &, const std::vector<std::string>&);
        ~nj_t();
        tree_t get_tree();
        
    private:
        void compute_q();
        void compute_r();
        void recompute_dists();

        void find_pair();
        void join_pair();
        void join_final();
        void make_tree();
        void flatten_tree();

        void clean_up();

        //these vectors need to stay in "lock step"
        //each entry needs to corrispond to a specific entry in another
        //these vectors are "flat"
        std::vector<node_t*> _tree;
        std::vector<node_t*> _unroot;
        std::vector<float> _r_vec;

        //these vectors are "square"
        std::vector<float> _dists;
        std::vector<float> _q_vec;

        size_t _row_size;
        size_t _i, _j;

        size_t _tree_size;
        node_t* _flat_tree;

        tree_t _final_tree;
};
