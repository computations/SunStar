//nj.h
//Ben Bettisworth
#pragma once

#include "tree.h"
#include <vector>
#include <string>

class nj_t{
    public:
        nj_t(const std::vector<double> &, const std::vector<std::string>&);
        tree_t get_tree();
        
    private:
        void compute_q();
        void compute_r();
        void recompute_dists();

        void find_pair();
        void join_pair();
        void join_final();
        void join_final_small();
        void make_tree();
        //void flatten_tree();

        void clean_up();

        //these vectors need to stay in "lock step"
        //each entry needs to corrispond to a specific entry in another
        //these vectors are "flat"
        std::vector<node_t*> _tree;
        std::vector<node_t*> _unroot;
        std::vector<double> _r_vec;

        //these vectors are "square"
        std::vector<double> _dists;
        std::vector<double> _q_vec;

        size_t _row_size;
        size_t _i, _j;

        size_t _tree_size;

        tree_t _final_tree;
};
