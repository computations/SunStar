//star.h
//Ben Bettisworth
//2016-06-29
//a class to implement the state stuff for the STAR algorithm
//for more details about this implementation STAR algorithm, please see:
//  http://www.dms.uaf.edu/~jrhodes/papers/STARandGeneralizations.pdf
#pragma once
#include "tree.h"
#include <vector>
#include <string>
#include <unordered_map>

class star_t{
    public:
        star_t(const std::vector<std::string>&);
        tree_t get_tree();
    private:
        void calc_average_distances(std::vector<tree_t>&);
        void make_tree();

        tree_t _final_tree;
        std::vector<double> _avg_dists;
        std::unordered_map<std::string, size_t> _label_map;
};
