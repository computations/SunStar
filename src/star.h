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
#include <functional>

class star_t{
    public:
        star_t(const std::vector<std::string>&);
        tree_t get_tree();
        tree_t get_tree(const std::function<double(size_t)>&);
        tree_t get_tree(const std::vector<double>&);
        size_t get_size();
        void set_outgroup(const std::string&);
        std::string get_first_label();
    private:
        void calc_average_distances();

        std::vector<double> _avg_dists;
        std::vector<tree_t> _tree_collection;
        std::unordered_map<std::string, size_t> _label_map;
};
