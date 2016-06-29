#include "star.h"
#include "tree.h"
#include "debug.h"
#include <vector>
using std::vector;
#include <string>
using std::string;
#include <iostream>

inline void init_array(float* v, size_t s){
    for(size_t i =0;i<s*s;++i){
        v[i]= 0.0;
    }
}

star_t::star_t(const vector<string>& newick_trees){
    vector<tree_t> trees;
    trees.reserve(newick_trees.size());
    for(auto &&s : newick_trees){
        trees.emplace_back(s);
    }
    calc_average_distances(trees);
}

void star_t::calc_average_distances(std::vector<tree_t>& tree_vector){
    auto label_map =  tree_vector.front().make_label_map();
    debug_print("front tree: %s", tree_vector.front().to_string().c_str());
    size_t row_size = label_map.size();
    debug_print("row_size: %lu", row_size);
    float* dists = new float[row_size*row_size];
    _avg_dists.resize(row_size*row_size, 0);
    init_array(dists, row_size);

    for(size_t i=0;i<tree_vector.size();++i){
        init_array(dists, row_size);
        tree_vector[i].calc_distance_matrix(label_map, dists);
        debug_print("current tree: %s", tree_vector[i].print_labels().c_str());
        for(size_t j=0;j<row_size; ++j){
            for(size_t k = 0; k<row_size; ++k){
                _avg_dists[j*row_size+k]+=dists[j*row_size+k];
                std::cerr<<dists[j*row_size+k]<<" ";
            }
            std::cerr<<"\n";
        }
    }
    for(size_t j=0;j<row_size*row_size; ++j){
        debug_print("_avg_dists no average: %f", _avg_dists[j]);
    }
    for(size_t i=0; i<row_size*row_size;++i){
        _avg_dists[i]/=(float)tree_vector.size();
    }
}
