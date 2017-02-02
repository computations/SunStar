//star.cpp
//Ben Bettisworth
//Implementation of the STAR algorithm.
//STAR is the algorithm where we compute a distance matrix of the various trees
//and then average the matrices entrywise. Then use NJ to make an ''average''
//tree.

#include "star.h"
#include "tree.h"
#include "nj.h"
#include "debug.h"

#include <vector>
using std::vector;

#include <string>
using std::string;

#include <unordered_map>
using std::unordered_map;

#include <iostream>

inline void init_array(double* v, size_t s){
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
    make_tree();
}

void star_t::calc_average_distances(std::vector<tree_t>& tree_vector){
    _label_map =  tree_vector.front().make_label_map();
    debug_print("front tree: %s", tree_vector.front().to_string().c_str());
    size_t row_size = _label_map.size();
    debug_print("row_size: %lu", row_size);
    double* dists = new double[row_size*row_size];
    _avg_dists.resize(row_size*row_size, 0);
    init_array(dists, row_size);

    for(size_t i=0;i<tree_vector.size();++i){
        init_array(dists, row_size);
        tree_vector[i].set_weights(1.0);
        debug_print("current tree with set weights: %s", tree_vector[i].to_string().c_str());
        tree_vector[i].calc_distance_matrix(_label_map, dists);
        debug_print("current tree: %s", tree_vector[i].print_labels().c_str());
        for(size_t j=0;j<row_size; ++j){
            for(size_t k = 0; k<row_size; ++k){
                _avg_dists[j*row_size+k]+=dists[j*row_size+k];
            }
        }
    }
    for(size_t j=0;j<row_size*row_size; ++j){
        debug_print("_avg_dists no average: %f", _avg_dists[j]);
    }
    for(size_t i=0; i<row_size*row_size;++i){
        _avg_dists[i]/=(double)tree_vector.size();
    }
    delete[] dists;
}

vector<string> invert_label_map(unordered_map<string, size_t> lm){
    vector<string> ret;
    ret.resize(lm.size());
    for(auto &&kv_pair : lm){
        ret[kv_pair.second] = kv_pair.first;
    }

    return ret;
}

void star_t::make_tree(){
    nj_t nj(_avg_dists, invert_label_map(_label_map));
    _final_tree = nj.get_tree();
}

tree_t star_t::get_tree(){
    return _final_tree;
}
