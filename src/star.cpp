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

std::vector<float> calc_average_distances(std::vector<tree_t> tree_vector){
    auto label_map =  tree_vector.front().make_label_map();
    debug_print("front tree: %s", tree_vector.front().to_string().c_str());
    size_t row_size = label_map.size();
    debug_print("row_size: %lu", row_size);
    float* dists = new float[row_size*row_size];
    std::vector<float> avg_dists;
    avg_dists.resize(row_size*row_size, 0);
    init_array(dists, row_size);

    for(size_t i=0;i<tree_vector.size();++i){
        init_array(dists, row_size);
        tree_vector[i].calc_distance_matrix(label_map, dists);
        debug_print("current tree: %s", tree_vector[i].print_labels().c_str());
        for(size_t j=0;j<row_size; ++j){
            for(size_t k = 0; k<row_size; ++k){
                avg_dists[j*row_size+k]+=dists[j*row_size+k];
                std::cerr<<dists[j*row_size+k]<<" ";
            }
            std::cerr<<"\n";
        }
    }
    for(size_t j=0;j<row_size*row_size; ++j){
        debug_print("avg_dists no average: %f", avg_dists[j]);
    }
    for(size_t i=0; i<row_size*row_size;++i){
        avg_dists[i]/=(float)tree_vector.size();
    }
    return avg_dists;
}
