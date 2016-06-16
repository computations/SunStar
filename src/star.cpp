#include "star.h"
#include "tree.h"
#include <vector>
using std::vector;
#include <string>
using std::string;

std::vector<float> calc_average_distances(std::vector<tree_t> tree_vector){
    auto label_map =  tree_vector.front().make_label_map();
    size_t row_size = label_map.size();
    float* dists = new float[row_size*row_size];
    std::vector<float> avg_dists;
    avg_dists.resize(row_size*row_size, 0);
    for(size_t i=0;i<tree_vector.size();++i){
        tree_vector[i].calc_distance_matrix(label_map, dists);
        for(size_t j=0;j<row_size*row_size; ++j){
            avg_dists[j]+=dists[j];
        }
    }
    for(size_t i=0; i<row_size*row_size;++i){
        avg_dists[i]/=(float)tree_vector.size();
    }
    return avg_dists;
}
