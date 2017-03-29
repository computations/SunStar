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

#include <functional>
using std::function;

inline void init_array(double* v, size_t s){
    for(size_t i =0;i<s*s;++i){
        v[i]= 0.0;
    }
}

star_t::star_t(const vector<string>& newick_trees){
    _tree_collection.reserve(newick_trees.size());
    for(auto &&s : newick_trees){
        _tree_collection.emplace_back(s);
    }
    _label_map = _tree_collection.front().make_label_map();
}

void star_t::calc_average_distances(){
    debug_print_map("label map", _label_map);
    debug_print("front tree: %s", _tree_collection.front().to_string().c_str());
    size_t row_size = _label_map.size();
    debug_print("row_size: %lu", row_size);
    double* dists = new double[row_size*row_size];
    _avg_dists.resize(row_size*row_size, 0);
    init_array(dists, row_size);

    for(size_t i=0;i<_tree_collection.size();++i){
        init_array(dists, row_size);
        debug_print("current tree with set weights: %s", _tree_collection[i].to_string(5).c_str());
        _tree_collection[i].calc_distance_matrix(_label_map, dists);
        debug_print("current tree: %s", _tree_collection[i].print_labels().c_str());
        for(size_t j=0;j<row_size; ++j){
            for(size_t k = 0; k<row_size; ++k){
                _avg_dists[j*row_size+k]+=dists[j*row_size+k];
            }
        }
    }
    debug_matrix("dists before average", _avg_dists, row_size);
    for(size_t i=0; i<row_size*row_size;++i){
        _avg_dists[i]/=(double)_tree_collection.size();
    }
    debug_matrix("dists after average", _avg_dists, row_size);
    delete[] dists;
}

vector<string> invert_label_map(unordered_map<string, size_t> lm){
    debug_print_map("label map before inversion", lm);
    vector<string> ret;
    ret.resize(lm.size());
    for(auto &&kv_pair : lm){
        ret[kv_pair.second] = kv_pair.first;
    }
    return ret;
}

tree_t star_t::get_tree(){
    calc_average_distances();
    return nj(_avg_dists, invert_label_map(_label_map));
}

tree_t star_t::get_tree(const function<double(size_t)>& f){
    double max = 0.0;
    size_t depth = get_size();
    for(size_t i =0;i<depth;++i){
        max+= f(i);
    }
    for(auto& t:_tree_collection){
        t.set_weights(f, max);
    }
    calc_average_distances();
    return nj(_avg_dists, invert_label_map(_label_map));
}

tree_t star_t::get_tree(const vector<double>& v){
    double max = 0.0;
    size_t depth = get_size();
    for(size_t i =0;i<depth;++i){
        max+= v[i];
    }
    for(auto& t:_tree_collection){
        t.set_weights(v, max);
        debug_string(t.to_string().c_str());
    }
    calc_average_distances();
    return nj(_avg_dists, invert_label_map(_label_map));
}

size_t star_t::get_size(){
    size_t max = 0;
    for(const auto& t:_tree_collection){
       size_t tmp = t.get_depth();
       if(max < tmp){
           max = tmp;
       }
    }
    return max;
}

void star_t::set_outgroup(const string& outgroup){
    for(auto& t:_tree_collection){
        t.set_outgroup(outgroup);
    }
}

string star_t::get_first_label(){
    string ret;
    for(auto kv : _label_map){
        ret = kv.first;
        if(!ret.empty() && ret != " ") break;
    }
    return ret;
}
