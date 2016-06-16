#include "debug.h"
#include "tree.h"
#include "newick.h"
#include <string>
using std::string;
#include <unordered_map>
using std::unordered_map;
#include<vector>
using std::vector;
#include <sstream>
using std::ostringstream;
#include <cassert>

tree_t::tree_t(const tree_t& t){
    _size = t._size;
    _tree = new node_t[_size];
}

//takes a tree specified by newick notation, and makes a tree_t
tree_t::tree_t(const string& newick){
    _tree = make_tree_from_newick(newick, _size);
}

tree_t::~tree_t(){
    if(_tree)
        delete[] _tree;
}

tree_t& tree_t::operator=(const tree_t& t){
    if(_tree)
        delete[] _tree;
    
    _size = t._size;
    _tree = new node_t[_size];
    for(size_t i=0;i<_size;++i){
        _tree[i] = t._tree[i];
    }
    return *this;
}

//we use a label to index map to make the matrix well ordered
//this is so we can do a blind average later on, and not have to worry about the ordering of the array
float* tree_t::calc_distance_matrix(const std::unordered_map<string, size_t>& label_map){
    size_t row_size = label_map.size();
    float* dists = new float[row_size*row_size];
    calc_distance_matrix(label_map, dists);
    return dists;
}

void tree_t::calc_distance_matrix(const std::unordered_map<string, size_t>& label_map, float* dists){
    size_t row_size = label_map.size();
    for(size_t i=0;i<_size;++i){
        //if the node has no children, then it is a leaf, and we need to find distances
        if(!_tree[i]._children){
            size_t matrix_index = label_map.at(_tree[i]._label);
            dists[matrix_index*row_size] = 0;
            for(size_t j=i+1;j<_size;++j){
                if(!_tree[j]._children){
                    size_t dest_matrix_index = label_map.at(_tree[j]._label);
                    dists[row_size*matrix_index+dest_matrix_index] = calc_distance(i,j);
                }
            }
        }
    }
    //need to reflect the matrix
    //because i is unsigned, when we go negative, we actually go up to 2^64-1.
    //fortunatly, that is always going to be larger than the size of the array
    //so, we check to see if i<_size, just like a regular for loop
    for(size_t i=_size-1; i<_size; --i){
        for(size_t j=i; j<_size; --j){
            dists[row_size*i+j] = dists[row_size*j+i];
        }
    }
}


//need to make a map of labels to indices, but the order doesnt really matter
//so, this is inteded to be called for on the first tree, and never again
//the return of this function is meant to be fed into the function
//  calc_distance_matrix();
//so that it can calculate ta specific matrix that is ''well ordered''
std::unordered_map<string, size_t> tree_t::make_label_map(){
    size_t label_index = 0;
    std::unordered_map<string, size_t> label_map;
    for(size_t i = 0; i<_size; ++i){
        if(!_tree[i]._children){
            label_map[_tree[i]._label] = label_index++;
        }
    }
    return label_map;
}

//calculate the distance between two nodes
//game plan:
//  make a list of parents for each node
//  compare those lists from the back (ie, root first)
//  when those lists diverge, thats the common parent
float tree_t::calc_distance(size_t src, size_t dst){

    if(src==dst) return 0.0;
    auto src_list = get_parents_of(src);
    auto dst_list = get_parents_of(dst);

    assert_string(src_list.back() == dst_list.back(), "leaves don't share a root");

    size_t src_index = src_list.size()-1;
    size_t dst_index = dst_list.size()-1;
    debug_print("src_list size: %lu | dst_list size: %lu", src_list.size(), dst_list.size());


    //walk through the list until the lists diverge
    while(true){
        debug_print("src_index: %lu | dst_index : %lu" , src_index, dst_index);
        debug_print("current src node: %lu | current dst node: %lu" , src_list[src_index], dst_list[dst_index]);
        assert_string((src_index < src_list.size() || dst_index < dst_list.size())
                , "Parent lists don't converge, but not same index");

        if(src_list[src_index] != dst_list[dst_index]){
            src_index++; dst_index++;//need to back up a step
            break;
        }
        src_index--; dst_index--;
    }

    assert_string((src_list[src_index] == dst_list[dst_index]), "parents don't equal each other");
    size_t common_parent = src_list[src_index];

    //calculating common parent can be faster, since I have a list of parents already, but this is fine
    return parent_distance(src, common_parent) + parent_distance(dst,common_parent);
}

vector<size_t> tree_t::get_parents_of(size_t index){
    vector<size_t> parent_list;
    parent_list.push_back(index);

    while(_tree[index]._parent != index){
        parent_list.push_back(_tree[index]._parent);
        index = _tree[index]._parent;
    }
    return parent_list;
}

float tree_t::parent_distance(size_t child, size_t parent){
    float distance = 0;
    while(child!=parent){
        distance+=_tree[child]._weight;
        child = _tree[child]._parent;
    }
    return distance;
}

string node_t::to_string(node_t* root){
    ostringstream ret;
    if(_children){
        ret<<"("<<(root+_lchild)->to_string(root)
            <<","<<(root+_rchild)->to_string(root)<<")"
            <<_label<<":"<<_weight;
    }
    else{
        ret<<_label<<":"<<_weight;
    }
    return ret.str();
}

string tree_t::to_string(){
    return _tree->to_string(_tree);
}

string tree_t::print_labels(){
    ostringstream ret;
    for(size_t i=0;i<_size;++i){
        ret<<_tree[i]._label<< "(" <<_tree[i]._lchild<<","<<_tree[i]._rchild<<")";
        if(i!=_size-1) ret<<" | ";
    }
    return ret.str();
}
