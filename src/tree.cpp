#include "debug.h"
#include "tree.h"
#include "newick.h"
#include <string>
using std::string;
#include <unordered_map>
using std::unordered_map;
#include<vector>
using std::vector;
#include <stack>
using std::stack;
#include <sstream>
using std::ostringstream;
#include <functional>
using std::function;
#include <cassert>

void node_t::set_weights(function<float(size_t)> w_func, size_t depth){
    _weight = w_func(depth);
    if(_children){
        _lchild->set_weights(w_func, depth+1);
        _rchild->set_weights(w_func, depth+1);
    }
}

tree_t::tree_t(const tree_t& t){
    _size = t._size;
    _tree = new node_t[_size];
    //need to update lineage
    //i feel dirty
    long int offset = _tree - t._tree;
    for(size_t i =0;i<_size;++i){
        _tree[i] = t._tree[i];
    }
    for(size_t i=0;i<_size;++i){
        _tree[i]._parent += offset;
        _tree[i]._lchild += offset;
        _tree[i]._rchild += offset;
    }
    _unroot = t._unroot;
    for(size_t i=0;i<t._unroot.size();++i){
        _unroot[i] += offset;
    }
}

tree_t::tree_t(node_t* tree, size_t size, const vector<node_t*>& unroots){
    
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
    //resize the tree array if we need to
    if(_tree && _size != t._size){
        _size = t._size;
        delete[] _tree;
        _tree = new node_t[_size];
    }
    
    for(size_t i=0;i<_size;++i){
        _tree[i] = t._tree[i];
    }
    long int offset = _tree - t._tree;
    for(size_t i =0;i<_size;++i){
        _tree[i] = t._tree[i];
    }
    for(size_t i=0;i<_size;++i){
        _tree[i]._parent += offset;
        _tree[i]._lchild += offset;
        _tree[i]._rchild += offset;
    }
    _unroot = t._unroot;
    for(size_t i=0;i<t._unroot.size();++i){
        _unroot[i] += offset;
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
            for(size_t j=0;j<_size;++j){
                if(!_tree[j]._children){
                    size_t dest_matrix_index = label_map.at(_tree[j]._label);
                    debug_print("calculating distance for (%lu,%lu), putting in: (%lu,%lu)", 
                            i, j, matrix_index,dest_matrix_index);
                    dists[row_size*matrix_index+dest_matrix_index] = calc_distance(_tree+i,_tree+j);
                }
            }
        }
    }
    //need to reflect the matrix
    //because i is unsigned, when we go negative, we actually go up to 2^64-1.
    //fortunatly, that is always going to be larger than the size of the array
    //so, we check to see if i<_size, just like a regular for loop
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
float tree_t::calc_distance(node_t* src, node_t* dst){
    debug_print("calculating distance between (%p, %p)", src, dst);
    if(src==dst){
        debug_string("src and dst are the same, returning zero");
        return 0.0;
    }
    auto src_list = get_parents_of(src);
    auto dst_list = get_parents_of(dst);

    size_t src_index = src_list.size()-1;
    size_t dst_index = dst_list.size()-1;


    //walk through the list until the lists diverge
    while(true){
        assert_string((src_index < src_list.size() || dst_index < dst_list.size())
                , "Parent lists don't converge, but not same index");

        if(src_list[src_index] != dst_list[dst_index]){
            src_index++; dst_index++;//need to back up a step
            break;
        }
        src_index--; dst_index--;
    }

    //calculating common parent can be faster, since I have a list of parents already, but this is fine
    float ret =  parent_distance(src, src_list[src_index]) + parent_distance(dst,dst_list[dst_index]);
    debug_print("returning the distance %f", ret);
    return ret;
}

vector<node_t*> tree_t::get_parents_of(node_t* cur_node){
    debug_print("getting the parents of %s", cur_node->_label.c_str());
    vector<node_t*> parent_list;
    parent_list.push_back(cur_node);

    while(cur_node->_parent != cur_node && cur_node->_parent!=0){
        debug_print("current node: %p", cur_node);
        parent_list.push_back(cur_node->_parent);
        cur_node = cur_node->_parent;
    }
    return parent_list;
}

float tree_t::parent_distance(node_t* child, node_t* parent){
    float distance = 0;
    while(child!=parent){
        distance+=child->_weight;
        child = child->_parent;
    }
    return distance;
}

string node_t::to_string(node_t* root){
    ostringstream ret;
    if(_children){
        ret<<"("<<_lchild->to_string(root)
            <<","<<_rchild->to_string(root)<<")"
            <<_label<<":"<<_weight;
    }
    else{
        ret<<_label<<":"<<_weight;
    }
    return ret.str();
}

string tree_t::to_string(){
    ostringstream ret;

    if(_unroot.size()>1)
        ret<<"(";

    for(size_t i=0;i<_unroot.size();++i){
        ret<<_unroot[i]->to_string(_tree);
        if(i!=_unroot.size()-1)
            ret<<",";
    }
    if(_unroot.size()>1)
        ret<<")";

    return ret.str();
}

string tree_t::print_labels(){
    ostringstream ret;
    for(size_t i=0;i<_size;++i){
        ret<<_tree[i]._label<< "(" << _tree[i]._parent <<"," <<_tree[i]._lchild<<","<<_tree[i]._rchild<<")";
        if(i!=_size-1) ret<<" | ";
    }
    return ret.str();
}

void tree_t::set_weights(function<float(size_t)> w_func){
    for(auto n : _unroot){
        n->set_weights(w_func, 0);
    }
}

void tree_t::set_weights(const vector<float>& w_vec){
    set_weights([&w_vec](size_t d){
            assert_string(d < w_vec.size(), "out of bounds for passed float vector");
            if(d == 0) return 0.0f;
            return w_vec[d-1];
    });
}
