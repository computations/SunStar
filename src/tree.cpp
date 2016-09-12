//tree.cpp
//Ben Bettisworth
//A class that implements a phylogentic tree. Features of the class:
//  - Handles the newick parsing (via newick.cpp)
//  - Packed to take advantage of cache locallity
//  - Handles labels

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

#include <queue>
using std::queue;

#include <sstream>
using std::ostringstream;

#include <functional>
using std::function;

#include <cassert>
#include <iostream>

size_t node_t::count_nodes(){
    size_t children = 0;
    if(_children){
        children += _lchild->count_nodes();
        children += _rchild->count_nodes();
    }
    return children+1;
}

void node_t::update_children(const unordered_map<node_t*, node_t*> node_map){
    if(_parent)
        _parent = node_map.at(_parent);
    if(_children){
        _lchild = node_map.at(_lchild);
        _rchild = node_map.at(_rchild);
        _lchild->update_children(node_map);
        _rchild->update_children(node_map);
    }
}

//TODO: fix this function
//currently sets the weights according to the function, but it needs to 
//ensure that the tree is ultrametric
void node_t::set_weights(function<float(size_t)> w_func, size_t depth){
    _weight = w_func(depth);
    if(_children){
        _lchild->set_weights(w_func, depth+1);
        _rchild->set_weights(w_func, depth+1);
    }
}

string node_t::sort(){
    if(_children){
        auto lchild_string = _lchild->sort();
        auto rchild_string = _rchild->sort();
        if(rchild_string < lchild_string){
            std::swap(_lchild, _rchild);
            std::swap(lchild_string, rchild_string);
        }
        if(lchild_string < _label || _label.empty()){
            return lchild_string;
        }
    }
    return _label;
}

tree_t::tree_t(const tree_t& t){
    make_flat_tree(t._unroot);
}

//Traverses the tree, and compresses it into an array
void tree_t::make_flat_tree(const vector<node_t*>& unroot){
    unordered_map<node_t*, node_t*> node_map;
    stack<node_t*> node_stack;
    queue<node_t*> node_q; //too many ueue to type each time
    
    for(size_t i=0;i<unroot.size(); ++i){
        node_stack.push(unroot[i]);
        node_q.push(unroot[i]);
    }
    debug_print("node_stack.size(): %lu", node_stack.size());

    while(!node_stack.empty()){
        node_t* cur = node_stack.top(); node_stack.pop();
        debug_string(cur->to_string().c_str());
        if(cur->_children){
            node_stack.push(cur->_lchild);
            node_stack.push(cur->_rchild);
            node_q.push(cur->_lchild);
            node_q.push(cur->_rchild);
        }
    }

    _size = node_q.size();
    _tree = new node_t[_size];
    size_t cur_index = 0;

    while(!node_q.empty()){
        debug_print("current_index: %lu, _size: %lu", cur_index, _size);
        node_t* cur = node_q.front(); node_q.pop();
        _tree[cur_index] = node_t(*cur);
        debug_string(_tree[cur_index]._label.c_str());
        node_map[cur] = _tree+cur_index;
        cur_index++;
    }

    for(auto &i:unroot){
        debug_string("updating children")
        _unroot.push_back(node_map.at(i));
        _unroot.back()->update_children(node_map);
        debug_string(_unroot.back()->to_string().c_str());
    }
    debug_print("_tree pointer: %p", _tree);
}

tree_t::tree_t(const vector<node_t*>& unroot){
    make_flat_tree(unroot);
}

//takes a tree specified by newick notation, and makes a tree_t
tree_t::tree_t(const string& newick){
    _tree = make_tree_from_newick(newick, _size);
}

tree_t::~tree_t(){
    debug_print("_tree: %p", _tree);
    if(_tree)
        delete[] _tree;
}

tree_t& tree_t::operator=(tree_t t){
    std::swap(_tree, t._tree);
    std::swap(_unroot, t._unroot);
    std::swap(_size, t._size);
    return *this;
}

//we use a label to index map to make the matrix well ordered
//this is so we can do a blind average later on, and not have to worry about 
//the ordering of the array
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
//  debug_print("_tree pointer: %p", _tree);
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

string node_t::to_string(){
    ostringstream ret;
    if(_lchild && _rchild){
        ret<<"("<<_lchild->to_string()
            <<","<<_rchild->to_string()<<")"
            <<_label<<":"<<_weight;
    }
    else{
        ret<<_label<<":"<<_weight;
    }
    return ret.str();
}

string tree_t::to_string() const{
    ostringstream ret;
    debug_print("_tree pointer: %p", _tree);

    if(_unroot.size()>1)
        ret<<"(";

    for(size_t i=0;i<_unroot.size();++i){
        ret<<_unroot[i]->to_string();
        if(i!=_unroot.size()-1)
            ret<<",";
    }
    if(_unroot.size()>1)
        ret<<")";

    return ret.str();
}

string tree_t::print_labels() const{
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

void tree_t::sort(){
    assert_string(_unroot.size() == 3, "the unroot is has a size different than expected");  
    vector<string> label_vector;
    label_vector.reserve(3);
    for(auto && n : _unroot){
        label_vector.push_back(n->sort());
    }

    //make sure the first element is the smallest
    if(label_vector[0] > label_vector[1]){
        std::swap(label_vector[0], label_vector[1]);
        std::swap(_unroot[0], _unroot[1]);
    }
    if(label_vector[0] > label_vector[2]){
        std::swap(label_vector[0], label_vector[2]);
        std::swap(_unroot[0], _unroot[2]);
    }
    //then swap the last two if needed
    if(label_vector[1] > label_vector[2]){
        std::swap(label_vector[1], label_vector[2]);
        std::swap(_unroot[1], _unroot[2]);
    }
}
