//newick.cpp
//Ben Bettisworth
//newick notation parser for gstar
//more about newick notation: 
//  https://en.wikipedia.org/wiki/Newick_format
#include "newick.h"
#include "debug.h"

#include <stack>
using std::stack;
#include <cctype>
using std::isdigit;
using std::isalpha;
#include <string>
using std::string;
#include <memory>
using std::shared_ptr;
#include <set>
#include <exception>


std::set<string> label_set;

void skip_whitespace(const string& s, size_t& idx){
    while(isspace(s[idx])!=0 && idx < s.size()) idx++;
}

size_t scan_nodes(const string& s){
    size_t n_nodes=1;
    for(size_t i=0;i<s.size();++i){
        if(s[i]==',') n_nodes++;
        else if(s[i]==')') n_nodes++;
    }
    return n_nodes;
}

string parse_label(const string& newick_string, size_t& idx){
    skip_whitespace(newick_string, idx);
    size_t start, end;
    start = end = idx;
    while(start<newick_string.size()){
        char cur = newick_string[end];
        if(cur<'0' || (cur >= '9'  && cur < 'A') 
                || (cur >'Z' && cur < 'a' && cur!='_') 
                || cur > 'z') break;
        end++;
    }
    idx = end;
    string label = newick_string.substr(start, end-start);
    auto rp = label_set.insert(label);
    if(!rp.second){
        throw std::logic_error("Duplicate taxa label: '" + label 
                + "' at character " + std::to_string(idx));
    }
    return newick_string.substr(start, end-start);
}

double parse_weight(const string& newick_string, size_t& idx){
    size_t start, end;
    start = end = idx;
    while(start<newick_string.size()){
        char cur = newick_string[end];
        if((cur < '0' || cur > '9') && cur!='.' && cur != '-' && cur!='e') break;
        end++;
    }
    idx = end;
    return stod(newick_string.substr(start, end-start));
}

node_t* parse_subtree(const string& newick, size_t& idx, node_t*& next_node){
    stack<node_t*> node_stack;
    node_stack.push(next_node); next_node++;
    node_stack.top()->_weight = 1.0;
    while(idx < newick.size()){
        skip_whitespace(newick, idx);
        char cur = newick[idx];
        debug_print("current character: '%c', current idx: %lu",cur, idx);
        if(cur == '(' || cur == ','){
            debug_string("found new node, pushing onto the stack");
            idx++;
            node_stack.push(next_node);
            next_node++;
            skip_whitespace(newick, idx);
            if(newick[idx] != '('){
                node_stack.top()->_label = parse_label(newick, idx);
            }
            node_stack.top()->_weight = 1.0;
        }
        else if(cur == ':'){
            debug_string("found ':', parsing weight");
            idx++;
            node_stack.top()->_weight = parse_weight(newick, idx);
        }
        else if(cur == ')'){
            idx++;
            debug_string("found ')', popping off the stack");
            node_t* tmp_l = node_stack.top(); node_stack.pop();
            node_t* tmp_r = node_stack.top(); node_stack.pop();
            node_stack.top()->_lchild = tmp_l;
            node_stack.top()->_rchild = tmp_r;
            tmp_l->_parent = node_stack.top();
            tmp_r->_parent = node_stack.top();
            debug_print("done popping off the stack, node_stack.size(): %lu", 
                    node_stack.size());
            if(node_stack.size() == 1){ break; }
        }
        else if(isalpha(cur)){
            debug_string("found a bare alpha, parsing label and breaking");
            node_stack.top()->_label = parse_label(newick, idx);
            node_stack.top()->_weight = 1.0;
            break;
        }
        else{
            assert_string(false, "could not recognize current character");
        }
    }
    debug_string("broke out of while loop, doing cleanup");
    if(newick[idx] == ':'){
        debug_string("found subtree weight, parsing");
        idx++;
        node_stack.top()->_weight = parse_weight(newick, idx);
    }
    debug_string("returning the top of the stack");
    return node_stack.top();
}

//make a tree from a string in newick notation, only a single tree per string.
//specifically, it takes a tree of the form:
//  ( (a:1.0 , b:2.0):2.0 , (c:1.0 , d:1.0):1.0 ):1.0; comments go here
//comments are only supported after the semicolon
std::vector<node_t*> make_tree_from_newick(const string& newick_string, size_t& tree_size){
    debug_string("starting newick parse");
    label_set.clear();

    tree_size = scan_nodes(newick_string);
    debug_print("tree size: %lu", tree_size);
    node_t* tree = new node_t[tree_size];

    size_t idx=0;
    node_t* next_node = tree;
    std::vector<node_t*> unroot;

    debug_string("starting newick parse | while loop");
    while(idx < newick_string.size() && newick_string[idx] != ';'){
        char cur = newick_string[idx];
        debug_print("current char: %c", cur);
        if(cur == '(' || cur == ','){
            idx++;
            unroot.push_back(parse_subtree(newick_string, idx, next_node));
        }
        else if(cur == ')'){
            idx++;
            assert_string(idx < newick_string.size() 
                    && newick_string[idx] == ';', "end of input was not a semicolon");
            break;
        }
        else{
            debug_string("there was an error parsing");
            return {nullptr};
        }
    }
    if(unroot.size() == 1){
        debug_string("found an unroot with size one, contracting the root");
        node_t* t = unroot.front();
        unroot.clear();
        unroot.push_back(t->_lchild);
        unroot.push_back(t->_rchild);
        t->_lchild->_parent = nullptr;
        t->_rchild->_parent = nullptr;
    }
    assert_string(unroot.size() == 2 || unroot.size() == 3, "Unroot is the wrong size");

    for(size_t i=0;i<tree_size;++i){
        tree[i]._children = tree[i]._lchild && tree[i]._rchild;
    }

    return unroot;
}

