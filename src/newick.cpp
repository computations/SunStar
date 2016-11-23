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
#include <string>
using std::string;
#include <memory>
using std::shared_ptr;

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
    size_t start, end;
    start = end = idx;
    while(start<newick_string.size()){
        char cur = newick_string[end];
        if(cur<'0' || (cur >= '9'  && cur < 'A') || (cur >'Z' && cur < 'a') 
                || cur > 'z') break;
        end++;
    }
    idx = end;
    return newick_string.substr(start, end-start);
}

float parse_weight(const string& newick_string, size_t& idx){
    size_t start, end;
    start = end = idx;
    while(start<newick_string.size()){
        char cur = newick_string[end];
        if((cur < '0' || cur > '9') && cur!='.' && cur != '-') break;
        end++;
    }
    idx = end;
    return stof(newick_string.substr(start, end-start));
}

//make a tree from a string in newick notation, only a single tree per string.
//specifically, it takes a tree of the form:
//  ( (a:1.0 , b:2.0):2.0 , (c:1.0 , d:1.0):1.0 ):1.0; comments go here
//comments are only supported after the semicolon
node_t* make_tree_from_newick(const string& newick_string, size_t& tree_size){
    debug_string("starting newick parse");

    tree_size = scan_nodes(newick_string);
    debug_print("tree size: %lu", tree_size);
    node_t* tree = new node_t[tree_size];

    stack<node_t*> node_stack;
    node_stack.push(tree);
    tree[0]._parent = tree;

    size_t idx=0;
    size_t current_node = 1;
    debug_string("starting newick parse| while loop");
    while(idx < newick_string.size()){
        skip_whitespace(newick_string, idx);
        char cur = newick_string[idx];
        debug_print("current character: %c, current idx: %lu",cur, idx);

        if(cur == '(' || cur == ','){
            debug_string("found new taxa, pushing new node");
            node_stack.push(tree+current_node);
            current_node++; idx++;
            node_stack.top()->_label = parse_label(newick_string, idx);
            debug_print("new taxa label: %s", node_stack.top()->_label.c_str());
            debug_string("setting weight on new taxa to 1.0");
            node_stack.top()->_weight = 1.0;
        }
        else if(cur == ':'){
            debug_string("found colon, parsing float");
            idx++;
            node_stack.top()->_weight = parse_weight(newick_string, idx);
            debug_print("new taxa weight: %f", node_stack.top()->_weight);
            
        }
        else if(cur == ')'){
            debug_string("found ')', finishing subtree");
            node_t* tmp_l = node_stack.top(); node_stack.pop();
            node_t* tmp_r = node_stack.top(); node_stack.pop();
            node_stack.top()->_lchild = tmp_l;
            node_stack.top()->_rchild = tmp_r;
            tmp_l->_parent = node_stack.top();
            tmp_r->_parent = node_stack.top();
            idx++;
        }
        else if(cur == ';'){
            break;
        }
        else{
            assert_string(false, "Didn't recongize the current character in while parsing newick notation");
        }
    }

    for(size_t i=0;i<tree_size;++i){
        tree[i]._children = tree[i]._lchild || tree[i]._rchild;
    }

    debug_string("done parsing, going into asserts");

    assert_string(node_stack.size() == 1, "There was a parse error");
    assert_string(node_stack.top() == tree, "the root wasn't the top of the stack");

    return tree;
}

