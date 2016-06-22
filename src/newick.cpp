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

inline size_t skip_whitespace(const string& s, size_t index){
    while(isspace(s[index])!=0 && index < s.size()) index++;
    return index;
}

size_t scan_nodes(const string& s){
    size_t n_nodes=1;
    for(size_t i=0;i<s.size();++i){
        if(s[i]==',') n_nodes++;
        else if(s[i]==')') n_nodes++;
    }
    return n_nodes;
}

//Currently, there needs to be a tree on each string line for this parser to work
//If you put multiple trees in, then it has undefined behavior
node_t* make_tree_from_newick(const string& newick_string, size_t& tree_size){
    debug_string("starting newick parse");
    tree_size = scan_nodes(newick_string);
    debug_print("tree size: %lu", tree_size);
    node_t* tree = new node_t[tree_size];
    stack<node_t*> node_stack;
    node_stack.push(tree);
    tree[0]._parent = tree;

    debug_print("tree size is : %lu", tree_size);

    size_t idx=0;
    size_t current_node = 1;
    debug_string("starting newick parse| while loop");
    while(idx < newick_string.size()){
        debug_print("current character: %c", newick_string[idx]);

        if(newick_string[idx] == '(' || newick_string[idx] == ','){
            debug_string("found new taxa, pushing new node");
            node_stack.push(tree+current_node);
            current_node++; idx++;
        }
        else if(newick_string[idx] == ')'){
            debug_string("closing a subtree");
            node_t* tmp1 = node_stack.top(); node_stack.pop();
            node_t* tmp2 = node_stack.top(); node_stack.pop();
            node_stack.top()->_lchild = tmp1;
            node_stack.top()->_rchild = tmp2;
            tmp1->_parent = node_stack.top();
            tmp2->_parent = node_stack.top();
            idx++;
        }
        else if(newick_string[idx] == ';'){
            debug_string("semicolon encountered, its over");
            break;
        }
        else{
            idx = skip_whitespace(newick_string, idx);
            size_t j = idx;
            while(j<newick_string.size() 
                    && newick_string[j] != ','
                    && newick_string[j] != ':'
                    && newick_string[j] != ';'
                    && newick_string[j] != ')'){ 
                j++;
            }
            node_stack.top()->_label = newick_string.substr(idx, j-idx);
            idx = j;
            if(newick_string[idx] == ':'){
                j = ++idx;
                debug_print("parsing number: current character: %c", newick_string[idx]);
                while(j < newick_string.size() 
                        && (isdigit(newick_string[j]) != 0 || newick_string[j]=='.')){
                    ++j;
                }
                debug_print("parsed number: %s" ,newick_string.substr(idx,j-idx).c_str());
                node_stack.top()->_weight = stof(newick_string.substr(idx, j-idx));
                idx = j;
            }
            else{
                node_stack.top()->_weight = 0.0;
            }
        }
    }

    for(size_t i=0;i<tree_size;++i){
        tree[i]._children = tree[i]._lchild || tree[i]._rchild;
    }

    debug_string("done parsing, going into asserts");

    assert_string(node_stack.size() == 1, "There was a parse error");
    assert_string(node_stack.top() == tree, "the root wasn't the top of the stack");

    //if the root isn't at the front, we need to make it that way
    return tree;
}
