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
    size_t n_nodes=0;
    for(size_t i=0;i<s.size();++i){
        if(s[i]==',') n_nodes++;
        else if(s[i]==')') n_nodes+=2;
    }
    return n_nodes;
}

//Currently, there needs to be a tree on each string line for this parser to work
//If you put multiple trees in, then it has undefined behavior
node_t* make_tree_from_newick(const string& newick_string, size_t& tree_size){
    tree_size = scan_nodes(newick_string);
    node_t* tree = new node_t[tree_size];
    stack<node_t*> node_stack;

    size_t idx=0;
    size_t current_node = 0;
    while(idx < newick_string.size()){
        if(newick_string[idx] == '(' || newick_string[idx] == ','){
            node_stack.push(tree+current_node);
            current_node++; idx++;
        }
        else if(newick_string[idx] == ')'){
            node_t* tmp1 = node_stack.top(); node_stack.pop();
            node_t* tmp2 = node_stack.top(); node_stack.pop();
            node_stack.top()->_lchild = tmp1 - tree;
            node_stack.top()->_rchild = tmp2 - tree;
            node_stack.top()->_children = true;
            idx++;
        }
        else if(newick_string[idx] == ';'){
            break;
        }
        else{
            idx = skip_whitespace(newick_string, idx);           
            size_t j = idx;
            while(j<newick_string.size() 
                    && newick_string[j] != ','
                    && newick_string[j] != ':'){ 
                j++;
            }
            node_stack.top()->_label = newick_string.substr(idx, j-idx);
            idx = j = j+1;
            while(j < newick_string.size() 
                    && (isdigit(newick_string[j] != 0 || newick_string[j]=='.'))){
                ++j;
            }
            if(idx != j)
                node_stack.top()->_weight = stof(newick_string.substr(idx, j-idx));
            idx = j;
        }
    }

    assert_string(node_stack.size()==1, "There was a parse error");
    assert_string(node_stack.top() == tree, "the root wasn't the top of the stack");

    //if the root isn't at the front, we need to make it that way
    return tree;
}
/*
node_t* make_tree_from_newick(const string& s, size_t& tree_size){
    stack<ptr_node_t*> node_stack;
    node_stack.push(new ptr_node_t);
    size_t index=0;

    while(index < s.size()){
        index = skip_whitespace(s, index);
        if(s[index] == '(' || s[index]==','){
            ptr_node_t* tmp = new ptr_node_t;
            index++;
            node_stack.push(tmp);
        }
        else if(s[index]==')'){
            ptr_node_t* tmp1 = node_stack.top(); node_stack.pop();
            ptr_node_t* tmp2 = node_stack.top(); node_stack.pop();
            node_stack.top()->_lchild=tmp1;
            node_stack.top()->_rchild=tmp2;
            index++;
        }
        else if(s[index]==';') break;
        else{
            index = skip_whitespace(s, index);
            size_t j=index;
            while(j<s.size() && s[j]!=',' && s[j]!=':') ++j;
            node_stack.top()->_label = s.substr(index, j-index);
            index = j = j+1;
            while(j < s.size() && (isdigit(s[j])!=0 || s[j]=='.')) ++j;
            if(index != j)
                node_stack.top()->_weight = stof( s.substr(index, j-index));
            index = j;
        }
    }
    
    //one final bit of setup to make the whole tree trivalant.
    //We make a special root taxon that gets joined to the top of the stack
    //Its "parent" is the one child it has
    debug_print("stack size at end: %lu", node_stack.size());

    node_stack.top()->_weight = 0;

    return convert_to_packed_tree(node_stack.top(), tree_size);
}*/
