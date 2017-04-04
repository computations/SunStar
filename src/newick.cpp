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

/*
 * Scan the string for number of nodes in the tree. Used to allocate memory at
 * the start of make_tree_from_newick. It does this by counting commas and
 * closing parens
 */
size_t scan_nodes(const string& s){
    size_t n_nodes=1;
    for(size_t i=0;i<s.size();++i){
        if(s[i]==',') n_nodes++;
        else if(s[i]==')') n_nodes++;
    }
    return n_nodes;
}

/* 
 * Pulls a valid label out of a string, starting at some index. Skips
 * whitespace. Valid identfies are any combination of digits, underscores and
 * alphas.
 *
 * Also checks if the labels on a tree are unique. If they are not, throws an
 * exception.
 *
 * Returns a string contianing the label
 */
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

/*
 * Parses the weight for a node. Weights from parsed trees are often discarded,
 * but we parse them anyways. Parsing is done via std::stod()
 */
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

/*
 * Parses the subtree of a newick tree. This abstraction is nessicary because
 * the newick tree could be rooted or unrooted, which will cause the top level
 * to have 2 or 3 taxa present. This means that we represent the tree root as a
 * vector of nodes. Therefore, its easier to parse an actually rooted tree with
 * one root as a subtree, then it is to deal with parsing the whole tree as an
 * unrooted tree.
 *
 * Returns a pointer to the node at the top of a subtree.
 */
node_t* parse_subtree(const string& newick, size_t& idx, node_t*& next_node){
    stack<node_t*> node_stack;
    node_stack.push(next_node); next_node++;
    node_stack.top()->_weight = 1.0;
    while(idx < newick.size()){
        skip_whitespace(newick, idx);
        char cur = newick[idx];
        debug_print("current character: '%c', current idx: %lu",cur, idx);
        /*
         * Found a new taxa, so we need to push a new node onto the stack,
         * ready to be processed. Then, we need to find the label, if it
         * exists. Set the wieght to 1 by default.
         */
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
        /*
         * Found a new weight, parse it and slap it on the node on the top of
         * the stack.
         */
        else if(cur == ':'){
            debug_string("found ':', parsing weight");
            idx++;
            node_stack.top()->_weight = parse_weight(newick, idx);
        }
        /*
         * Closing paren means that we have a finished subtree. Pop the top two
         * items off the stack, and make them the children of the new top. If
         * the stack only has one item left on it, then we are done, and we
         * break out of the loop.
         */
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
        /*
         * Handle the corner case where the subtree is just an label. This
         * occurs when the tree is something like (a,b); In this case, we can
         * just make a node with that label and push it on the stack. We don't
         * need to parse the weight because there is clean up to handle that
         * outside the while loop
         */
        else if(isalpha(cur)){
            debug_string("found a bare alpha, parsing label and breaking");
            node_stack.top()->_label = parse_label(newick, idx);
            node_stack.top()->_weight = 1.0;
            break;
        }
        /*
         * The character is something unexpected, time to give up!
         */
        else{
            throw std::runtime_error("Error in parsing, did not expect current character");
        }
    }
    /*
     * We broke out when we found either a singleton subtree, or a finished
     * subtree. There might still be a weight after the subtree ends, so parse
     * that if its there
     */
    debug_string("broke out of while loop, doing cleanup");
    if(newick[idx] == ':'){
        debug_string("found subtree weight, parsing");
        idx++;
        node_stack.top()->_weight = parse_weight(newick, idx);
    }
    debug_string("returning the top of the stack");
    /*
     * Now that is done, we return the top of the stack
     */
    return node_stack.top();
}

/* 
 * Make a tree from a string in newick notation, only a single tree per string.
 * specifically, it takes a tree of the form:
 *  ( (a:1.0 , b:2.0):2.0 , (c:1.0 , d:1.0):1.0 ):1.0; comments go here
 * comments are only supported after the semicolon
 */
std::vector<node_t*> make_tree_from_newick(const string& newick_string, size_t& tree_size){
    debug_string("starting newick parse");
    label_set.clear();

    tree_size = scan_nodes(newick_string);
    debug_print("tree size: %lu", tree_size);
    node_t* tree = new node_t[tree_size];

    size_t idx=0;
    node_t* next_node = tree;
    std::vector<node_t*> unroot;

    /*
     * Here we need to be able to parse at least 3 subtrees. This while loop
     * will parse as many as we can find. Every iteration produces a new
     * subtree, via parse_subtree(). Each subtree just gets pushed onto
     * the unroot, and the unroot is returned.
     */
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
            throw std::runtime_error("Could not recognize the next character");
        }
    }
    /*
     * There was an early version of this code that would produce "strictly"
     * rooted trees, by which I mean that there was only one node in the
     * unroot. Later, lots of code changed, and it became easier to make the
     * assumption that I was always dealing with trees that had an unroot with
     * more than one node in it. This chekc is here to ensure that is true.
     */
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

