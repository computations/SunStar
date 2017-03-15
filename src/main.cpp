#include "tree.h"
#include "star.h"
#include "debug.h"
#include "nj.h"
#include "gstar.h"
#include <iostream>
using std::cout;
using std::endl;
#include <vector>
using std::vector;
#include <string>
using std::string;
#include <fstream>
using std::ifstream;
#include <getopt.h>

void print_usage(){
    std::cout<<
    "Usage: gstar [options]\n"<<
    "Application Options:\n"<<
    "    -f, --filename      Filename for the set of gene trees in Newick notation\n"<<
    "    -t, --trials        Number of trials\n"<<
    "    -l, --logfile       File to log the sequences to (default schedule.log)\n"<<
    "    -o, --outgroup      Taxa label of the outgroup of the gene trees\n";
}

bool check_rooted(const vector<string>& nstrings){
    for(auto& t : nstrings){
        tree_t tmp(t);
        if(!tmp.is_rooted()) return false;
    }
    return true;
}

int main(int argc, char** argv){
    int c;
    //argument buffers
    std::string filename;
    std::string outgroup;
    std::string logfile="schedule.log";
    size_t trials=0;

    while(true){
        static struct option long_options[] = 
        {
            {"filename",    required_argument,  0,   'f'},
            {"outgroup",    required_argument,  0,   'o'},
            {"logfile",    required_argument,  0,   'l'},
            {"trials",    required_argument,  0,   't'},
            {0,0,0,0}
        };
        int option_index = 0;
        c = getopt_long(argc, argv, "f:o:t:l:", long_options, &option_index);
        if(c==-1){
            break;
        }
        switch(c){
            case 'f':
                filename = string(optarg);
                break;
            case 'o':
                outgroup = string(optarg);
                break;
            case 'l':
                logfile = string(optarg);
                break;
            case 't':
                trials = std::stoi(optarg);
                break;
            default:
                break;
        }
    }
    if(filename.empty()){
        print_usage();
        return 1;
    }
    ifstream newick_string_file(filename.c_str());
    std::string line;
    std::vector<std::string> newick_strings;

    while(std::getline(newick_string_file, line)){
        newick_strings.emplace_back(line);
    }

    if(outgroup.empty() && !check_rooted(newick_strings)){
        std::cout<<"Error, trees are not rooted, or an outgroup has not been specified"<<std::endl;
        return 1;
    }

    auto trees = gstar(newick_strings, logfile, trials, outgroup);
    for(const auto& kv:trees){
        std::cout<<kv.first<<kv.second<<std::endl;
    }
    
    return 0;
}
