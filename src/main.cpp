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
    "    -f, --filename      Filename for the set of gene trees in Newick notation\n";
}

int main(int argc, char** argv){
    int c;
    std::string filename;
    while(true){
        static struct option long_options[] = 
        {
            {"filename",    required_argument,  0,   'f'},
            {0,0,0,0}
        };
        int option_index = 0;
        c = getopt_long(argc, argv, "f:", long_options, &option_index);
        if(c==-1){
            break;
        }
        switch(c){
            case 'f':
                filename = string(optarg);
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

    auto trees = gstar(newick_strings);
    for(const auto& kv:trees){
        std::cout<<kv.first<<kv.second<<std::endl;
    }
    
    return 0;
}
