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
#include <algorithm>
//for std::sort
#include <getopt.h>

//dumb hack to get the progress bar crap to work
bool __PROGRESS_BAR_FLAG__=true;

void print_usage(){
    std::cout<<
"Usage: gstar [options]\n"<<
"Application Options:\n"<<
"    -f, --filename [FILE]\n"<<
"           Filename of a file that contains a list of newick strings,\n"<<
"           representing gene trees to summarize\n"<<
"    -t, --trials [NUMBER]\n"<<
"           Number of trials. This will induce the random schedule\n"<<
"    -r, --required-ratio [RATIO]\n"<<
"           The minimum ratio required for a tree to be output at the end.\n"<<
"           Trees below this threshold will be silenced. Intended to be bewteen\n"<<
"           0 and 1\n"<<
"    -l, --logfile [FILE]\n"<<
"           File to log the sequences to (defaults to schedule.log)\n"<<
"    -o, --outgroup [STRING]\n"<<
"           Taxa label of the outgroup of the gene trees\n"<<
"    -s, --silent\n"<<
"           Silence the progress bar, only output results\n";
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
    double threshold=0;

    while(true){
        static struct option long_options[] = 
        {
            {"silent",    no_argument, 0, 's'},
            {"require-ratio", required_argument, 0, 'r'},
            {"filename",    required_argument,  0,   'f'},
            {"outgroup",    required_argument,  0,   'o'},
            {"logfile",    required_argument,  0,   'l'},
            {"trials",    required_argument,  0,   't'},
            {0,0,0,0}
        };
        int option_index = 0;
        c = getopt_long(argc, argv, "sf:o:t:l:p:", long_options, &option_index);
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
            case 'p':
                threshold = std::stod(optarg);
                break;
            case 's':
                turn_off_progress();
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

    if(!newick_string_file){
        std::cerr<<"Failed to open newick string file"<<std::endl;
        return 1;
    }

    while(std::getline(newick_string_file, line)){
        newick_strings.emplace_back(line);
    }

    if(outgroup.empty() && !check_rooted(newick_strings)){
        std::cout<<"Error, trees are not rooted, or an outgroup has not been specified"<<std::endl;
        return 1;
    }

    auto trees = gstar(newick_strings, logfile, trials, outgroup);
    //sort the trees
    
    auto pc_lambda = [](auto lhs, auto rhs){
        return lhs.second>rhs.second;
    };

    std::sort(trees.begin(), trees.end(), pc_lambda);

    for(const auto& kv:trees){
        if(kv.second<threshold) break;
        std::cout<<kv.first<<kv.second<<std::endl;
    }
    
    return 0;
}
