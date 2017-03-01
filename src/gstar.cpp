//gstar.cpp
//Ben Bettisworth
//2017-02-20
#include "gstar.h"
#include "debug.h"
#include <unordered_map>
using std::unordered_map;
#include <string>
using std::string;
#include <vector>
using std::vector;
#include <utility>
#include <random>
#include <fstream>
using std::ofstream;

vector<std::pair<string, double>> make_return_vector(
        const unordered_map<string, int>& counts, size_t trials){
    vector<std::pair<string, double>> ret;
    for(auto&& kv : counts){
        double ratio = ((double)kv.second)/(trials+1);
        ret.push_back(std::make_pair(kv.first, ratio));
    }
    return ret;
}

void write_sequence_to_file(const vector<double>& s, string newick_string,
        ofstream& outfile){
    
    outfile<<"This newick string: "<<newick_string<<" using the sequence:\n\t";
    for(size_t i = 0; i < s.size(); ++i){
        outfile<<s[i];
        if(i!=s.size()-1){
            outfile<<',';
        }
    }
    outfile<<"\n";
}

/*
 * Since we can set the weights on the tree to various weights, as long as we
 * follow a schedule, we can try and find an some error if we play with the
 * schedule a bit. Specifically, we want to set some some weights on a set of
 * gene trees to zero, and make a species tree out of that set. If we find that
 * changing the schedule changes the resulting tree, we can view that as having
 * something like less confidence in the result.
 *
 * This function returns a vector of pairs. The first of the pair is the newick
 * string associated with the tree. The second in the pair is the "support" for
 * the tree as a ratio. The ratio is the number of times the tree was produce
 * over the total trials.
 */
vector<std::pair<string, double>> gstar(const vector<string>& newick_strings, 
        string outgroup){
    star_t star(newick_strings);
    if(!outgroup.empty()){
        star.set_outgroup(outgroup);
    }
    size_t max_depth = star.get_size();
    vector<double> schedule(max_depth, 1.0);
    unordered_map<string, int> counts;

    {
        string s = star.get_tree(schedule).sort().clear_weights().to_string();
        counts[s] +=1;
    }

    for(size_t i=0;i<max_depth;++i){
        schedule[i] = 0.0;
        string s = star.get_tree(schedule).sort().clear_weights().to_string();
        counts[s] +=1;
    }

    return make_return_vector(counts, max_depth);
}

vector<std::pair<string, double>> gstar(const vector<string>& newick_strings,
        const string& filename, size_t trials, string outgroup){
    if(trials==0){
        return gstar(newick_strings, outgroup);
    }
    star_t star(newick_strings);
    if(!outgroup.empty()){
        star.set_outgroup(outgroup);
    }
    size_t max_depth = star.get_size();
    vector<double> schedule (max_depth, 0.0);
    unordered_map<string, int> counts;

    ofstream outfile(filename.c_str());

    std::mt19937 gen((std::random_device())());
    std::normal_distribution<> d(1.0,1.0);

    for(size_t i = 0; i < trials; i++){
        double tmp = 0.0;
        //Need to randomize the schedule
        while((tmp = d(gen)) < 0){}
        debug_print("setting schedule[0]: %f", tmp);
        schedule[0] = tmp;
        for(size_t k = 1; k < schedule.size(); ++k){
            while((tmp = d(gen)) < 0){ }
            debug_print("setting schedule[%lu]: %f, tmp: %f", k, schedule[k-1] + tmp, tmp);
            schedule[k] = schedule[k-1]+tmp;
        }
        string s = star.get_tree(schedule).sort().clear_weights().to_string();
        write_sequence_to_file(schedule, s, outfile);
        counts[s]+=1;
    }
    return make_return_vector(counts, trials);
}
