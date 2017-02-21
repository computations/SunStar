//gstar.cpp
//Ben Bettisworth
//2017-02-20
#include "gstar.h"
#include <unordered_map>
using std::unordered_map;
#include <string>
using std::string;
#include <vector>
using std::vector;
#include <utility>

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
vector<std::pair<string, double>> gstar(const vector<string>& newick_strings){
    star_t star(newick_strings);
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

    vector<std::pair<string, double>> ret;
    for(auto&& kv : counts){
        double ratio = ((double)kv.second)/(max_depth+1);
        ret.push_back(std::make_pair(kv.first, ratio));
    }
    return ret;
}
