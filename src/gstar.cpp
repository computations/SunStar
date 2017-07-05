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
        double ratio = ((double)kv.second)/(trials);
        ret.push_back(std::make_pair(kv.first, ratio));
    }
    return ret;
}

void write_sequence_to_file(const vector<double>& s, string newick_string,
        ofstream& outfile){
    
    outfile<<"{\"tree\":\""<<newick_string<<"\",\"weights\": [";
    for(size_t i = 0; i < s.size(); ++i){
        outfile<<s[i];
        if(i!=s.size()-1){
            outfile<<',';
        }
    }
    outfile<<"]}\n";
}

/*
 * An implementation of the Dirichlet Distribution. Produces a random (math)
 * vector of size len, with the property that the vectors are located on a
 * simplex in dimension len. 
 *
 *  len:    The length of the vector to be returned. Also controls what alpha
 *  should be to be uniform.
 * 
 *  alpha:  Controls how tight the distribution is around the center of the
 *  simplex. If alpha = len, then the distribution is uniform on the simplex.
 *  If alpha is less than len, then the distribution is concentrated near the
 *  edges of the simplex. If alpha is greater than len, the distribution is
 *  concentrated towards the center.
 *
 *  beta:   Don't use this.
 */
vector<double> dirichlet(size_t len, double alpha, double beta=1.0){
#ifndef DEBUG
    std::mt19937 gen((std::random_device())());
    std::gamma_distribution<double> gd(alpha*(beta/(len*beta)), 1.0);
#else
    static std::mt19937 gen((std::random_device())());
    static std::gamma_distribution<double> gd(alpha*(beta/(len*beta)), 1.0);
#endif
    vector<double> ret;
    ret.reserve(len);
    double total = 0.0;
    for(size_t i=0;i<len;++i){
        double tmp = gd(gen);
        total+=tmp;
        ret.push_back(tmp);
    }
    for(auto&& v: ret){
        v /= total;
    }
    return ret;
}

/*
 * Helper function when I just want a uniform vector.
 */
auto dirichlet(size_t len){
    return dirichlet(len, (double)len);
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
vector<std::pair<string, double>> gstar_with_default_schedule(
        star_t& star, const string& logfile, const string& outgroup){

    size_t max_depth = star.get_size();
    vector<double> schedule(max_depth, 0.0);
    unordered_map<string, int> counts;

    ofstream outfile(logfile.c_str());
    outfile<<"using root: '"<<outgroup<<"'"<<std::endl;

    size_t trials = (1<<max_depth) - 1;
    print_progress(0ul, trials);

    for(size_t i=0;i<trials;++i){
        if(i % 100 == 0) {print_progress(i,trials);}
        schedule[0]+=1.0;
        for(size_t k = 0;k<schedule.size()-1;++k){
            if(schedule[k]>1.0){
                schedule[k]=0.0;
                schedule[k+1]+=1.0;
            }
        }
        if(schedule.back()> 1.0){ schedule.back()=0.0;}

        string s = star.get_tree(schedule).set_outgroup(outgroup).
            sort().clear_weights().to_string();
        write_sequence_to_file(schedule, s, outfile);
        counts[s] +=1;
    }

    print_progress(trials, trials);
    finish_progress();
    return make_return_vector(counts, trials);
}

/*
 * A function that does the randomized schedule for GSTAR. In this case, random
 * means that the schedule is pulled from the Dirichlet distribution. This is
 * intended to make a more "uniform" distribution than just a uniform
 * distribution on the reals. 
 *
 * Currently, this is by far the preferred method of using GSTAR. There isn't a
 * lot of reason to use the default schedule
 */
vector<std::pair<string, double>> gstar_with_random_schedule(
        star_t& star, const string& logfile, size_t trials, 
        const string& outgroup){
    size_t max_depth = star.get_size();
    unordered_map<string, int> counts;

    ofstream outfile(logfile.c_str());
    outfile<<"using root: '"<<outgroup<<"'"<<std::endl;

    print_progress(0ul, trials);
    for(size_t i = 0; i < trials; i++){
        if(i % 100 == 0) {print_progress(i,trials);}
        auto schedule = dirichlet(max_depth);
        string s = star.get_tree(schedule).set_outgroup(outgroup).
            sort().clear_weights().to_string();
        write_sequence_to_file(schedule, s, outfile);
        counts[s]+=1;
    }
    print_progress(trials, trials);
    finish_progress();
    return make_return_vector(counts, trials);
}

vector<std::pair<string, double>> gstar(const vector<string>& newick_strings,
        size_t trials, string filename, string outgroup){
    star_t star(newick_strings);
    if(!outgroup.empty()){
        star.set_outgroup(outgroup);
    }
    else{
        outgroup = star.get_first_label();
    }
    if(trials==0){
        return gstar_with_default_schedule(star, filename, outgroup);
    }
    else{
        return gstar_with_random_schedule(star, filename, trials, outgroup);
    }
}
