#include "tree.h"
#include "star.h"
#include "debug.h"
#include <iostream>
using std::cout;
using std::endl;
#include <vector>
using std::vector;


int main(){
    tree_t t1("((a:1.0,b:1.0),(c:1.0,d:1.0));");
    tree_t t2("(((a:1.0,d:1.0),c:1.0),d:1.0);");

    cout<<t1.to_string()<<endl;
    cout<<t1.print_labels()<<endl;
    cout<<t2.to_string()<<endl;
    cout<<t2.print_labels()<<endl;

    vector<tree_t> tree_vector;
    tree_vector.push_back(t1);
    tree_vector.push_back(t2);

    auto dists = calc_average_distances(tree_vector);

    for(size_t i = 0; i < 4; ++i){
        for(size_t j = 0; j < 4; ++j){
            cout<<dists[i*4+j]<<" ";
        }
        cout<<endl;
    }

    t1.set_weights([](size_t d){ return d*2.56;});
    cout<<t1.to_string()<<endl;
    t1.set_weights({1,2,3,4,5,6,7,8});
    cout<<t1.to_string()<<endl;

    return 0;
}
