#include "tree.h"
#include "debug.h"
#include <iostream>
using std::cout;
using std::endl;


int main(){
    tree_t t("(a:1.0,b:1.0);");
    auto lm = t.make_label_map();
    cout<<t.to_string()<<endl;
    cout<<t.print_labels()<<endl;
    size_t row_size = lm.size();
    size_t matrix_size  = row_size*row_size;
    float* dists = new float[matrix_size];
    for(size_t i=0; i<matrix_size; ++i){
        dists[i] = 0.0;
    }
    t.calc_distance_matrix(lm, dists);
    for(size_t i=0;i<row_size;++i){
        for(size_t j=0;j<row_size;++j){
            cout<<dists[i*row_size+j]<<" ";
        }
        cout<<'\n';
    }
    return 0;
}
