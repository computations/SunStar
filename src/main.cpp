#include "tree.h"
#include "debug.h"
#include <iostream>
using std::cout;
using std::cin;


int main(){
    tree_t t("(a,b);");
    auto lm = t.make_label_map();
    t.calc_distance_matrix(lm);
    return 0;
}
