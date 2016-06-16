#include "tree.h"
#include "debug.h"
#include <iostream>
using std::cout;
using std::endl;


int main(){
    tree_t t("(a,b)k;");
    auto lm = t.make_label_map();
    cout<<t.to_string()<<endl;

    return 0;
}
