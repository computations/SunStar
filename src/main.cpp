#include "tree.h"
#include "star.h"
#include "debug.h"
#include "nj.h"
#include <iostream>
using std::cout;
using std::endl;
#include <vector>
using std::vector;
#include <string>
using std::string;


int main(){

    vector<string> tree_strings = 
        {"((a:1.0,b:1.0):1.0,(c:1.0,d:1.0):1.0);",
         "(((a:1.0,b:1.0):1.0,c:1.0):1.0,d:1.0);"};

    //tree_t tests
    tree_t t1(tree_strings[0]);
    tree_t t2(tree_strings[1]);

    cout<<t1.to_string()<<endl;
    cout<<t1.print_labels()<<endl;
    cout<<t2.to_string()<<endl;
    cout<<t2.print_labels()<<endl;

    t1.set_weights([](size_t d){ return d*2.56;});
    cout<<t1.to_string()<<endl;
    t1.set_weights({1,2,3,4,5,6,7,8});
    cout<<t1.to_string()<<endl;

    //star tests
    star_t star(tree_strings);
    auto t = star.get_tree();

    cout<<t.print_labels()<<endl;
    cout<<t.to_string()<<endl;
    t.sort();
    cout<<t.to_string()<<endl;

    return 0;
}
