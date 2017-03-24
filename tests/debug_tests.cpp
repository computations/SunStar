#include "catch.hpp"
#include "../src/debug.h"
bool __PROGRESS_BAR_FLAG__=false;

TEST_CASE("debug, print clock", "[debug]"){
    print_clock;
    fprintf(stderr, "\n");
}

TEST_CASE("debug, debug_print", "[debug]"){
    debug_print("here is a 2.0: %f", 2.0);
}

TEST_CASE("debug, debug_string", "[debug]"){
    debug_string("here is a test string");
}

TEST_CASE("debug, print_progress", "[.][debug]"){
    print_progress(0ul,100ul);
    fprintf(stdout,"\n");
    print_progress(50ul,100ul);
    fprintf(stdout,"\n");
    print_progress(100ul,100ul);
    fprintf(stdout,"\n");
}
