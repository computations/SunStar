#include "catch.hpp"
#include "../src/gstar.cpp"

#include <cstdio> //for remove

const double epsilon = 1e-9;


TEST_CASE("gstar, functionality test","[gstar]"){
    std::string ns = "(a,b);";
    auto c = gstar({ns});
}

TEST_CASE("gstar, functionality test, larger class", "[gstar]"){
    std::string ns1 = "(a,(b,(c,d)));";
    std::string ns2 = "((a,b),(c,d));";
    auto c = gstar({ns1,ns2});
}

TEST_CASE("gstar, functional test, larger class", "[gstar]"){
    std::string ns = "(((a,b),(c,d)),(((e,f),g),h));";
    auto c = gstar({ns});
}

TEST_CASE("gstar, testing if a file gets made", "[gstar][random]"){
    std::string ns = "(((a,b),(c,d)),(((e,f),g),h));";
    std::string fname = "schedule.log";
    auto c = gstar({ns}, 10, fname);
    std::ifstream ifile(fname.c_str());
    REQUIRE(ifile.good());
    ifile.close();
    std::remove(fname.c_str());
}

TEST_CASE("gstar, functional, with rerooting", "[reroot]"){
    std::string ns = "((a,b),(c,d));";
    std::string fname = "schedule.log";
    auto c = gstar({ns}, 10, fname, "a");
    std::remove(fname.c_str());
}

TEST_CASE("gstar, small regression", "[gstar][regression]"){
    std::string s1 = "((a,((b,c),k)),e);";
    std::string s2 = "((b,((a,c),k)),e);";
    auto trees = gstar({s1,s2});
    std::unordered_map<std::string, double> tree_map;
    for(auto &t : trees){
        tree_map[t.first] = t.second;
    }
    std::vector<std::pair<std::string, double>>ans {
        {"((((a,e),k),c),b);", 0.3 + 1.0/6.0},
        {"(((a,(e,k)),c),b);", 0.1 + 1.0/6.0},
        {"(((a,e),(c,k)),b);", 0.1 + 1.0/6.0}};
    for(auto& a : ans){
        REQUIRE(tree_map.count(a.first));
        REQUIRE(tree_map[a.first] - a.second < epsilon);
        REQUIRE(-(tree_map[a.first] - a.second) < epsilon);
    }
}

TEST_CASE("gstar, checking random_schedule", "[gstar][random]"){
    std::string s1 = "((a,((b,c),k)),e);";
    std::string s2 = "((b,((a,c),k)),e);";
    auto trees = gstar({s1,s2}, 10);
}

TEST_CASE("dirichlet random numbers 1, dim 3", "[gstar][dirichlet]"){
    auto dv = dirichlet(3);
    double total = 0.0;
    for(auto&& v : dv){
        total+=v;
    }
    REQUIRE( (total-1.0) < epsilon);
    REQUIRE(-(total-1.0) < epsilon);
}

TEST_CASE("dirichlet random numbers 1, dim 1000", "[gstar][dirichlet]"){
    auto dv = dirichlet(1000);
    double total = 0.0;
    for(auto&& v : dv){
        total+=v;
    }
    REQUIRE( (total-1.0) < epsilon);
    REQUIRE(-(total-1.0) < epsilon);
}

TEST_CASE("dirichlet random numbers 2, average test, dim 3", "[gstar][dirichlet]"){
    const size_t TRIALS = 1e5;
    const int DIM = 3;
    const double EXPECTED = 1.0/DIM;
    vector<double> avg_vec(DIM,0.0);
    for(size_t i =  0; i <TRIALS; ++i){
        auto dv = dirichlet(DIM);
        for(size_t j = 0; j< dv.size(); ++j){
            avg_vec[j]+=dv[j];
        }
    }
    for(auto&& i:avg_vec){
        i/=TRIALS;
    }
    //We could run the number of trials required to get to the "real" epsilon
    //here. But random trials has a sublinear growth in accuracy, so that takes
    //forever. We are talking 1e12 or more. Therefore, we relax the constraints
    //for this test.
    for(auto&& i:avg_vec){
        REQUIRE( (i - EXPECTED) < epsilon*1e6);
        REQUIRE(-(i - EXPECTED) < epsilon*1e6);
    }
}

TEST_CASE("dirichlet random numbers 2, average test, dim 100", "[gstar][dirichlet]"){
    const int TRIALS = 1e5;
    const int DIM = 100;
    const double EXPECTED = 1.0/DIM;
    vector<double> avg_vec(DIM,0.0);
    for(size_t i =  0; i <TRIALS; ++i){
        auto dv = dirichlet(DIM);
        for(size_t j = 0; j< dv.size(); ++j){
            avg_vec[j]+=dv[j];
        }
    }
    for(auto&& i:avg_vec){
        i/=TRIALS;
    }
    for(auto&& i:avg_vec){
        REQUIRE( (i - EXPECTED) < epsilon*1e6);
        REQUIRE(-(i - EXPECTED) < epsilon*1e6);
    }
}
