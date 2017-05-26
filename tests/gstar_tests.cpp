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

TEST_CASE("gstar, testing if a file gtes made", "[gstar][random]"){
    std::string ns = "(((a,b),(c,d)),(((e,f),g),h));";
    std::string fname = "schedule.log";
    auto c = gstar({ns}, fname, 10);
    std::ifstream ifile(fname.c_str());
    REQUIRE(ifile.good());
    ifile.close();
    std::remove(fname.c_str());
}

TEST_CASE("gstar, functional, with rerooting", "[reroot]"){
    std::string ns = "((a,b),(c,d));";
    auto c = gstar({ns}, "a");
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
