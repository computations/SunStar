#include "catch.hpp"
#include "../src/gstar.cpp"

#include <cstdio> //for remove


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
