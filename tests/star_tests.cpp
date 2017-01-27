#include "catch.hpp"
#include "../src/star.cpp"
#include <fstream>

TEST_CASE("star, one tree","[star]"){
    std::string t1="(a:1.0,b:1.0);";
    std::vector<std::string> vt;
    vt.push_back(t1);
    star_t s(vt);
    auto star_tree = s.get_tree();
    star_tree.sort();
    REQUIRE(star_tree.to_string()== t1);
}

TEST_CASE("star, two cloned trees", "[star]"){
    std::string t1="(a:1.0,b:1.0);";
    std::vector<std::string> vt;
    vt.push_back(t1);
    vt.push_back(t1);
    star_t s(vt);
    auto star_tree = s.get_tree();
    star_tree.sort();
    REQUIRE(star_tree.to_string()== t1);
}

TEST_CASE("star, two different trees", "[star]"){
    std::string t1 = "(a:1.0,b:1.0);";
    std::string t2 = "(a:2.0,b:0.5);";
    std::vector<std::string> vt;
    vt.push_back(t1);
    vt.push_back(t2);
    star_t s(vt);
    auto star_tree = s.get_tree();
    star_tree.set_weights([](size_t) -> float {return 1.0;});
    star_tree.sort();
    REQUIRE(star_tree.to_string()== t1);
}

TEST_CASE("star, massive trees from ASTRID","[star]"){
    std::string astrid_tree_string = "(((Tree_Shrew,((Rabbit,Pika),(Squirrel,(Guinea_Pig,(Kangaroo_Rat,(Rat,Mouse)))))),((Mouse_Lemur,Galagos),(Tarsier,(Marmoset,(Macaque,(Orangutan,(Gorilla,(Human,Chimpanzee)))))))),((Shrew,Hedgehog),((Megabat,Microbat),((Alpaca,(Pig,(Dolphin,Cow))),(Horse,(Cat,Dog))))),(((Armadillos,Sloth),(Lesser_Hedgehog_Tenrec,(Elephant,Hyrax))),((Wallaby,Opossum),(Platypus,Chicken))));";
    std::ifstream tree_file("tests/song_mammals.424.gene.tre");
    std::string line;
    std::vector<std::string> vt;
    while(std::getline(tree_file, line)){
        vt.push_back(line);
    }
    star_t s(vt);
    auto star_tree = s.get_tree();
    star_tree.set_weights([](size_t) -> float {return 0.0;});
    star_tree.sort();
    REQUIRE(astrid_tree_string == star_tree.to_string());
}
