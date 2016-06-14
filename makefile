CXX=g++
IFLAGS=
CFLAGS=-Wall -Wextra -O2 -std=c++14

all: gstar

gstar: main.o tree.o newick.o
	$(CXX) $(CFLAGS) -o $@ $^

main.o: main.cpp
	$(CXX) $(CFLAGS) -o $@ $^ -c

tree.o: tree.cpp
	$(CXX) $(CFLAGS) -o $@ $^ -c

newick.o: newick.cpp
	$(CXX) $(CFLAGS) -o $@ $^ -c

clean:
	rm -rf *.o gstar
