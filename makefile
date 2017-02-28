CXX=clang++
CFLAGS=-Wall -Wextra -std=c++14 
DFLAGS=
IFLAGS=

OBJDIR=obj
SRCDIR=src
TSTDIR=tests

TEST_SOURCES := $(shell find $(TSTDIR) -name '*cpp')
RELEASE_OBJS := $(addprefix $(OBJDIR)/,main.o tree.o newick.o star.o nj.o gstar.o)
TEST_OBJS := $(addprefix $(OBJDIR)/, $(TEST_SOURCES:$(TSTDIR)/%.cpp=%.o))

all: debug

debug: CFLAGS+= -DDEBUG -DEMIT_DEBUG -g -O0
debug: gstar

release: CFLAGS+= -DRELEASE -Ofast -march=native
release: gstar

gstar: $(RELEASE_OBJS)
	$(CXX) $(CFLAGS) -o $@ $^ $(DFLAGS)

gstar_tests: $(TEST_OBJS)
	$(CXX) $(CFLAGS) -o $@ $^ $(DFLAGS)

%o: CFLAGS+=-c

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CXX) $(CFLAGS) -o $@ $^ $(DFLAGS)

$(OBJDIR)/%.o: $(TSTDIR)/%.cpp | $(OBJDIR)
	$(CXX) $(CFLAGS) -o $@ $^ $(DFLAGS)

$(OBJDIR):
	mkdir $@

run: debug
	./gstar

tests: CFLAGS+= -DDEBUG -g -O0
tests: gstar_tests

run-tests: tests
	./gstar_tests 

verbose-tests: CFLAGS+= -DEMIT_DEBUG
verbose-tests: tests

clean:
	rm -rf obj gstar gstar_tests *.log
