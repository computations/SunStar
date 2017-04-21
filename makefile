CXX=clang++
CFLAGS=-Wall -Wextra -std=c++14 
DFLAGS=
IFLAGS=

OBJDIR=obj
SRCDIR=src
TSTDIR=tests
DOCDIR=doc

DFLAGS+= -DGIT_REV=$(shell git describe --tags --always)

TEST_SOURCES := $(shell find $(TSTDIR) -name '*cpp')
RELEASE_OBJS := $(addprefix $(OBJDIR)/,main.o tree.o newick.o star.o nj.o gstar.o)
TEST_OBJS := $(addprefix $(OBJDIR)/, $(TEST_SOURCES:$(TSTDIR)/%.cpp=%.o))

all: release

debug: CFLAGS+= -DDEBUG -DEMIT_DEBUG -g -O0
debug: sunstar

release: CFLAGS+= -DRELEASE -O3 -march=native
release: sunstar

sunstar: $(RELEASE_OBJS)
	$(CXX) $(CFLAGS) -o $@ $^ $(DFLAGS)

sunstar_tests: $(TEST_OBJS)
	$(CXX) $(CFLAGS) -o $@ $^ $(DFLAGS)

%o: CFLAGS+=-c

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CXX) $(CFLAGS) -o $@ $^ $(DFLAGS)

$(OBJDIR)/%.o: $(TSTDIR)/%.cpp | $(OBJDIR)
	$(CXX) $(CFLAGS) -o $@ $^ $(DFLAGS)

$(OBJDIR):
	mkdir $@

run: debug
	./sunstar

tests: CFLAGS+= -DDEBUG -g -O0
tests: sunstar_tests

run-tests: tests
	./sunstar_tests 

verbose-tests: CFLAGS+= -DEMIT_DEBUG
verbose-tests: tests

docs:
	make -C $(DOCDIR)

clean:
	rm -rf obj sunstar sunstar_tests *.log
