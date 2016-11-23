CXX=g++
CFLAGS=-Wall -Wextra -std=c++14
DFLAGS=
IFLAGS=

OBJDIR=obj
SRCDIR=src
TSTDIR=tests

TEST_SOURCES := $(shell find $(TSTDIR) -name '*cpp')
RELEASE_OBJS := $(addprefix $(OBJDIR)/,main.o tree.o newick.o star.o nj.o)
TEST_OBJS := $(addprefix $(OBJDIR)/, $(TEST_SOURCES:$(TSTDIR)/%.cpp=%.o))


all: debug

debug: CFLAGS+= -DDEBUG -g -O0
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

tests: gstar_tests

run-tests: tests
	./gstar_tests -d yes

clean:
	rm -rf obj gstar gstar_tests *.log
