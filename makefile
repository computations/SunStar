CXX=g++
CFLAGS=-Wall -Wextra -O2 -std=c++14
IFLAGS=

OBJDIR=obj
SRCDIR=src

OBJS := $(addprefix $(OBJDIR)/,main.o tree.o newick.o)

all: gstar

gstar: $(OBJS)
	$(CXX) $(CFLAGS) -o $@ $^

%o: CFLAGS+=-c

$(OBJDIR)/main.o: $(SRCDIR)/main.cpp | $(OBJDIR)
	$(CXX) $(CFLAGS) -o $@ $^

$(OBJDIR)/tree.o: $(SRCDIR)/tree.cpp | $(OBJDIR)
	$(CXX) $(CFLAGS) -o $@ $^

$(OBJDIR)/newick.o: $(SRCDIR)/newick.cpp | $(OBJDIR)
	$(CXX) $(CFLAGS) -o $@ $^

$(OBJDIR):
	mkdir $@

run: gstar
	./gstar

clean:
	rm -rf obj gstar
