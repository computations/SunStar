CXX=g++
CFLAGS=-Wall -Wextra -std=c++14
DFLAGS=
IFLAGS=

OBJDIR=obj
SRCDIR=src

OBJS := $(addprefix $(OBJDIR)/,main.o tree.o newick.o)

all: debug

debug: CFLAGS+= -DDEBUG -g
debug: gstar

release: CFLAGS+= -DRELEASE -O2
release: gstar

gstar: $(OBJS)
	$(CXX) $(CFLAGS) -o $@ $^ $(DFLAGS)

%o: CFLAGS+=-c

$(OBJDIR)/main.o: $(SRCDIR)/main.cpp | $(OBJDIR)
	$(CXX) $(CFLAGS) -o $@ $^ $(DFLAGS)

$(OBJDIR)/tree.o: $(SRCDIR)/tree.cpp | $(OBJDIR)
	$(CXX) $(CFLAGS) -o $@ $^ $(DFLAGS)

$(OBJDIR)/newick.o: $(SRCDIR)/newick.cpp | $(OBJDIR)
	$(CXX) $(CFLAGS) -o $@ $^ $(DFLAGS)

$(OBJDIR):
	mkdir $@

run: debug
	./gstar

clean:
	rm -rf obj gstar
