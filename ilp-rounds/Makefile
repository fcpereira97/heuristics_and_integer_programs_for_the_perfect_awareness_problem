SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

#------------------------------------------------------------
#
# When you adapt this makefile to compile your CPLEX programs
# please copy this makefile and set CPLEXDIR and CONCERTDIR to
# the directories where CPLEX and CONCERT are installed.
#
#------------------------------------------------------------

CPLEXDIR = /home/felipe/ibm/ILOG/CPLEX_Studio1210/cplex
CONCERTDIR = /home/felipe/ibm/ILOG/CPLEX_Studio1210/concert

# ---------------------------------------------------------------------
# Compiler selection 
# ---------------------------------------------------------------------

#CCC = g++ -O0
CCC = g++ -O0 -Wall
CC  = gcc -O0
JAVAC = javac 

# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------

#CCOPT = -m64 -O -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD
CCOPT = -m64 -O -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD
COPT  = -m64 -fPIC -fno-strict-aliasing
JOPT  = -classpath $(CPLEXDIR)/lib/cplex.jar -O

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXJARDIR   = $(CPLEXDIR)/lib/cplex.jar
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

# For dynamic linking
CPLEXBINDIR   = $(CPLEXDIR)/bin/$(SYSTEM)
CPLEXLIB      = cplex$(dynamic:yes=12100)
run           = $(dynamic:yes=LD_LIBRARY_PATH=$(CPLEXBINDIR))

CCLNDIRS  = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR) $(dynamic:yes=-L$(CPLEXBINDIR))
CLNDIRS   = -L$(CPLEXLIBDIR) $(dynamic:yes=-L$(CPLEXBINDIR))
CCLNFLAGS = -lconcert -lilocplex -l$(CPLEXLIB) -lm -lpthread -ldl
CLNFLAGS  = -l$(CPLEXLIB) -lm -lpthread -ldl
JAVA      = java  -d64 -Djava.library.path=$(CPLEXDIR)/bin/x86-64_linux -classpath $(CPLEXJARDIR):

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

EXDIR         = $(CPLEXDIR)/examples
EXINC         = $(EXDIR)/include
EXDATA        = $(EXDIR)/data
EXSRCC        = $(EXDIR)/src/c
EXSRCCX       = $(EXDIR)/src/c_x
EXSRCCPP      = $(EXDIR)/src/cpp
EXSRCJAVA     = $(EXDIR)/src/java

CFLAGS  = $(COPT)  -I$(CPLEXINCDIR)
CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) 
JCFLAGS = $(JOPT)


#compile = g++ -DIL_STD
#cflags  = -Wall -c -O3
#xcflags  = $(cflags) -isystem$(cplex)/include -isystem$(concert)/include
#lflags  = -L$(cplex)/lib/x86-64_linux/static_pic -L$(concert)/lib/x86-64_linux/static_pic -lconcert -lilocplex -lcplex -lpthread -ldl
HEADERS =  my_struct.h load_input.h warm_start.h preprocessing_by_collapse.h nodal_heuristic.h
SOURCES =  my_struct.cpp load_input.cpp warm_start.cpp preprocessing_by_collapse.cpp nodal_heuristic.cpp pli-pa.cpp 
OBJECTS =  my_struct.o load_input.o warm_start.o preprocessing_by_collapse.o nodal_heuristic.o pli-pa.o

pli-pa: $(OBJECTS)
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o pli-pa $(OBJECTS) $(CCLNFLAGS)
	rm -f *.o

pli-pa.o: $(HEADERS) $(SOURCES)
	$(CCC) -c $(CCFLAGS) $(SOURCES)

clean:
	rm -f pli-pa *.o
