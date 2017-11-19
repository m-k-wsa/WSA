SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

#------------------------------------------------------------
#
# When you adapt this makefile to compile your CPLEX programs
# please copy this makefile and set CPLEXDIR and CONCERTDIR to
# the directories where CPLEX and CONCERT are installed.
#
#------------------------------------------------------------

CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio1271/cplex
CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio1271/concert
# ---------------------------------------------------------------------
# Compiler selection 
# ---------------------------------------------------------------------

CCC = g++ -O0 -std=c++0x
CC  = gcc -O0
JAVAC = javac 

# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------

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

CCLNDIRS  = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR)
CLNDIRS   = -L$(CPLEXLIBDIR)
CCLNFLAGS = -lconcert -lilocplex -lcplex -lm -lpthread 
CLNFLAGS  = -lcplex -lm -lpthread
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

CFLAGS  = -std=c++0x $(COPT)  -I$(CPLEXINCDIR)
CCFLAGS = -std=c++0x $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) 
JCFLAGS = $(JOPT)


#------------------------------------------------------------
#  make all      : to compile the examples. 
#  make execute  : to compile and execute the examples.
#------------------------------------------------------------





# ------------------------------------------------------------

clean :
	/bin/rm -rf *.o *~ *.class
	/bin/rm -rf $(C_EX) $(CX_EX) $(CPP_EX)
	/bin/rm -rf *.mps *.ord *.sos *.lp *.sav *.net *.msg *.log *.clp
	find . -type f | xargs -n 5 touch
	rm -rf $(OBJS)

# ------------------------------------------------------------
#
# The examples
#

all : wsa

wsa: wsa.o task.o wa.o wa_kill.o exp.o util.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o wsa task.o wa.o wa_kill.o exp.o wsa.o util.o $(CCLNFLAGS)

wsa.o: task.hpp wa.hpp 
	$(CCC) -c $(CCFLAGS) wsa.cpp -o wsa.o


task.o : task.cpp util.hpp
	$(CCC) -c $(CCFLAGS) task.cpp -o task.o

wa.o : wa.hpp
	$(CCC) -c $(CCFLAGS) wa.cpp -o wa.o

wa_kill.o: wa.hpp 
	$(CCC) -c $(CCFLAGS) wa_kill.cpp -o wa_kill.o

exp.o : exp.hpp util.hpp
	$(CCC) -c $(CCFLAGS) exp.cpp -o exp.o

util.o : util.hpp
	$(CCC) -c $(CCFLAGS) util.cpp -o util.o



#blend: blend.o
#	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o blend blend.o $(CCLNFLAGS)
#blend.o: $(EXSRCCPP)/blend.cpp
#	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/blend.cpp -o blend.o



# Local Variables:
# mode: makefile
# End:
