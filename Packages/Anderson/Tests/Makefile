# ARPACK++ v1.2 2/20/2000
# c++ interface to ARPACK code.
# examples/superlu/complex directory makefile.

# including other makefiles.

include ./Makefile.inc

# defining cscmat directory.

CSCMAT_DIR = ./matrices/complex
INC_DIR = ../include

# compiling and linking all examples.

all: unitTest XXZTest #lcompshf lcompgre lcompgsh

# compiling and linking each complex problem.

unitTest:	unitTest.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o unitTest unitTest.o $(SUPERLU_LIB) $(ALL_LIBS)

XXZTest:	XXZTest.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o XXZTest XXZTest.o $(SUPERLU_LIB) $(ALL_LIBS)

# defining cleaning rule.

.PHONY:	clean
clean:
	rm -f *~ *.o core unitTest XXZTest

# defining pattern rules.

%.o:	%.cpp
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -I$(INC_DIR) -I$(SUPERLU_DIR) -c $<

