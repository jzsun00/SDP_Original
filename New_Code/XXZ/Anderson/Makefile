# ARPACK++ v1.2 2/20/2000
# c++ interface to ARPACK code.
# examples/superlu/complex directory makefile.

# including other makefiles.

include ./Makefile.inc

SRC_DIR := ../../Basics ../../Spin_1D ..
TEST_DIR := ./_Tests
EXE_DIR := ./_Exes

# defining cscmat directory.

CSCMAT_DIR = ./matrices/complex
INC_DIR = ./include

# compiling and linking all examples.
SRC_FILES := $(foreach dir,$(SRC_DIR),$(wildcard $(dir)/*.cpp))
HEADER_FILES := $(foreach dir,$(SRC_DIR),$(wildcard $(dir)/*.hpp))
TEST_FILES := $(wildcard $(TEST_DIR)/*.cpp)
OBJ_FILES := $(patsubst %.cpp,%.o,$(filter-out $(TEST_DIR)/%.cpp,$(SRC_FILES)))
EXE_FILES := $(patsubst $(TEST_DIR)/%.cpp, $(EXE_DIR)/%, $(TEST_FILES))

#all: XXZTest #lcompshf lcompgre lcompgsh

# compiling and linking each complex problem.

#XXZTest:	XXZTest.o
#	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o XXZTest XXZTest.o $(SUPERLU_LIB) $(ALL_LIBS)

# defining cleaning rule.

.PHONY: all clean

all: $(EXE_FILES)

$(EXE_DIR)/%: $(TEST_DIR)/%.cpp $(OBJ_FILES)
	@mkdir -p $(EXE_DIR)
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) $^ -o $@ $(SUPERLU_LIB) $(ALL_LIBS)

# defining pattern rules.

%.o:	%.cpp $(HEADER_FILES)
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -I$(INC_DIR) -I$(SUPERLU_DIR) -c $< -o $@

clean:
	rm -rf $(EXE_DIR) $(SRC_DIR)/*.o

