# Include provided Makefile.inc for ARPACK++.
include ./Makefile.inc

# Directories
SRC_DIR := ..
BASIC_DIR := ../../Basics
SPIN_DIR := ../../Spin_1D
HEADER_DIRS := .. ../../Basics ../../Spin_1D
TEST_DIR := ./_Tests
EXE_DIR := ./_Exes
OBJ_DIR := ../_Objs
BASIC_OBJ_DIR := ../../Basics/_Objs
SPIN_OBJ_DIR := ../../Spin_1D/_Objs

# defining cscmat directory.
CSCMAT_DIR = ./matrices/complex
INC_DIR = ./include

# compiling and linking all examples.
SRC_FILES := $(foreach dir,$(SRC_DIR),$(wildcard $(dir)/*.cpp))
HEADER_FILES := $(foreach dir,$(HEADER_DIRS),$(wildcard $(dir)/*.hpp))
TEST_FILES := $(wildcard $(TEST_DIR)/*.cpp)

# Object files
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o,$(filter-out $(TEST_DIR)/%.cpp,$(SRC_FILES)))
BASIC_OBJ_FILES := $(wildcard $(BASIC_OBJ_DIR)/*.o)
SPIN_OBJ_FILES := $(wildcard $(SPIN_OBJ_DIR)/*.o)

# Executables
EXE_FILES := $(patsubst $(TEST_DIR)/%.cpp, $(EXE_DIR)/%.exe, $(TEST_FILES))

# Targets
.PHONY: all clean

all: $(EXE_FILES) $(OBJ_FILES) clean

$(EXE_DIR)/%.exe: $(TEST_DIR)/%.cpp $(OBJ_FILES) $(BASIC_OBJ_FILES) $(SPIN_OBJ_FILES)
	@mkdir -p $(EXE_DIR)
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) $^ -o $@ $(SUPERLU_LIB) $(ALL_LIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADER_FILES)
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -I$(INC_DIR) -I$(SUPERLU_DIR) -c $< -o $@

clean:
	rm -rf $(SRC_DIR)/*~ $(TEST_DIR)/*~

