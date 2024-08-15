# Makefile

# Directories
SRC_DIR := .
BASIC_DIR := ../Basics
SPIN_DIR := ../Spin_1D
FERMI_DIR := ../Fermion
HARDCORE_DIR := ../HardCore
LA_DIR := ../../LinearAlgebra
HEADER_DIRS = . ../Basics ../Spin_1D ../Fermion ../HardCore ../../LinearAlgebra
TEST_DIR := ./_Tests
EXE_DIR := ./_Exes
OBJ_DIR := ./_Objs
BASIC_OBJ_DIR := ../Basics/_Objs
SPIN_OBJ_DIR := ../Spin_1D/_Objs
FERMI_OBJ_DIR := ../Fermion/_Objs
HARDCORE_OBJ_DIR := ../HardCore/_Objs
LA_OBJ_DIR := ../LinearAlgebra/_Objs

# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++20 -Wall -Wextra -pedantic -fopenmp
LDFLAGS := -lsuperlu  -larpack -lopenblas -llapacke -lgfortran

# Source files
SRC_FILES := $(foreach dir,$(SRC_DIR),$(wildcard $(dir)/*.cpp))
HEADER_FILES := $(foreach dir,$(HEADER_DIRS),$(wildcard $(dir)/*.hpp))
TEST_FILES := $(wildcard $(TEST_DIR)/*.cpp)

# Object files
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o,$(filter-out $(TEST_DIR)/%.cpp,$(SRC_FILES)))
BASIC_OBJ_FILES := $(wildcard $(BASIC_OBJ_DIR)/*.o)
SPIN_OBJ_FILES := $(wildcard $(SPIN_OBJ_DIR)/*.o)
FERMI_OBJ_FILES := $(wildcard $(FERMI_OBJ_DIR)/*.o)
HARDCORE_OBJ_FILES := $(wildcard $(HARDCORE_OBJ_DIR)/*.o)
LA_OBJ_FILES := $(wildcard $(LA_OBJ_DIR)/*.o)

# Executables
EXE_FILES := $(patsubst $(TEST_DIR)/%.cpp, $(EXE_DIR)/%.exe, $(TEST_FILES))

# Targets
.PHONY: all clean

all: $(EXE_FILES) $(OBJ_FILES) clean

$(EXE_DIR)/%.exe: $(TEST_DIR)/%.cpp $(OBJ_FILES) $(BASIC_OBJ_FILES) $(SPIN_OBJ_FILES) $(FERMI_OBJ_FILES) $(HARDCORE_OBJ_FILES) $(LA_OBJ_FILES)
	@mkdir -p $(EXE_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADER_FILES)
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(TEST_DIR)/%~ $(SRC_DIR)/%~
