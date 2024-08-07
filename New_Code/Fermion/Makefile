# Makefile

# Directories
SRC_DIR := .
BASIC_DIR := ../Basics
LA_DIR := ../LinearAlgebra
TEST_DIR := ./_Tests
EXE_DIR := ./_Exes
OBJ_DIR := ./_Objs
BASIC_OBJ_DIR := ../Basics/_Objs
LA_OBJ_DIR := ../LinearAlgebra/_Objs

# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++20 -Wall -Wextra -pedantic
LDFLAGS := -lopenblas -llapacke -lgfortran

# Source files
SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
HEADER_FILES := $(wildcard $(SRC_DIR)/*.hpp) $(wildcard $(BASIC_DIR)/*.hpp) $(wildcard $(LA_DIR)/*.hpp)
TEST_FILES := $(wildcard $(TEST_DIR)/*.cpp)

# Object files
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRC_FILES))
BASIC_OBJ_FILES := $(wildcard $(BASIC_OBJ_DIR)/*.o)
LA_OBJ_FILES := $(wildcard $(LA_OBJ_DIR)/*.o)

# Executables
EXE_FILES := $(patsubst $(TEST_DIR)/%.cpp, $(EXE_DIR)/%.exe, $(TEST_FILES))

# Targets
.PHONY: all clean

all: $(EXE_FILES) $(OBJ_FILES) clean

$(EXE_DIR)/%.exe: $(TEST_DIR)/%.cpp $(OBJ_FILES) $(BASIC_OBJ_FILES) $(LA_OBJ_FILES)
	@mkdir -p $(EXE_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADER_FILES)
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(TEST_DIR)/*~ $(SRC_DIR)/*~
