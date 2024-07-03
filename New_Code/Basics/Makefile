# Makefile

# Directories
SRC_DIR := .
TEST_DIR := ./_Tests
EXE_DIR := ./_Exes
OBJ_DIR := ./_Objs

# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++20 -Wall -Wextra -pedantic
LDFLAGS := #-lblas -llapack -lgfortran

# Source files
SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
HEADER_FILES := $(wildcard $(SRC_DIR)/*.hpp)
TEST_FILES := $(wildcard $(TEST_DIR)/*.cpp)

# Object files
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRC_FILES))

# Executables
EXE_FILES := $(patsubst $(TEST_DIR)/%.cpp, $(EXE_DIR)/%.exe, $(TEST_FILES))

# Targets
all: $(EXE_FILES) $(OBJ_FILES) clean

$(EXE_DIR)/%.exe: $(TEST_DIR)/%.cpp $(OBJ_FILES)
	@mkdir -p $(EXE_DIR)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADER_FILES)
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: all clean

clean:
	rm -rf $(SRC_DIR)/*.cpp~ $(SRC_DIR)/*.hpp~ $(TEST_DIR)/*.cpp~
