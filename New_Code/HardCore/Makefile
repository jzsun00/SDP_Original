# Makefile

# Directories
SRC_DIR := ../Basics .
TEST_DIR := ./_Tests
EXE_DIR := ./_Exes

# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++11 -Wall -Wextra
LDFLAGS := -lblas -llapack -lgfortran

# Source files
SRC_FILES := $(foreach dir,$(SRC_DIR),$(wildcard $(dir)/*.cpp))
HEADER_FILES := $(foreach dir,$(SRC_DIR),$(wildcard $(dir)/*.hpp))
#SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
#HEADER_FILES := $(wildcard $(SRC_DIR)/*.hpp)
TEST_FILES := $(wildcard $(TEST_DIR)/*.cpp)

# Object files
#OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp, $(SRC_DIR)/%.o, $(SRC_FILES))
OBJ_FILES := $(patsubst %.cpp,%.o,$(filter-out $(TEST_DIR)/%.cpp,$(SRC_FILES)))
# Executables
EXE_FILES := $(patsubst $(TEST_DIR)/%.cpp, $(EXE_DIR)/%, $(TEST_FILES))

# Targets
.PHONY: all clean

all: $(EXE_FILES)

$(EXE_DIR)/%: $(TEST_DIR)/%.cpp $(OBJ_FILES)
	@mkdir -p $(EXE_DIR)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@

#$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADER_FILES)
#	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.cpp $(HEADER_FILES)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(EXE_DIR) $(SRC_DIR)/*.o
