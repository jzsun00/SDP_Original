# Compiler
CC = g++

# Compiler flags
CFLAGS = -std=c++11 -Wall

# Source directory
SRC_DIR = .

# Include directory
INC_DIR = .

# Build directory
BUILD_DIR = .

# Executable directory
BIN_DIR = .

# Source files
SOURCES = $(wildcard $(SRC_DIR)/*.cpp)

# Header files
HEADERS = $(wildcard $(INC_DIR)/*.hpp)

# Executable names
EXECUTABLES = $(patsubst $(SRC_DIR)/%.cpp,$(BIN_DIR)/%,$(SOURCES))

# Default target
all: $(EXECUTABLES)

# Rule to build each executable
$(BIN_DIR)/%: $(BUILD_DIR)/%.o $(HEADERS)
	$(CC) $(CFLAGS) $< -o $@

# Rule to compile .cpp files into .o files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS) | $(BUILD_DIR)
	$(CC) $(CFLAGS) -I$(INC_DIR) -c $< -o $@

# Create build directory
#$(BUILD_DIR):
#    mkdir -p $(BUILD_DIR)

# Clean rule
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

.PHONY: all clean

