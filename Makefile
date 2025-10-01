# Makefile for Signal Analysis Project

# Compiler and flags
CC = gcc
CFLAGS = -Wall -O2 -Isrc/
LDFLAGS = -lfftw3 -lm

# Directories
SRC_DIR = src
BUILD_DIR = build
TARGET = $(BUILD_DIR)/analyzer

# Source files and object files
SOURCES = $(wildcard $(SRC_DIR)/*.c)
OBJECTS = $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o, $(SOURCES))

# Default target
all: $(TARGET)

# Linking the executable
$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CC) $(OBJECTS) -o $(TARGET) $(LDFLAGS)
	@echo "Linking complete. Executable is at $(TARGET)"

# Compiling object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@

# Clean target
clean:
	@echo "Cleaning up build files..."
	rm -rf $(BUILD_DIR)
	@echo "Cleanup complete."

# Phony targets
.PHONY: all clean