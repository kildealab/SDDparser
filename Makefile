CXX = g++

CXXFLAGS = -std=c++17 -Wall -Wextra -O2 -Iinclude

TARGET = SDDparser

SRC_DIR = source
OBJ_DIR = objects
INC_DIR = include

# Recursively find all files matching a pattern under a directory,
# at any depth - so adding source/karyogram/ (or any future
# subdirectory) needs no Makefile changes.
rwildcard = $(wildcard $1$2) $(foreach d,$(wildcard $1*),$(call rwildcard,$d/,$2))


SOURCES = $(call rwildcard,$(SRC_DIR)/,*.cpp)

OBJECTS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SOURCES))

HEADERS = $(call rwildcard,$(INC_DIR)/,*.h)


$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECTS) -lcairo


$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@


clean:
	rm -rf $(OBJ_DIR)
	rm -f $(TARGET)
