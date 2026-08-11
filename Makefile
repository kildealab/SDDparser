CXX = g++

CXXFLAGS = -std=c++17 -Wall -Wextra -O2 -Iinclude

TARGET = SDDparser

SRC_DIR = source
OBJ_DIR = objects

SOURCES = $(SRC_DIR)/main.cpp \
          $(SRC_DIR)/SDDparser.cpp

OBJECTS = $(OBJ_DIR)/main.o \
          $(OBJ_DIR)/SDDparser.o


$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECTS)


$(OBJ_DIR)/main.o: $(SRC_DIR)/main.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@


$(OBJ_DIR)/SDDparser.o: $(SRC_DIR)/SDDparser.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@


clean:
	rm -rf $(OBJ_DIR)/*.o $(TARGET)
