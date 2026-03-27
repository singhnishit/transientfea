CXX      = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra -Iinclude
SRC      = src/linalg.cpp src/beam.cpp src/forces.cpp src/solver.cpp src/output.cpp src/main.cpp
OBJ      = $(SRC:.cpp=.o)
TARGET   = beam_fea

.PHONY: all clean run visualize

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^

run: all
	mkdir -p results
	./$(TARGET)

visualize: run
	python3 visualize.py

clean:
	rm -f $(TARGET) src/*.o results/*.csv results/*.gif results/*.png
