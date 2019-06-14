BINARY=projection
CXX=g++
CXXFLAGS = -O3 -Wall
DEPS = img_out.hpp auxiliary.hpp

LDFLAGS=

OBJ=main.o img_out.o

default: $(BINARY)


%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)


$(BINARY): $(OBJ)
	$(CXX) -o $@ $^ $(LDFLAGS)


.PHONY: clean clean_all

clean:
	rm -f *.o $(BINARY)
