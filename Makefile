CXXFLAGS=-std=c++11 -Wall -pedantic

rainfall: rainfall.cpp
	g++ $(CXXFLAGS) -o rainfall rainfall.cpp

clean:
	rm rainfall
