CXXFLAGS=-std=c++11 -Wall -pedantic -pthread

all: rainfall_seq rainfall_pt

rainfall_seq: rainfall_seq.cpp
	g++ $(CXXFLAGS) -o rainfall_seq rainfall_seq.cpp

rainfall_pt: rainfall_pt.cpp
	g++ $(CXXFLAGS) -o rainfall_pt rainfall_pt.cpp

clean:
	rm rainfall_seq
	rm rainfall_pt
