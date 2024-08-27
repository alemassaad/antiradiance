CXX = g++
CXXFLAGS = -Wall -std=c++17 -g
OBJS = main.o mesh.o radiosity.o file_io.o emissive_config.o
TARGET = antiradiance

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

main.o: main.cpp mesh.h radiosity.h file_io.h emissive_config.h
	$(CXX) $(CXXFLAGS) -c main.cpp

mesh.o: mesh.cpp mesh.h
	$(CXX) $(CXXFLAGS) -c mesh.cpp

radiosity.o: radiosity.cpp radiosity.h mesh.h
	$(CXX) $(CXXFLAGS) -c radiosity.cpp

file_io.o: file_io.cpp file_io.h mesh.h
	$(CXX) $(CXXFLAGS) -c file_io.cpp

emissive_config.o: emissive_config.cpp emissive_config.h mesh.h
	$(CXX) $(CXXFLAGS) -c emissive_config.cpp

clean:
	rm -f $(TARGET) $(OBJS)
	rm -f output.off output_checker
	rm -f objects/*.off
	rm -f ./timing_results.log
	rm -f ./form_factor_debug.txt
	