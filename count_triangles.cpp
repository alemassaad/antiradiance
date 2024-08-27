#include <iostream>
#include <string>
#include "mesh.h"      
#include "file_io.h"   

int main() {
    Mesh mesh;
    std::string filename = "room_radiosity.obj";

    if (loadObjFile(filename, mesh)) { 
            std::cout << "Total number of triangles: " << mesh.indices.size() << std::endl;
    } else {
            std::cerr << "Failed to load mesh from file: " << filename << std::endl;
    }

    return 0;
}
