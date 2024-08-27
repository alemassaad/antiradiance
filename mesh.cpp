#include "mesh.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <thread>
#include <mutex>

bool Mesh::load(const std::string&filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return false;
    }

    size_t numVertices = 0;
    size_t numTriangles = 0;

    // one pass to count the nb of vertexe and indices ONLY
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string prefix;
        iss >> prefix;

        if (prefix=="v") {
            numVertices++;
        } else if (prefix=="f") {
            numTriangles++;
        }
    }
    vertices.reserve(numVertices);
    indices.reserve(numTriangles);

    file.clear();
    file.seekg(0, std::ios::beg);
    // second pass to load  actual data
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string prefix;
        iss >> prefix;

        if (prefix=="v") {  
            Vertex vertex;
            iss >> vertex.x >> vertex.y >> vertex.z;
            vertices.push_back(vertex);
            std::cout << "Loaded vertex: " << vertex.x << "," << vertex.y << ", " << vertex.z << std::endl;
        } else if (prefix=="f") {  
            Triangle triangle; 
            std::string vertex1, vertex2, vertex3;
            iss >> vertex1 >> vertex2 >> vertex3;

            triangle.vtxi = std::stoi(vertex1.substr(0, vertex1.find('/')))-1;
            triangle.vtxj = std::stoi(vertex2.substr(0, vertex2.find('/')))-1;
            triangle.vtxk = std::stoi(vertex3.substr(0, vertex3.find('/')))-1;

            triangle.emissivity =0.0f;
            indices.push_back(triangle);
        }
    }

    if (vertices.empty() || indices.empty()) {
        std::cerr << "Error: Mesh file does not contain any vertices or faces." << std::endl;
        return false;
    }

    std::cout << "Mesh successfully loaded with " << vertices.size() << " vertices and " << indices.size() << " triangles." << std::endl;

    return true;
}

bool Mesh::save(const std::string& filename, const std::vector<Vertex>& colors) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return false;
    }

    file << "OFF\n";
    file << vertices.size() << " " << indices.size() <<" 0\n";

    for (const auto& vertex : vertices) {
        file << vertex.x << " " << vertex.y << " " << vertex.z << "\n";
    }

    for (size_t i=0; i<indices.size(); ++i) {
        const auto& triangle = indices[i];
        const auto& color = colors[i];

        file << "3 " << triangle.vtxi << " " << triangle.vtxj << " " << triangle.vtxk << " "
             << color.x << " " << color.y << " " << color.z << "\n";
    }

    file.close();
    std::cout << "Mesh successfully saved to " << filename << std::endl;

    return true;
}
