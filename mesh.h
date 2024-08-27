#pragma once
#include <vector>
#include <string>
#include <mutex>

struct Vertex {
    double x, y, z;
    Vertex(double x, double y, double z) : x(x), y(y), z(z) {}
    Vertex() : x(0), y(0), z(0) {}
};

struct Triangle {
    unsigned int vtxi, vtxj, vtxk;
    double emissivity; 
};

class Mesh {
public:
    std::vector<Vertex> vertices;
    std::vector<Triangle> indices;

    bool load(const std::string& filename);
    
    bool save(const std::string &filename, const std::vector<Vertex>& colors) const; 

};
