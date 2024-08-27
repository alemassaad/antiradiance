#include "file_io.h"
#include <fstream>
#include <iostream>
#include <sstream>

bool loadObjFile(const std::string &filename, Mesh& mesh) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return false;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string prefix;
        iss >> prefix;

        if (prefix=="v") {
            Vertex v;
            iss >> v.x >> v.y >> v.z;
            mesh.vertices.push_back(v);
        } else if (prefix=="f") {
            std::string vertex1, vertex2, vertex3;
            iss >> vertex1 >> vertex2 >> vertex3;
            Triangle t;
            try {
                t.vtxi = std::stoi(vertex1)-1;
                t.vtxj = std::stoi(vertex2)-1;
                t.vtxk = std::stoi(vertex3)-1;
                if (t.vtxi >= mesh.vertices.size() || t.vtxj >= mesh.vertices.size() || t.vtxk >= mesh.vertices.size()) {
                    std::cerr << "Invalid face index in OBJ file: " << t.vtxi << " " << t.vtxj << " " << t.vtxk << std::endl;
                    continue;
                }

                //computing normal vector for the triangle
                const Vertex& v1 = mesh.vertices[t.vtxi];
                const Vertex& v2 = mesh.vertices[t.vtxj];
                const Vertex& v3 = mesh.vertices[t.vtxk];
                Vertex normal = {
                    (v2.y-v1.y) * (v3.z-v1.z) - (v2.z-v1.z) * (v3.y-v1.y),
                    (v2.z-v1.z) * (v3.x-v1.x) - (v2.x-v1.x) * (v3.z-v1.z),
                    (v2.x-v1.x) * (v3.y-v1.y) - (v2.y-v1.y) * (v3.x-v1.x)
                };

                double normalLengthSquared = normal.x*normal.x + normal.y*normal.y + normal.z*normal.z;

                if (normalLengthSquared < 1e-10) {
                    std::cerr << "Degenerate or nearly collinear triangle found and skipped:" << t.vtxi << " " << t.vtxj << " " << t.vtxk << std::endl;
                    continue;
                }

                // filtering out invalid triangles (w small area)

                mesh.indices.push_back(t);

            } catch (const std::exception& e) {
                std::cerr << "Error parsing face indices in OBJ file: " << e.what() << std::endl;
            }
        }
    }

    file.close();
    return true;
}

void saveOffFile(const std::string& filename, const Mesh & mesh, const std::vector<Vertex>& colors) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    file << "OFF\n";
    file << mesh.vertices.size() << " " << mesh.indices.size() << " 0\n";

    for (const auto& v : mesh.vertices) {
        file << v.x << " " << v.y << " " << v.z << "\n";
    }

    for (size_t i = 0; i < mesh.indices.size(); ++i) {
        const auto& t = mesh.indices[i];
        const auto& color = colors[i];

        file << "3 " << t.vtxi << " " << t.vtxj << " " << t.vtxk << " "
             << std::clamp(color.x, 0.0, 1.0) << " "
             << std::clamp(color.y, 0.0, 1.0) << " "
             << std::clamp(color.z, 0.0, 1.0) << "\n";
    }

    file.close();

    std::ofstream logFile("mesh_save_log.txt");
    if (logFile.is_open()) {
        logFile << "Mesh saved to " << filename << "\n";
        logFile << "Vertices: " << mesh.vertices.size() << "\n";
        logFile << "Indices: " << mesh.indices.size() << "\n";
        
        for (size_t i=0; i<mesh.indices.size(); ++i) {
            const auto& t = mesh.indices[i];
            const auto& color = colors[i];
            logFile << "Triangle " << i << ": Indices (" << t.vtxi << ", " << t.vtxj << ", " << t.vtxk << "), Color: ("
                    << color.x << ", " << color.y << ", " << color.z << ")\n";
        }
    }
    logFile.close();

    std::cout << "Mesh successfully saved to " << filename << std::endl;
}

void validateOffFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    std::string line;
    std::getline(file, line);
    if (line != "OFF") {
        std::cerr << "Invalid OFF header" << std::endl;
        return;
    }

    std::getline(file, line);
    std::istringstream header(line);
    int numVertices, numFaces, numEdges;
    header >> numVertices >> numFaces >> numEdges;
    if (numVertices <= 0 || numFaces <= 0 || numEdges!=0) {
        std::cerr << "Invalid OFF header values" << std::endl;
        return;
    }

    int vertexCount = 0;
    int faceCount = 0;
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        if (vertexCount<numVertices) {
            double x, y, z;
            ss >> x >> y >> z;
            if (ss.fail()) {
                std::cerr << "Invalid vertex format at line " << vertexCount+3 << std::endl;
                return;
            }
            vertexCount++;
        } else if (faceCount<numFaces) {
            int v1, v2, v3;
            double r, g, b;
            std::string prefix;
            ss >> prefix >> v1 >> v2 >> v3 >> r >> g>> b;
            if (prefix!="3" || ss.fail() || r<0.0f || r>1.0f || g<0.0f || g>1.0f || b<0.0f || b>1.0f) {
                std::cerr << "Invalid face or color format at line " << vertexCount+faceCount+3 << std::endl;
                return;
            }
            faceCount++;
        } else {
            std::cerr << "Extra data after expected vertices and faces" << std::endl;
            return;
        }
    }

    if (vertexCount!=numVertices) {
        std::cerr << "Vertex count mismatch: expected " << numVertices << " but found " << vertexCount << std::endl;
        return;
    }

    if (faceCount!=numFaces) {
        std::cerr << "Face count mismatch: expected " << numFaces << " but found " << faceCount << std::endl;
        return;
    }

    std::cout << "OFF file validation successful" << std::endl;
    file.close();
}
