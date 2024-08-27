#include "mesh.h"
#include <iostream>
#include <limits>

enum EmissiveConfig {
    LEGO_BRICK,
    SMALL_ROOM, 
    BEAR
};

void setEmissiveTriangles(Mesh &mesh, EmissiveConfig config) {
    int emissiveCount = 0;
    int nonEmissiveCount = 0;

    if (config==LEGO_BRICK) {
        for (size_t i=0; i<mesh.indices.size(); ++i) {
            const auto& triangle = mesh.indices[i];

            double y1 = mesh.vertices[triangle.vtxi].y;
            double y2 = mesh.vertices[triangle.vtxj].y;
            double y3 = mesh.vertices[triangle.vtxk].y;

            if (y1>1.41f && y2>1.41f && y3>1.41f) {
                mesh.indices[i].emissivity = 1.0f;
                mesh.vertices[triangle.vtxi] = Vertex(1.0f, 1.0f, 0.0f); //R+G=Y
                mesh.vertices[triangle.vtxj] = Vertex(1.0f, 1.0f, 0.0f); 
                mesh.vertices[triangle.vtxk] = Vertex(1.0f, 1.0f, 0.0f); 
                emissiveCount++;
            } else {
                mesh.indices[i].emissivity = 0.0f;
                nonEmissiveCount++;
            }
        }
    } else if (config==SMALL_ROOM) {
        double max_y = -std::numeric_limits<double>::infinity();

        for (const auto& vertex : mesh.vertices) {
            if (vertex.y>max_y) {
                max_y = vertex.y;
            }
        }

        double threshold = max_y-1.0f; // ceiling height - light buffer

        for (size_t i=0; i<mesh.indices.size(); ++i) {
            const auto& triangle = mesh.indices[i];

            double y1 = mesh.vertices[triangle.vtxi].y;
            double y2 = mesh.vertices[triangle.vtxj].y;
            double y3 = mesh.vertices[triangle.vtxk].y;

            if (y1>threshold && y2>threshold && y3>threshold) {
                mesh.indices[i].emissivity = 1.0f;
                emissiveCount++;
            } else {
                mesh.indices[i].emissivity = 0.0f;
                nonEmissiveCount++;
            }
        }
    } else if (config==BEAR) {
        double midpoint_y = 10.36f/2; // bear's height=10.36

        for (size_t i=0; i<mesh.indices.size(); ++i) {
            const auto& triangle = mesh.indices[i];

            double y1 = mesh.vertices[triangle.vtxi].y;
            double y2 = mesh.vertices[triangle.vtxj].y;
            double y3 = mesh.vertices[triangle.vtxk].y;

            if (y1>midpoint_y && y2>midpoint_y && y3>midpoint_y) {
                mesh.indices[i].emissivity = 1.0f;
                mesh.vertices[triangle.vtxi] = Vertex(0.0f, 1.0f, 0.0f); // Green
                mesh.vertices[triangle.vtxj] = Vertex(0.0f, 1.0f, 0.0f); 
                mesh.vertices[triangle.vtxk] = Vertex(0.0f, 1.0f, 0.0f); 
                emissiveCount++;
            } else {
                mesh.indices[i].emissivity = 0.0f;
                nonEmissiveCount++;
            }
        }
    }

    std::cout << "Summary:  " << emissiveCount << "emissive triangles, " << nonEmissiveCount << "non-emissive triangles. " << std::endl;
}
