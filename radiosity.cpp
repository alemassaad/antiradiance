#include "radiosity.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <limits>
#include <unordered_set>
#include <vector>
#include <fstream>

Radiosity::Radiosity(Mesh& mesh) 
    : mesh(mesh), 
      radiance(mesh.indices.size()), 
      antiradiance(mesh.indices.size()),
      formFactors(mesh.indices.size(), std::vector<double>(mesh.indices.size(), 0.0)),
      reflectionMatrix(mesh.indices.size(), std::vector<double>(mesh.indices.size(), 0.0)),
      goThroughMatrix(mesh.indices.size(), std::vector<double>(mesh.indices.size(), 0.0)) 
{
    std::cout << "Starting form factor calculations..." << std::endl;
    calculateFormFactors();
    discretizeHemisphere();
    initializeGoThroughMatrix(); 
    std::cout << "Form factors, hemisphere discretization, identity matrix, and go-through matrix initialization completed." << std::endl;
}

double Radiosity::computeTriangleArea(const Triangle& triangle) const {
    const Vertex& v1 = mesh.vertices[triangle.vtxi];
    const Vertex& v2 = mesh.vertices[triangle.vtxj];
    const Vertex& v3 = mesh.vertices[triangle.vtxk];

    Vertex edge1 = {v2.x-v1.x, v2.y-v1.y, v2.z-v1.z};
    Vertex edge2 = {v3.x-v1.x, v3.y-v1.y, v3.z-v1.z};

    Vertex crossProduct = {
        edge1.y*edge2.z - edge1.z*edge2.y,
        edge1.z*edge2.x - edge1.x*edge2.z,
        edge1.x*edge2.y - edge1.y*edge2.x
    };

    double area = safeSqrt(lengthSquared(crossProduct)) * 0.5;

    return area;
}

void Radiosity::validateFormFactors() {
    size_t numTriangles = mesh.indices.size();
    std::ofstream debugFile("form_factor_debug.txt");

    for (size_t i=0; i<numTriangles; ++i) {
        double sumF = 0.0;
        for (size_t j=0; j<numTriangles; ++j) {
            sumF += formFactors[i][j];

            double area_i = computeTriangleArea(mesh.indices[i]);
            double area_j = computeTriangleArea(mesh.indices[j]);
            double reciprocalDifference = std::abs(area_i*formFactors[i][j] - area_j*formFactors[j][i]);

            if (reciprocalDifference > 1e-2) {
                std::cerr << "Reciprocity checkfailed between triangles" << i << " and " << j << " with error" << reciprocalDifference << std::endl;
                debugFile << "Reciprocity checkfailed between triangles" << i << " and " << j << " with error" << reciprocalDifference << std::endl;
                debugFile << "Triangle "<< i << ": Area = " << area_i << ", Form Factor = " << formFactors[i][j] << std::endl;
                debugFile << "Triangle " << j << ": Area = " << area_j << ", Form Factor = " << formFactors[j][i] << std::endl;
                debugFile << "Reciprocal difference: " << reciprocalDifference << std::endl;
            }
        }
        if (sumF > 1.0 + 1e-2) {
            std::cerr << "Energy conservation violated for triangle " << i << ": sum of form factors = " << sumF << std::endl;
            debugFile << "Energy conservation violated for triangle " << i << ": sum of form factors = " << sumF << std::endl;
        }
    }

    debugFile.close();
}

void Radiosity::calculateFormFactors() {
    size_t numTriangles = mesh.indices.size();
    const int numSamples = 20; 
    const double omegaBin = 4.0*M_PI /(numSamples*numSamples); // solid angle of bin

    for (size_t i = 0; i < numTriangles; ++i) {
        const Vertex& v1 = mesh.vertices[mesh.indices[i].vtxi];
        const Vertex& v2 = mesh.vertices[mesh.indices[i].vtxj];
        const Vertex& v3 = mesh.vertices[mesh.indices[i].vtxk];

        Vertex normal1 = computeNormal(v1, v2, v3);
        if (lengthSquared(normal1)<1e-10) {
            std::cerr << "Skipping degenerate triangle in form factor calculation: Triangle " << i << std::endl;
            continue;
        }
        normal1 = normalize(normal1);

        for (size_t j=0; j<numTriangles; ++j) {
            if (i!=j) {
                double formFactor = 0.0;

                const Vertex& u1 = mesh.vertices[mesh.indices[j].vtxi];
                const Vertex& u2 = mesh.vertices[mesh.indices[j].vtxj];
                const Vertex& u3 = mesh.vertices[mesh.indices[j].vtxk];

                Vertex normal2= computeNormal(u1, u2, u3);
                if (lengthSquared(normal2)<1e-10) { //degenerate
                    continue; 
                }
                normal2 = normalize(normal2);

                double area_j = computeTriangleArea(mesh.indices[j]);

                for (int k=0; k<numSamples; ++k) {
                    for (int l=0; l<numSamples; ++l) {
                        Vertex p1 = sampleTriangle(v1, v2, v3, k, l, numSamples);
                        Vertex p2 = sampleTriangle(u1, u2, u3, l, k, numSamples);

                        Vertex direction = { p1.x-p2.x, p1.y-p2.y, p1.z-p2.z };
                        double distanceSquared = lengthSquared(direction);

                        if (distanceSquared<1e-20) {
                            continue;
                        }

                        double distance = safeSqrt(distanceSquared);
                        direction.x /= distance;
                        direction.y /= distance;
                        direction.z /= distance;

                        double cosTheta1 = std::max(0.0, dotProduct(normal1, direction));
                        double cosTheta2 = std::max(0.0, -dotProduct(normal2, direction));

                        double dFormFactor = cosTheta1*cosTheta2*area_j /(distanceSquared*omegaBin);
                        formFactor += dFormFactor;
                    }
                }

                formFactor /= (numSamples*numSamples);

                formFactors[i][j] = formFactor;
                reflectionMatrix[i][j] = formFactor;
                goThroughMatrix[i][j] = 0.0;
            } else {
                formFactors[i][j] = 0.0;
                reflectionMatrix[i][j] = 0.0;
                goThroughMatrix[i][j] = 1.0;
            }
        }
        if (i%100 ==0 || i == numTriangles-1) {
            std::cout << "Processed form factors for " << i << " out of " << numTriangles << " triangles" << std::endl;
        }
    }

    validateFormFactors();
}

void Radiosity::discretizeHemisphere() {
    size_t numTriangles = mesh.indices.size();
    const int numThetaSamples = 20;
    const int numPhiSamples = 40;

    std::unordered_set<size_t> processedTriangles;

    for (size_t i=0; i<numTriangles; ++i) {
        if (processedTriangles.find(i)!=processedTriangles.end()) continue;

        std::vector<Vertex> hemisphereDirections;

        Vertex normal = computeNormal(
            mesh.vertices[mesh.indices[i].vtxi],
            mesh.vertices[mesh.indices[i].vtxj],
            mesh.vertices[mesh.indices[i].vtxk]
        );

        if (lengthSquared(normal)<1e-10) {
            std::cerr << "Warning: Triangle " << i << " has a zero-length normal vector!" << std::endl;
            continue;  // skip if normal is v v small
        }

        normal = normalize(normal);

        for (int t=0; t<numThetaSamples; ++t) {
            double theta = (t+0.5)*M_PI /(2*numThetaSamples);
            for (int p=0; p<numPhiSamples; ++p) {
                double phi = (p+0.5)*2*M_PI /numPhiSamples;

                Vertex direction = {
                    static_cast<float>(std::sin(theta) * std::cos(phi)),
                    static_cast<float>(std::sin(theta) * std::sin(phi)),
                    static_cast<float>(std::cos(theta))
                };

                double dotProduct = direction.x*normal.x + direction.y*normal.y + direction.z*normal.z;
                if (dotProduct<0) {
                    direction.x = -direction.x;
                    direction.y = -direction.y;
                    direction.z = -direction.z;
                }

                hemisphereDirections.push_back(direction);
            }
        }

        std::cout << "Discretized " << hemisphereDirections.size() << " hemisphere directions for triangle " << i << std::endl;
        processedTriangles.insert(i);
    }

    std::cout << "Hemisphere discretization completed for  all triangles. " << std::endl;
}

void Radiosity::initializeGoThroughMatrix() {
    size_t numTriangles = mesh.indices.size();
    for (size_t i=0; i<numTriangles; ++i) {
        for (size_t j=0; j<numTriangles; ++j) {
            if (i==j) {
                goThroughMatrix[i][j] = 1.0;
            } else {
                goThroughMatrix[i][j] = 0.0;
            }
        }
    }

    std::cout << "Go-through matrix J initialized as an identity matrix." << std::endl;
}

void Radiosity::computeRadiosityAntiradiance() {
    size_t numTriangles = mesh.indices.size();
    std::vector<Vertex> E(numTriangles);

    for (size_t i=0; i<numTriangles; ++i) {
        double emissivity = mesh.indices[i].emissivity;
        E[i] = {static_cast<float>(emissivity), static_cast<float>(emissivity), static_cast<float>(emissivity)};
    }

    std::vector<Vertex> L = E;
    std::vector<Vertex> A(numTriangles, {0.0f, 0.0f, 0.0f});
    const int maxIterations = 1000;
    const double convergenceThreshold = 1e-10;

    for (int iter=0; iter<maxIterations; ++iter) {
        std::vector<Vertex> newL(numTriangles, {0.0f, 0.0f, 0.0f});
        std::vector<Vertex> newA(numTriangles, {0.0f, 0.0f, 0.0f});
        bool converged = true;

        for (size_t i=0; i<numTriangles; ++i) {
            Vertex L_new = E[i];

            for (size_t j=0; j<numTriangles; ++j) {
                if (i!=j) {
                    double reflectionTerm = reflectionMatrix[i][j];
                    L_new.x += static_cast<float>(reflectionTerm * (L[j].x-A[j].x));
                    L_new.y += static_cast<float>(reflectionTerm * (L[j].y-A[j].y));
                    L_new.z += static_cast<float>(reflectionTerm * (L[j].z-A[j].z));
                }
            }

            newL[i] = L_new;

            if (std::abs(newL[i].x-L[i].x) > convergenceThreshold ||
                std::abs(newL[i].y-L[i].y) > convergenceThreshold ||
                std::abs(newL[i].z-L[i].z) > convergenceThreshold) {
                converged=false;
            }
        }

        for (size_t i=0; i<numTriangles; ++i) {
            Vertex A_new = {0.0f, 0.0f, 0.0f};

            for (size_t j=0; j<numTriangles; ++j) {
                if (i!=j) {
                    double goThroughTerm = goThroughMatrix[i][j];
                    A_new.x += static_cast<float>(goThroughTerm * (newL[j].x-A[j].x));
                    A_new.y += static_cast<float>(goThroughTerm * (newL[j].y-A[j].y));
                    A_new.z += static_cast<float>(goThroughTerm * (newL[j].z-A[j].z));
                }
            }

            newA[i] = A_new;
        }

        L=newL;
        A=newA; 

        if (converged) {break;}

    }
    for (size_t i=0; i<L.size(); ++i) {
        L[i].x = std::min(1.0, std::max(0.0, static_cast<double>(L[i].x)));
        L[i].y = std::min(1.0, std::max(0.0, static_cast<double>(L[i].y)));
        L[i].z = std::min(1.0, std::max(0.0, static_cast<double>(L[i].z)));
    }
    radiance = L;
    antiradiance = A;
}

const std::vector<std::vector<double>>& Radiosity::getFormFactors() const {
    return formFactors;
}
const std::vector<Vertex>& Radiosity::getRadiance() const {
    return radiance;
}
const std::vector<Vertex>&Radiosity::getAntiradiance() const {
    return antiradiance;
}

Vertex Radiosity::sampleTriangle(const Vertex& v1, const Vertex &v2, const Vertex& v3, int k, int l, int numSamples) {
    double u = static_cast<double>(k)/numSamples;
    double v = static_cast<double>(l)/numSamples;
    if (u+v > 1.0) {
        u = 1.0-u;
        v = 1.0-v;
    }
    return {
        static_cast<float>((1-u-v)*v1.x + u*v2.x + v*v3.x),
        static_cast<float>((1-u-v)*v1.y + u*v2.y + v*v3.y),
        static_cast<float>((1-u-v)*v1.z + u*v2.z + v*v3.z)
    };
}

Vertex Radiosity::computeNormal(const Vertex& v1, const Vertex& v2, const Vertex & v3) {
    Vertex edge1 = { v2.x-v1.x, v2.y-v1.y, v2.z-v1.z };
    Vertex edge2 = { v3.x-v1.x, v3.y-v1.y, v3.z-v1.z };

    Vertex normal = {
        edge1.y*edge2.z - edge1.z*edge2.y,
        edge1.z*edge2.x - edge1.x*edge2.z,
        edge1.x*edge2.y - edge1.y*edge2.x
    };

    double length = safeSqrt(lengthSquared(normal));
    if (length>1e-10) {
        normal.x = static_cast<float>(normal.x/length);
        normal.y = static_cast<float>(normal.y/length);
        normal.z = static_cast<float>(normal.z/length);
    } else {
        normal = {0.0f, 0.0f, 0.0f};  // for zero-length normal vectors
    }
    return normal;
}

Vertex Radiosity::normalize(const Vertex& v) {
    double length = std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
    if (length>1e-10) {  //no dividing by v small numbers close to 0
        return {static_cast<float>(v.x/length), static_cast<float>(v.y/length), static_cast<float>(v.z/length)};
    } else {
        return {0.0f, 0.0f, 0.0f};  
    }
}

double Radiosity::lengthSquared(const Vertex& v) const {
    return dotProduct(v, v);}

double Radiosity::safeDivision(double numerator, double denominator) const {
    if (std::abs(denominator)>1e-10) {  
        return numerator/denominator;
    } else {
        return 0.0; 
    }
}

double Radiosity::safeSqrt(double value) const {return std::sqrt(std::max(value,0.0));}

double Radiosity::dotProduct(const Vertex & v1, const Vertex& v2) const {
    return v1.x*v2.x + v1.y*v2.y +v1.z*v2.z;
}
