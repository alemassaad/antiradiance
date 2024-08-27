#pragma once
#include "mesh.h"
#include <vector>
#include <unordered_set>

class Radiosity {
public:
    Radiosity(Mesh& mesh);

    void computeRadiosityAntiradiance();

    const std::vector<std::vector<double>>& getFormFactors() const;
    const std::vector<Vertex>& getRadiance() const;
    const std::vector<Vertex>& getAntiradiance() const;

private:
    Mesh& mesh;
    std::vector<Vertex>radiance;
    std::vector<Vertex> antiradiance;
    std::vector<Vertex> combined;
    std::vector<std::vector<double>> formFactors;
    std::vector<std::vector<double>> reflectionMatrix;
    std::vector<std::vector<double>> goThroughMatrix;

    double computeTriangleArea(const Triangle& triangle) const;
    void validateFormFactors();
    void calculateFormFactors();
    void discretizeHemisphere();
    void initializeGoThroughMatrix(); 

    Vertex sampleTriangle(const Vertex& v1, const Vertex &v2, const Vertex& v3, int k, int l, int numSamples);
    Vertex computeNormal(const Vertex& v1, const Vertex& v2, const Vertex& v3);

    Vertex normalize(const Vertex& v);

    double lengthSquared(const Vertex& v) const;
    double safeDivision(double numerator, double denominator) const;
    double safeSqrt(double value) const;
    
    double dotProduct(const Vertex& v1, const Vertex& v2) const;

};
