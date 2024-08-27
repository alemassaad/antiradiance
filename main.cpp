#include "mesh.h"
#include "radiosity.h"
#include "emissive_config.h"
#include <iostream>
#include <string>

int main(int argc, char** argv) {
    if (argc<3) {
        std::cerr << "Usage:" << argv[0] << " <input_file> <config_type>" << std::endl;
        return 1;
    }

    Mesh mesh;
    if (!mesh.load(argv[1])) {
        std::cerr << "Failed to load mesh from " << argv[1] << std::endl;
        return 1;
    }

    EmissiveConfig config = (std::string(argv[2])=="lego_brick") ? LEGO_BRICK : SMALL_ROOM;
    setEmissiveTriangles(mesh, config); 

    Radiosity radiosity(mesh);
    radiosity.computeRadiosityAntiradiance();

    const auto& combined = radiosity.getRadiance();  // final combined radiance (L - A)

    std::string input_filename = argv[1];
    std::string output_filename = input_filename.substr(0, input_filename.find_last_of("."))+"_output.off";

    if (!mesh.save(output_filename.c_str(), combined)) {
        std::cerr << "Failed to save mesh to " << output_filename << std::endl;
        return 1; 
    }

    std::cout << "Saved output to " << output_filename << std::endl;
    return 0;
}
