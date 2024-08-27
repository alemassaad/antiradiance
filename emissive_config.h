#ifndef EMISSIVE_CONFIG_H
#define EMISSIVE_CONFIG_H
#include "mesh.h"

enum EmissiveConfig {
    LEGO_BRICK,
    SMALL_ROOM,
    BEAR
};

void setEmissiveTriangles(Mesh& mesh, EmissiveConfig config);

#endif
