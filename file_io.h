#ifndef FILE_IO_H
#define FILE_IO_H

#include "mesh.h"
#include <string>

bool loadObjFile(const std::string & filename, Mesh& mesh);
void saveOffFile(const std::string& filename, const Mesh& mesh, const std::vector<Vertex>&colors);
void validateOffFile(const std::string& filename);

#endif
