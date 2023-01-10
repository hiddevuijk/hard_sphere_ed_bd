#ifndef GUARD_READ_POSITIONS_H
#define GUARD_READ_POSITIONS_H

#include <fstream>
#include <string>
#include <vector>
#include <sstream>

#include "vec3.h"

std::vector<Vec3> read_positions(std::string positions_data_name)
{

  std::vector<Vec3> positions;
  std::string line;

  std::fstream file;
  file.open(positions_data_name);

  if (!file.is_open()) {
    std::cerr << "ERROR: \n";
    std::cerr << positions_data_name;
    std::cerr << " could not be opened.\n";
    exit(1);
  }

  Vec3 temp;
  while (std::getline(file, line)) {
    std::stringstream str_stream(line);
    str_stream >> temp.x;
    str_stream >> temp.y;
    str_stream >> temp.z;
    positions.push_back(temp);
  }

  file.close();

  return positions;
}


#endif
