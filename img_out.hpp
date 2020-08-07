
#ifndef _PGM_OUT_H
#define _PGM_OUT_H


#include <iostream>
#include <vector>
#include <fstream>
#include <limits>
#include <string>
#include <algorithm>

void write_image( std::string fn, unsigned char *data, int width, int height, bool invert = false );


#endif
