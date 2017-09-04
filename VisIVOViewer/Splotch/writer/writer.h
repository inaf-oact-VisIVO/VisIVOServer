#ifndef WRITE_TGA_H
#define WRITE_TGA_H

#include <string>
#include "cxxsupport/arr.h"
#include "cxxsupport/paramfile.h"
#include "kernel/colour.h"

void write_tga(paramfile &params, const arr2<COLOUR> &pic,
  const std::string &frame_name);

void write_tga_rle(paramfile &params, const arr2<COLOUR> &pic,
  const std::string &frame_name);

void write_ppm_ascii(paramfile &params, const arr2<COLOUR> &pic,
  const std::string &frame_name);

void write_ppm_bin(paramfile &params, const arr2<COLOUR> &pic,
  const std::string &frame_name);

#endif
