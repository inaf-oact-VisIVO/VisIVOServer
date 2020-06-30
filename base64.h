//
//  base64.h
//  VisIVOServer_New072012
//
//  Created by Fabio Vitello on 29/09/14.
//
//

#include <string>

std::string base64_encode(unsigned char const* , unsigned int len);
std::string base64_decode(std::string const& s);