#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cstdlib>

namespace Utils{

//Some string manipulation functions for saving files. pad_int(1234,5) returns "01234".
inline std::string pad_int(int arg, int padcount) {
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(padcount) << arg;
    return ss.str();
}

//Returns a file name in the form of "prefix00###suffix". For example "image0032.bmp"
inline std::string getFilename(std::string prefix, int num, int padcount, std::string suffix) {
    return prefix + pad_int(num, padcount) + suffix;
}

inline float randf(float min, float max) {
    return min + float(rand())*(max-min)/RAND_MAX;
}


}
#endif
