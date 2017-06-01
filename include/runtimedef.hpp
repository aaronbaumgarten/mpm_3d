//
// Created by aaron on 5/28/17.
// runtimedef.hpp
//

#ifndef MPM_V2_RUNTIMEDEF_HPP
#define MPM_V2_RUNTIMEDEF_HPP

#include <stdlib.h>
#include <string>
#include <vector>

class RunTimeDef{
public:
    std::string fullpath; //long path to file
    std::string filename; //name of file
    std::string filepath; //directory of file for access
    std::vector<double> fp64_props; //double properties
    std::vector<int> int_props; //integer properties
    std::vector<std::string> str_props; //string properties
    void *handle; //.so file handle
};

#endif //MPM_V2_RUNTIMEDEF_HPP
