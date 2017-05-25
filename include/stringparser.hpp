//
// Created by aaron on 5/14/17.
// stringparse.hpp
//

#ifndef MPM_V2_STRINGPARSE_HPP
#define MPM_V2_STRINGPARSE_HPP

#include <stdlib.h>
#include <string>
#include <vector>

class StringParser{
public:
    static std::string stringRemoveSpaces(std::string); //remove spaces from string
    static std::vector<std::string> stringSplitString(std::string, char); //split string by character
    static std::string stringRemoveComments(std::string); //remove everything after first "#" in string
    static std::string stringRemoveBraces(std::string); //remove braces from string ("{","}")
    static int stringFindStringID(std::vector<std::string>, std::string); //find id of string from list
    static std::string stringRemoveQuotes(std::string); //remove quotes from string(",')
};

#endif //MPM_V2_STRINGPARSE_HPP
