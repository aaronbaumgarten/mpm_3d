//
// Created by aaron on 5/14/17.
// stringparse.cpp
//

#include <stdlib.h>
#include <string>
#include <vector>
#include <regex>
#include <algorithm>

#include "stringparser.hpp"

static std::string StringParser::stringRemoveSpaces(std::string s){
    s.erase(std::remove_if(s.begin(), s.end(), isspace), s.end());
    return s;
}

static std::vector<std::string> StringParser::stringSplitString(std::string s,char delim){
    //credit to Evan Teran
    //found on stack exchange
    std::stringstream ss;
    std::string tmp;
    std::vector<std::string> elems;
    ss.str(s);
    while (std::getline(ss, tmp, delim)){
        elems.push_back(tmp);
    }
    return elems;
}

static std::string StringParser::stringRemoveComments(std::string s) {
    std::vector<std::string> svec;
    svec = StringParser::stringSplitString(s,'#');
    if (svec.size()>0) {
        return svec[0];
    } else {
        return std::string();
    }
}

static std::string StringParser::stringRemoveBraces(std::string s){
    std::vector<std::string> svec;
    std::stringstream ss;
    //std::cout << s << std::endl;
    svec = StringParser::stringSplitString(s,'{');
    for(size_t i=0;i<svec.size();i++){
        ss << svec[i];
    }
    //std::cout << ss.str() << std::endl;
    svec = StringParser::stringSplitString(ss.str(),'}');
    ss.str("");
    ss.clear();
    //std::cout << ss.str() << std::endl;
    for(size_t i=0;i<svec.size();i++){
        ss << svec[i];
    }
    //std::cout << ss.str() << std::endl;
    return ss.str();
}

static int StringParser::stringffFindStringID(std::vector<std::string> svec, std::string s){
    for (size_t i=0;i<svec.size();i++){
        if(svec[i].compare(s) == 0){
            return i;
        }
    }
    return -1;
}

static std::string StringParser::stringRemoveQuotes(std::string s){
    std::vector<std::string> svec;
    std::stringstream ss;
    //std::cout << s << std::endl;
    svec = StringParser::stringSplitString(s,'\'');
    for(size_t i=0;i<svec.size();i++){
        ss << svec[i];
    }
    //std::cout << ss.str() << std::endl;
    svec = StringParser::stringSplitString(ss.str(),'\"');
    ss.str("");
    ss.clear();
    //std::cout << ss.str() << std::endl;
    for(size_t i=0;i<svec.size();i++){
        ss << svec[i];
    }
    //std::cout << ss.str() << std::endl;
    return ss.str();
}