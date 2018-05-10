//
// Created by aaron on 5/9/18.
// parser.hpp
//

#ifndef MPM_V3_PARSER_HPP
#define MPM_V3_PARSER_HPP

#include <stdlib.h>
#include <string>
#include <vector>
#include <regex>
#include <algorithm>

class Parser{
public:
    //remove spaces from string
    static std::string removeSpaces(std::string s){
        s.erase(std::remove_if(s.begin(), s.end(), isspace), s.end());
        return s;
    }

    //split string by character
    static std::vector<std::string> splitString(std::string s,char delim) {
        //credit to Evan Teran
        //found on stack exchange
        std::stringstream ss;
        std::string tmp;
        std::vector<std::string> elems;
        ss.str(s);
        while (std::getline(ss, tmp, delim)) {
            elems.push_back(tmp);
        }
        return elems;
    }

    //remove everything after first "#" in string
    static std::string removeComments(std::string s) {
        std::vector<std::string> svec;
        svec = splitString(s,'#');
        if (svec.size()>0) {
            return svec[0];
        } else {
            return std::string();
        }
    }

    //remove braces from string ("{","}")
    static std::string removeBraces(std::string s){
        std::vector<std::string> svec;
        std::stringstream ss;
        //std::cout << s << std::endl;
        svec = splitString(s,'{');
        for(size_t i=0;i<svec.size();i++){
            ss << svec[i];
        }
        //std::cout << ss.str() << std::endl;
        svec = splitString(ss.str(),'}');
        ss.str("");
        ss.clear();
        //std::cout << ss.str() << std::endl;
        for(size_t i=0;i<svec.size();i++){
            ss << svec[i];
        }
        //std::cout << ss.str() << std::endl;
        return ss.str();
    }

    //find id of string from list
    static int findStringID(std::vector<std::string> svec, std::string s){
        for (size_t i=0;i<svec.size();i++){
            if(svec[i].compare(s) == 0){
                return i;
            }
        }
        return -1;
    }

    //remove quotes from string(",')
    static std::string removeQuotes(std::string s){
        std::vector<std::string> svec;
        std::stringstream ss;
        //std::cout << s << std::endl;
        svec = splitString(s,'\'');
        for(size_t i=0;i<svec.size();i++){
            ss << svec[i];
        }
        //std::cout << ss.str() << std::endl;
        svec = splitString(ss.str(),'\"');
        ss.str("");
        ss.clear();
        //std::cout << ss.str() << std::endl;
        for(size_t i=0;i<svec.size();i++){
            ss << svec[i];
        }
        //std::cout << ss.str() << std::endl;
        return ss.str();
    }

    //add "/" if string ends without it
    static std::string makeDirectory(std::string s) {
        if (s.empty() || s.back() == '/'){
            return s;
        } else {
            return s + "/";
        }
    }
};

#endif //MPM_V3_PARSER_HPP
