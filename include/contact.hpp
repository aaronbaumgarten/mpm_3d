//
// Created by aaron on 12/23/16.
// contact.hpp
//

#include <stdlib.h>
#include <vector>
#include <string>

#ifndef MPM_3D_CONTACT_HPP
#define MPM_3D_CONTACT_HPP

class job_t;

class Contact{
public:
    size_t id;
    std::string contact_filename;
    int use_builtin;
    std::vector<int> bodyIDs;

    void *handle;

    void (*contact_init)(job_t*, size_t);
    void (*resolve_contact)(job_t*, size_t);

    std::vector<double> fp64_props;
    std::vector<int> int_props;
    size_t num_fp64_props;
    size_t num_int_props;

    Contact();
    Contact(std::string,size_t,std::vector<int>,std::vector<double>,std::vector<int>);
    ~Contact();
    void setContact(std::string,size_t,std::vector<int>,std::vector<double>,std::vector<int>);
};

#endif //MPM_3D_CONTACT_HPP
