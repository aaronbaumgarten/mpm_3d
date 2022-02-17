//
// Created by aaron on 5/17/18.
// config.cpp
//

#include <functional>
#include <map>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <stdexcept>
#include <Eigen/Core>

#include "parser.hpp"
#include "job.hpp"
#include "config.hpp"

#include "mpm_objects.hpp"
#include "registry.hpp"

#include <Eigen/Core>
#include "mpm_vector.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensor.hpp"
#include "mpm_tensorarray.hpp"

// Type of a function called when parsing a section.
using section_callback = std::function<void(Job *, std::ifstream &, Configurator &)>;

// Type of a function called when parsing the value of a key-value pair.
using value_callback = std::function<void(std::string)>;

// We we will be calling these when we configure the job, so forward declare
// them and add to the section table.
static void parseJobSection(Job *job, std::ifstream &fin, Configurator &c);
static void parseSerializerSection(Job *job, std::ifstream &fin, Configurator &c);
static void parseDriverSection(Job *job, std::ifstream &fin, Configurator &c);
static void parseSolverSection(Job *job, std::ifstream &fin, Configurator &c);
static void parseGridSection(Job *job, std::ifstream &fin, Configurator &c);
static void parseBodySection(Job *job, std::ifstream &fin, Configurator &c);
static void parseContactSection(Job *job, std::ifstream &fin, Configurator &c);

// A section in the configuration file is roughly something that has a braced
// list of key-value pairs denoted by a header, e.g.:
//
// my_section_header
// {
//     key = value
// }
//
// This is the list of known sections and associated callbacks to parse the
// contents.
const std::map<const std::string, section_callback> sections{
    { "job",        parseJobSection },
    { "serializer", parseSerializerSection },
    { "driver",     parseDriverSection },
    { "solver",     parseSolverSection },
    { "grid",       parseGridSection },
    { "body",       parseBodySection },
    { "contact",    parseContactSection },
};

/*----------------------------------------------------------------------------*/
//

template <typename T = std::string>
T stringTo(std::string arg) {
    return arg;
}

template <>
int stringTo<int>(std::string arg) {
  return std::stoi(arg);
}

template <>
double stringTo<double>(std::string arg) {
  return std::stod(arg);
}

template <typename T>
static std::vector<T> parsePropertyList(std::string &s) {
    std::vector<T> propertyList{};
    auto const tokens = Parser::splitString(s, ',');
    for (auto && token: tokens) {
        propertyList.push_back(stringTo<T>(token));
    }
    return propertyList;
}

// Helper for iterating over the key value pairs and
// applying the callback from table. Used by the parseNamedSection functions.
static void handleKeyValuePairs(const std::map<std::string, std::string> &inputs,
                                const std::map<std::string, value_callback> &knownKeyValuePairs);

// Helper for creating mapping of key value pairs present in a section. The
// keys are single-valued, and a later line with the same key will overwrite an
// earlier one (a warning is generated for this).
static std::map<std::string, std::string> parseSectionToKeyValuePairs(std::ifstream &fin);

/*----------------------------------------------------------------------------*/
//
void Configurator::init(std::string fileIN){
    serializer_registry = Registry<Serializer>();
    driver_registry = Registry<Driver>();
    solver_registry = Registry<Solver>();

    grid_registry = Registry<Grid>();
    body_registry = Registry<Body>();
    points_registry = Registry<Points>();
    nodes_registry = Registry<Nodes>();
    contact_registry = Registry<Contact>();

    material_registry = Registry<Material>();
    boundary_registry = Registry<Boundary>();

    file = fileIN;
    return;
}


/*----------------------------------------------------------------------------*/
//
void Configurator::setMainPath(std::string program){
    std::vector<std::string> svec;
    svec = Parser::splitString(program,'/');
    std::string filepath = "";
    for (int i=0; i<(svec.size()-1);i++){
        filepath += svec[i];
        filepath += "/";
    }
    mainpath = filepath;
    return;
}

/*----------------------------------------------------------------------------*/
//
void Configurator::checkConfigFile(std::string filename) {
    std::map<std::string, int> check;
    for (auto &&kv: sections) {
        check[kv.first] = 0;
    }

    int id;
    std::string line;
    std::ifstream fin(filename);
    if (fin.is_open()){
        while (std::getline(fin,line)){
            line = Parser::removeComments(line);
            line = Parser::removeSpaces(line);

            // Found a header we know about, mark it in the checked table.
            if (sections.count(line) == 1) {
                check[line]++;
            }
        }

        // If we're missing a section inform the user. This is not a hard
        // error, since some sims don't care about things like contact.
        for (auto &&kv: check) {
            if (kv.second == 0) {
                std::cout << "Could not find \"" << kv.first << "\" in configuration file." << std::endl;
            }
        }
        fin.close();
    } else {
        std::cout << "ERROR: Unable to open file: " << filename << std::endl;
    }
    return;
}

/*----------------------------------------------------------------------------*/
//
int Configurator::configureJob(Job* job){
    //declare variables and strings and vectors
    std::string line;

    //declare fstream from filename
    std::ifstream fin(file);

    //read config file and configure
    if (fin.is_open()) {
        //if open, read lines
        while (std::getline(fin,line)){
            //remove spaces and comments from line
            line = Parser::removeComments(line);
            line = Parser::removeSpaces(line);

            //check if line matches header
            if (line.size()>0){
                // Match to something in the section map.
                if (sections.count(line) == 1) {
                    auto const configure = sections.at(line);
                    configure(job, fin, *this);
                } else {
                    std::cerr << "Unknown header \"" << line
                        << "\", ignoring it.\n";
                }
            }
        }
        //close file
        fin.close();

        std::cout << "Job Configured.\n" << std::endl;
    } else {
        std::cout << "ERROR: Unable to open file: " << file << std::endl;
        return 0;
    }

    return 1;
}

/*----------------------------------------------------------------------------*/
//
static void parseJobSection(Job *job, std::ifstream &fin, Configurator &c) {
    auto const inputs = parseSectionToKeyValuePairs(fin);

    const std::map<std::string, value_callback> knownKeyValuePairs{
        { "dt",      [&](std::string value) { job->dt = std::stod(value); } },
        { "t",       [&](std::string value) { job->t = std::stod(value); job->t0 = job->t; } },
        { "TYPE",    [&](std::string value) { job->assignJobType(std::stoi(value)); } },
        { "threads", [&](std::string value) { job->thread_count = (std::stoi(value)); } },
    };

    handleKeyValuePairs(inputs, knownKeyValuePairs);
    return;
}

/*----------------------------------------------------------------------------*/
//
static void parseSerializerSection(Job *job, std::ifstream &fin, Configurator &c) {
    auto const inputs = parseSectionToKeyValuePairs(fin);
    MPMObject tmp;

    const std::map<std::string, value_callback> knownKeyValuePairs{
        { "class",          [&](std::string value) { job->serializer = c.serializer_registry.get_object(value); } },
        { "properties",     [&](std::string value) { tmp.fp64_props = parsePropertyList<double>(value); } },
        { "int-properties", [&](std::string value) { tmp.int_props = parsePropertyList<int>(value); } },
        { "str-properties", [&](std::string value) { tmp.str_props = parsePropertyList<std::string>(value); } },
    };

    handleKeyValuePairs(inputs, knownKeyValuePairs);

    if (job->serializer == nullptr) {
        throw std::runtime_error("Serializer was nullptr -- was the class set in the config file?");
    } else {
        job->serializer->fp64_props = tmp.fp64_props;
        job->serializer->int_props = tmp.int_props;
        job->serializer->str_props = tmp.str_props;
    }

    return;
}

/*----------------------------------------------------------------------------*/
//
static void parseDriverSection(Job *job, std::ifstream &fin, Configurator &c) {
    auto const inputs = parseSectionToKeyValuePairs(fin);
    MPMObject tmp;

    const std::map<std::string, value_callback> knownKeyValuePairs{
        { "class",          [&](std::string value) { job->driver = c.driver_registry.get_object(value); } },
        { "properties",     [&](std::string value) { tmp.fp64_props = parsePropertyList<double>(value); } },
        { "int-properties", [&](std::string value) { tmp.int_props = parsePropertyList<int>(value); } },
        { "str-properties", [&](std::string value) { tmp.str_props = parsePropertyList<std::string>(value); } },
    };

    handleKeyValuePairs(inputs, knownKeyValuePairs);

    if (job->driver == nullptr) {
        throw std::runtime_error("Driver was nullptr -- was the class set in the config file?");
    } else {
        job->driver->fp64_props = tmp.fp64_props;
        job->driver->int_props = tmp.int_props;
        job->driver->str_props = tmp.str_props;
    }

    std::cout << "Driver Configured: " << job->driver->object_name << std::endl;
    return;
}

/*----------------------------------------------------------------------------*/
//
static void parseSolverSection(Job *job, std::ifstream &fin, Configurator &c) {
    auto const inputs = parseSectionToKeyValuePairs(fin);
    MPMObject tmp;

    const std::map<std::string, value_callback> knownKeyValuePairs{
        { "class",          [&](std::string value) { job->solver = c.solver_registry.get_object(value); } },
        { "properties",     [&](std::string value) { tmp.fp64_props = parsePropertyList<double>(value); } },
        { "int-properties", [&](std::string value) { tmp.int_props = parsePropertyList<int>(value); } },
        { "str-properties", [&](std::string value) { tmp.str_props = parsePropertyList<std::string>(value); } },
    };

    handleKeyValuePairs(inputs, knownKeyValuePairs);

    if (job->solver == nullptr) {
        throw std::runtime_error("Solver was nullptr -- was the class set in the config file?");
    } else {
        job->solver->fp64_props = tmp.fp64_props;
        job->solver->int_props = tmp.int_props;
        job->solver->str_props = tmp.str_props;
    }

    std::cout << "Solver Configured: " << job->solver->object_name << std::endl;
    return;
}

/*----------------------------------------------------------------------------*/
//
static void parseGridSection(Job *job, std::ifstream &fin, Configurator &c) {
    auto const inputs = parseSectionToKeyValuePairs(fin);
    MPMObject tmp;

    const std::map<std::string, value_callback> knownKeyValuePairs{
        { "class",          [&](std::string value) { job->grid = c.grid_registry.get_object(value); } },
        { "properties",     [&](std::string value) { tmp.fp64_props = parsePropertyList<double>(value); } },
        { "int-properties", [&](std::string value) { tmp.int_props = parsePropertyList<int>(value); } },
        { "str-properties", [&](std::string value) { tmp.str_props = parsePropertyList<std::string>(value); } },
    };

    handleKeyValuePairs(inputs, knownKeyValuePairs);

    if (job->grid == nullptr) {
        throw std::runtime_error("Grid was nullptr -- was the class set in the config file?");
    } else {
        job->grid->fp64_props = tmp.fp64_props;
        job->grid->int_props = tmp.int_props;
        job->grid->str_props = tmp.str_props;

        //assign grid dimension from job type
        if (job->JOB_TYPE == job->JOB_1D) {
            job->grid->GRID_DIM = 1;
        } else if (job->JOB_TYPE == job->JOB_2D) {
            job->grid->GRID_DIM = 2;
        } else if (job->JOB_TYPE == job->JOB_3D) {
            job->grid->GRID_DIM = 3;
        } else if (job->JOB_TYPE == job->JOB_2D_OOP) {
            job->grid->GRID_DIM = 2; //this is important, job->DIM =/= job->grid->GRID_DIM
        } else if (job->JOB_TYPE == job->JOB_AXISYM) {
            job->grid->GRID_DIM = 2; //this is important, job->DIM =/= job->grid->GRID_DIM
        } else {
            std::cerr << "Job doesn't have defined type for input " << job->JOB_TYPE << "." << std::endl;
        }
    }

    std::cout << "Grid Configured: " << job->grid->object_name << std::endl;
    return;
}

/*----------------------------------------------------------------------------*/
//
static void parseBodySection(Job *job, std::ifstream &fin, Configurator &c) {
    auto const inputs = parseSectionToKeyValuePairs(fin);
    MPMObject tmp;
    MPMObject point_tmp;
    MPMObject node_tmp;
    MPMObject material_tmp;
    MPMObject boundary_tmp;

    job->bodies.push_back(std::unique_ptr<Body>(nullptr));
    job->activeBodies.push_back(0); //set body to inactive
    const size_t id = job->bodies.size() - 1;

    std::string name;
    std::string filename;

    const std::map<std::string, value_callback> knownKeyValuePairs{
        { "class",              [&](std::string value) {
                job->bodies[id] = c.body_registry.get_object(value);
                job->bodies[id]->id = id;
            }
        },
        { "name",               [&](std::string value) { name = value; } },
        { "point-file",         [&](std::string value) { filename = value; } },
        { "properties",         [&](std::string value) { tmp.fp64_props = parsePropertyList<double>(value); } },
        { "int-properties",     [&](std::string value) { tmp.int_props = parsePropertyList<int>(value); } },
        { "str-properties",     [&](std::string value) { tmp.str_props = parsePropertyList<std::string>(value); } },

        { "point-class",        [&](std::string value) { point_tmp.object_name = value; } },
        { "point-props",        [&](std::string value) { point_tmp.fp64_props = parsePropertyList<double>(value); } },
        { "point-int-props",    [&](std::string value) { point_tmp.int_props = parsePropertyList<int>(value); } },
        { "point-str-props",    [&](std::string value) { point_tmp.str_props = parsePropertyList<std::string>(value); } },

        { "node-class",         [&](std::string value) { node_tmp.object_name = value; } },
        { "node-props",         [&](std::string value) { node_tmp.fp64_props = parsePropertyList<double>(value); } },
        { "node-int-props",     [&](std::string value) { node_tmp.int_props = parsePropertyList<int>(value); } },
        { "node-str-props",     [&](std::string value) { node_tmp.str_props = parsePropertyList<std::string>(value); } },

        { "material-class",     [&](std::string value) { material_tmp.object_name = value; } },
        { "material-props",     [&](std::string value) { material_tmp.fp64_props = parsePropertyList<double>(value); } },
        { "material-int-props", [&](std::string value) { material_tmp.int_props = parsePropertyList<int>(value); } },
        { "material-str-props", [&](std::string value) { material_tmp.str_props = parsePropertyList<std::string>(value); } },

        { "boundary-class",     [&](std::string value) { boundary_tmp.object_name = value; } },
        { "boundary-props",     [&](std::string value) { boundary_tmp.fp64_props = parsePropertyList<double>(value); } },
        { "boundary-int-props", [&](std::string value) { boundary_tmp.int_props = parsePropertyList<int>(value); } },
        { "boundary-str-props", [&](std::string value) { boundary_tmp.str_props = parsePropertyList<std::string>(value); } },
    };

    handleKeyValuePairs(inputs, knownKeyValuePairs);

    if (job->bodies[id] == nullptr) {
        throw std::runtime_error("Body was nullptr -- was the class set in the config file?");
    } else {
        //if body is valid
        //job->activeBodies[id] = 1;
        job->bodies[id]->name = name;
        job->bodies[id]->fp64_props = tmp.fp64_props;
        job->bodies[id]->int_props = tmp.int_props;
        job->bodies[id]->str_props = tmp.str_props;

        if (point_tmp.object_name.empty()) {
            std::cerr << "No Points!" << std::endl;
        } else {
            //if a point object is defined
            job->bodies[id]->points = c.points_registry.get_object(point_tmp.object_name);
            if (job->bodies[id]->points == nullptr) {
                throw std::runtime_error("Points was nullptr -- was the class set in the config file?");
            } else {
                job->bodies[id]->points->file = filename;
                job->bodies[id]->points->fp64_props = point_tmp.fp64_props;
                job->bodies[id]->points->int_props = point_tmp.int_props;
                job->bodies[id]->points->str_props = point_tmp.str_props;

                //setup job and points
                job->bodies[id]->points->readFromFile(job, job->bodies[id].get(), filename);
                job->activeBodies[id] = 1;

                //initialize non-active material/boundary
                job->bodies[id]->activeBoundary = 0;
                job->bodies[id]->activeMaterial = 0;
            }
        }

        if (node_tmp.object_name.empty()) {
            std::cerr << "No Nodes!" << std::endl;
        } else {
            //if a node object is defined
            job->bodies[id]->nodes = c.nodes_registry.get_object(node_tmp.object_name);
            if (job->bodies[id]->nodes == nullptr){
                throw std::runtime_error("Nodes was nullptr -- was the class set in the config file?");
            } else {
                job->bodies[id]->nodes->fp64_props = node_tmp.fp64_props;
                job->bodies[id]->nodes->int_props = node_tmp.int_props;
                job->bodies[id]->nodes->str_props = node_tmp.str_props;
            }
        }

        if (material_tmp.object_name.empty()) {
            std::cerr << "No Material!" << std::endl;
        } else {
            //if a material object is defined
            job->bodies[id]->material = c.material_registry.get_object(material_tmp.object_name);
            if (job->bodies[id]->material == nullptr) {
                throw std::runtime_error("Material was nullptr -- was the class set in the config file?");
            } else {
                job->bodies[id]->material->fp64_props = material_tmp.fp64_props;
                job->bodies[id]->material->int_props = material_tmp.int_props;
                job->bodies[id]->material->str_props = material_tmp.str_props;

                job->bodies[id]->activeMaterial = 1;
            }
        }

        if (boundary_tmp.object_name.empty()){
            std::cerr << "No Boundary!" << std::endl;
        } else {
            //if a boundary object is defined
            job->bodies[id]->boundary = c.boundary_registry.get_object(boundary_tmp.object_name);
            if (job->bodies[id]->boundary == nullptr) {
                throw std::runtime_error("Boundary was nullptr -- was the class set in the config file?");
            } else {
                job->bodies[id]->boundary->fp64_props = boundary_tmp.fp64_props;
                job->bodies[id]->boundary->int_props = boundary_tmp.int_props;
                job->bodies[id]->boundary->str_props = boundary_tmp.str_props;

                job->bodies[id]->activeBoundary = 1;
            }
        }
    }

    std::cout << "Body [" << id << ", " << job->bodies[id]->name
        << "] Configured: " << job->bodies[id]->object_name << std::endl;
    return;
}

/*----------------------------------------------------------------------------*/
//
static void parseContactSection(Job *job, std::ifstream &fin, Configurator &c) {
    auto const inputs = parseSectionToKeyValuePairs(fin);
    MPMObject tmp;

    job->contacts.push_back(std::unique_ptr<Contact>(nullptr));
    job->activeContacts.push_back(0); //set body to inactive
    const size_t id = job->contacts.size() - 1;

    std::string name;

    const std::map<std::string, value_callback> knownKeyValuePairs{
        { "class",          [&](std::string value) {
                job->contacts[id] = c.contact_registry.get_object(value);
                job->contacts[id]->id = id;
            }
        },
        { "name",           [&](std::string value) { name = value; } },
        { "properties",     [&](std::string value) { tmp.fp64_props = parsePropertyList<double>(value); } },
        { "int-properties", [&](std::string value) { tmp.int_props = parsePropertyList<int>(value); } },
        { "str-properties", [&](std::string value) { tmp.str_props = parsePropertyList<std::string>(value); } },
    };

    handleKeyValuePairs(inputs, knownKeyValuePairs);

    if (job->contacts[id] == nullptr) {
        throw std::runtime_error("Contacts[" + std::to_string(id) + "] was nullptr -- was the class set in the config file?");
    } else {
        job->activeContacts[id] = 1;
        job->contacts[id]->name = name;
        job->contacts[id]->fp64_props = tmp.fp64_props;
        job->contacts[id]->int_props = tmp.int_props;
        job->contacts[id]->str_props = tmp.str_props;
    }

    std::cout << "Grid Configured: " << job->grid->object_name << std::endl;
    return;
}

/*----------------------------------------------------------------------------*/
//
static std::map<std::string, std::string> parseSectionToKeyValuePairs(std::ifstream &fin) {
    std::map<std::string, std::string> keyValuePairs;
    std::string line;

    std::getline(fin, line);
    line = Parser::removeComments(line);
    line = Parser::removeSpaces(line);
    if (line.compare("{") == 0) {
        while (std::getline(fin, line)){
            line = Parser::removeComments(line);
            line = Parser::removeSpaces(line);

            if (line.compare("}") == 0) {
                break;
            }

            line = Parser::removeBraces(line);
            line = Parser::removeQuotes(line);
            auto const tokens = Parser::splitString(line, '=');
            if (tokens.size() > 1) {
                if (tokens.size() != 2) {
                    std::cerr << "Unexpected split of line \"" << line
                        << "\" into more than 2 tokens, but continuing anyways.\n";
                }

                auto const key = tokens[0];
                auto const value = tokens[1];
                if (keyValuePairs.count(key) >= 1) {
                    std::cerr << "Already have seen key \"" << key
                        << "\", but overriding old seen value \""
                        << keyValuePairs[key] << "\" with new value \""
                        << value << "\".";
                }

                keyValuePairs[key] = value;
            }
        }
    }

    return keyValuePairs;
}

/*----------------------------------------------------------------------------*/
//
static void handleKeyValuePairs(const std::map<std::string, std::string> &inputs,
                                const std::map<std::string, value_callback> &knownKeyValuePairs) {

    // Iterate through all the given key-value input pairs and check if we know
    // how to handle them via the map knownKeyValuePairs. If so, run the
    // callback function which evalues the value string -- this will usually be
    // a lambda set by the parsing function.

    // Unknown keys will generate a warning, but not an error. Defining a key
    // multiple times will simply use the last set value and ignore earlier
    // ones (warnings are generated for that case when creating the inputs map,
    // not here).
    for (auto && kv : inputs) {
        if (knownKeyValuePairs.count(kv.first) == 1) {
            auto parseValue = knownKeyValuePairs.at(kv.first);
            parseValue(kv.second);
        } else {
            std::cerr << "Unknown key \"" << kv.first << "\" (with value \""
                << kv.second << "\") provided. Ignoring.\n";
        }
    }

    return;
}

/*----------------------------------------------------------------------------*/
