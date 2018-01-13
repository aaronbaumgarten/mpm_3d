//
// Created by aaron on 10/18/17.
// minimal_vtk.cpp
//

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <sys/stat.h>
#include <math.h>

#include "runtimedef.hpp"
#include "stringparser.hpp"

#include "job.hpp"
#include "config.hpp"

#include "serializer.hpp"
#include "driver.hpp"
#include "solver.hpp"

#include "grid.hpp"
#include "contact.hpp"

#include "body.hpp"
#include "nodes.hpp"
#include "points.hpp"

#include "material.hpp"
#include "boundary.hpp"

//file storing and writing variables
std::string frameDirectory;
std::string outputDirectory;
std::string outputName;

//job writing and sampling
size_t sampledFrames;
double sampleRate;
double t_last_frame;

//frame writing variables
std::string pfilename;
std::string nfilename;
std::ofstream pfile;
std::ofstream nfile;
size_t plen;
size_t nlen;
Body* currentBody;

extern "C" int serializerWriteFrame(Job* job); //initialize writing of frame from job
extern "C" void serializerWriteScalarArray(Eigen::VectorXd& scalarArray, std::string scalarName); //call to functions for dumping state into frame file
extern "C" void serializerWriteVectorArray(Eigen::MatrixXd& vectorArray, std::string vectorName); //pass name of vector
extern "C" void serializerWriteTensorArray(Eigen::MatrixXd& tensorArray, std::string tensorName);

extern "C" void serializerInit(Job* job); //initialize serializer
extern "C" void serializerSetMainPath(Job* job, std::string program); //set path to main program
extern "C" std::string serializerSaveState(Job* job); //save job state (output string for output directory)
extern "C" int serializerLoadState(Job* job, std::string fullpath); //load state from given full path

/*----------------------------------------------------------------------------*/

void serializerInit(Job* job){
    if (job->serializer.str_props.size() < 3 || job->serializer.fp64_props.size() < 1){
        std::cout << job->serializer.str_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 4 properties defined ({sampleRate}, {frameDirectory, outputDirectory, outputName}).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        frameDirectory = StringParser::stringMakeDirectory(job->serializer.str_props[0]);
        outputDirectory = StringParser::stringMakeDirectory(job->serializer.str_props[1]);
        outputName = job->serializer.str_props[2];

        sampleRate = job->serializer.fp64_props[0];
        sampledFrames = 0;
        t_last_frame = 0;

        printf("Serializer properties (frameDirectory = %s, outputDirectory = %s, outputName = %s, sampleRate = %g).\n",
               frameDirectory.c_str(), outputDirectory.c_str(), outputName.c_str(), sampleRate);
    }

    std::cout << "Serializer Initialized." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/

int serializerWriteFrame(Job* job){
    if ((job->t - t_last_frame) >= (1.0/sampleRate) || sampledFrames == 0){ //job->t >= sampledFrames/sampleRate){
        Eigen::VectorXd gamma_dot, speed, pressure, displacement;
        t_last_frame = sampledFrames/sampleRate;
        sampledFrames += 1;
        //write frame for each body
        for (size_t b=0;b<job->bodies.size();b++) {
            currentBody = &(job->bodies[b]);
            gamma_dot.resize(currentBody->points.x.rows());
            speed.resize(currentBody->points.x.rows());
            pressure.resize(currentBody->points.x.rows());
            displacement.resize(currentBody->points.x.rows());

            double p, v, u, gdp;
            Eigen::MatrixXd T = job->jobTensor<double>();
            Eigen::MatrixXd L = job->jobTensor<double>();
            Eigen::MatrixXd D = job->jobTensor<double>();
            Eigen::VectorXd tmpVec;
            for (size_t i=0;i<currentBody->points.x.rows();i++){
                tmpVec = job->bodies[b].points.T.row(i);
                T = job->jobTensor<double>(tmpVec.data());
                tmpVec = job->bodies[b].points.L.row(i);
                L = job->jobTensor<double>(tmpVec.data());
                D = 0.5 * (L + L.transpose());
                u = job->bodies[b].points.u.row(i).norm();
                v = job->bodies[b].points.x_t.row(i).norm();
                p = -T.trace() / T.rows();
                gdp = (D - D.trace() / D.rows() * job->jobTensor<double>(Job::IDENTITY)).norm() * std::sqrt(2.0);

                gamma_dot(i) = gdp;
                speed(i) = v;
                pressure(i) = p;
                displacement(i) = u;
            }

            //open point file
            std::stringstream ss;
            ss << "fpd." << job->bodies[b].id << "." << job->bodies[b].name << "." << std::setw(10) << std::setfill('0') << (sampledFrames-1) << ".vtk";
            pfilename = ss.str();
            pfile = std::ofstream(frameDirectory+pfilename,std::ios::trunc);

            ss.str(std::string());
            ss.clear();

            //set length of frame data
            plen = job->bodies[b].points.x.rows();

            if (pfile.is_open()){
                pfile << "# vtk DataFile Version 3.0\n";
                pfile << "Frame: " << (sampledFrames-1) << ", Time: " << job->t << "\n";
                pfile << "ASCII\n";
                pfile << "DATASET UNSTRUCTURED_GRID\n";

                pfile << "POINTS " << plen << " double\n";
                for (size_t i=0;i<plen;i++){
                    //vtk files require x,y,z
                    for (size_t pos = 0; pos < 3; pos++){
                        if (pos < job->bodies[b].points.x.cols() && (job->bodies[b].points.active(i) != 0) && std::isfinite(job->bodies[b].points.x(i,pos))){
                            pfile << job->bodies[b].points.x(i,pos) << " ";
                        } else {
                            pfile << "0 ";
                        }
                    }
                    pfile << "\n";
                }

                pfile << "CELLS " << plen << " " << 2*plen << "\n";
                for (size_t i=0;i<plen;i++){
                    pfile << "1 " << i << "\n";
                }

                pfile << "CELL_TYPES " << plen << "\n";
                for (size_t i=0;i<plen;i++){
                    pfile << "1\n";
                }

                pfile << "POINT_DATA " << plen << "\n";
                //scalars, vectors and tensors here
            } else {
                std::cerr << "Could not open point frame: " << frameDirectory+pfilename << " !" << std::endl;
            }

            serializerWriteScalarArray(displacement,"displacement");
            serializerWriteScalarArray(speed,"speed");
            serializerWriteScalarArray(pressure,"pressure");
            serializerWriteScalarArray(gamma_dot,"gamma_dot");

            pfile.close();
        }

        return 1;
    } else {
        //do not write frame
        return 0;
    }
}

/*----------------------------------------------------------------------------*/

void serializerWriteScalarArray(Eigen::VectorXd& scalarArray, std::string name){
    //check length of array vs. length of open files and write to correct file
    if (pfile.is_open() && scalarArray.rows() == plen){
        //write to point file
        pfile << "SCALARS " << name << " double 1\n";
        pfile << "LOOKUP_TABLE default\n";
        for (size_t i = 0; i < plen; i++){
            if (currentBody->points.active(i) == 1 && std::isfinite(scalarArray(i))) {
                pfile << scalarArray(i) << "\n";
            } else {
                pfile << "0" << "\n";
            }
        }
    }

    if (nfile.is_open() && scalarArray.rows() == nlen){
        //write to point file
        nfile << "SCALARS " << name << " double 1\n";
        nfile << "LOOKUP_TABLE default\n";
        for (size_t i = 0; i < nlen; i++){
            if (currentBody->nodes.active(i) == 1 && std::isfinite(scalarArray(i))) {
                nfile << scalarArray(i) << "\n";
            } else {
                nfile << "0" << "\n";
            }
        }
    }
    return;
}

/*----------------------------------------------------------------------------*/

void serializerWriteVectorArray(Eigen::MatrixXd& vectorArray, std::string name){
    //check length of array vs. length of open files and write to correct file
    if (pfile.is_open() && vectorArray.rows() == plen){
        //write to point file
        pfile << "VECTORS " << name << " double\n";
        for (size_t i = 0; i < plen; i++){
            //vtk format requires x,y,z
            for (size_t pos = 0; pos < 3; pos++){
                if (pos < vectorArray.cols() && (currentBody->points.active(i) == 1) && std::isfinite(vectorArray(i,pos))){
                    pfile << vectorArray(i,pos) << " ";
                } else {
                    pfile << "0 ";
                }
            }
            pfile << "\n";
        }
    }

    if (nfile.is_open() && vectorArray.rows() == nlen){
        //write to point file
        nfile << "VECTORS " << name << " double\n";
        for (size_t i = 0; i < nlen; i++){
            //vtk format requires x,y,z
            for (size_t pos = 0; pos < 3; pos++){
                if (pos < vectorArray.cols() && (currentBody->nodes.active(i) == 1) && std::isfinite(vectorArray(i,pos))){
                    nfile << vectorArray(i,pos) << " ";
                } else {
                    nfile << "0 ";
                }
            }
            nfile << "\n";
        }
    }
    return;
}

/*----------------------------------------------------------------------------*/

void serializerWriteTensorArray(Eigen::MatrixXd& tensorArray, std::string name){
    //check length of array vs. length of open files and write to correct file
    if (pfile.is_open() && tensorArray.rows() == plen){
        //write to point file
        pfile << "TENSORS " << name << " double\n";
        for (size_t i = 0; i < plen; i++){
            //vtk format requires x,y,z
            //brute force this one
            bool finite = true;
            for (size_t pos=0;pos<tensorArray.cols();pos++){
                if (!std::isfinite(tensorArray(i,pos))){
                    finite = false;
                    break;
                }
            }
            if (currentBody->points.active(i) == 0 || !finite){
                pfile << "0 0 0\n";
                pfile << "0 0 0\n";
                pfile << "0 0 0\n";
                pfile << "\n";
            } else if (tensorArray.cols() == 1) {
                pfile << tensorArray(i,0) << " 0 0\n";
                pfile << "0 0 0\n";
                pfile << "0 0 0\n";
                pfile << "\n";
            } else if (tensorArray.cols() == 4) {
                pfile << tensorArray(i,0) << " " << tensorArray(i,1) << " 0\n";
                pfile << tensorArray(i,2) << " " << tensorArray(i,3) << " 0\n";
                pfile << "0 0 0\n";
                pfile << "\n";
            } else if (tensorArray.cols() == 9) {
                pfile << tensorArray(i,0) << " " << tensorArray(i,1) <<  " " << tensorArray(i,2) << "\n";
                pfile << tensorArray(i,3) << " " << tensorArray(i,4) <<  " " << tensorArray(i,5) << "\n";
                pfile << tensorArray(i,6) << " " << tensorArray(i,7) <<  " " << tensorArray(i,8) << "\n";
                pfile << "\n";
            }
        }
    }

    if (nfile.is_open() && tensorArray.rows() == nlen){
        //write to point file
        nfile << "TENSORS " << name << " double\n";
        for (size_t i = 0; i < nlen; i++){
            //vtk format requires x,y,z
            for (size_t pos = 0; pos < 3; pos++){
                //brute force this one
                if (currentBody->nodes.active(i) == 0){
                    nfile << "0 0 0\n";
                    nfile << "0 0 0\n";
                    nfile << "0 0 0\n";
                    nfile << "\n";
                } else if (tensorArray.cols() == 1) {
                    nfile << tensorArray(i,0) << " 0 0\n";
                    nfile << "0 0 0\n";
                    nfile << "0 0 0\n";
                    nfile << "\n";
                } else if (tensorArray.cols() == 4) {
                    nfile << tensorArray(i,0) << " " << tensorArray(i,1) << " 0\n";
                    nfile << tensorArray(i,2) << " " << tensorArray(i,3) << " 0\n";
                    nfile << "0 0 0\n";
                    nfile << "\n";
                } else if (tensorArray.cols() == 9) {
                    nfile << tensorArray(i,0) << " " << tensorArray(i,1) <<  " " << tensorArray(i,2) << "\n";
                    nfile << tensorArray(i,3) << " " << tensorArray(i,4) <<  " " << tensorArray(i,5) << "\n";
                    nfile << tensorArray(i,6) << " " << tensorArray(i,7) <<  " " << tensorArray(i,8) << "\n";
                    nfile << "\n";
                }
            }
        }
    }
    return;
}

/*----------------------------------------------------------------------------*/

void saveStandardProps(RunTimeDef* rtd, std::ofstream &ffile){

    ffile << rtd->filename << "\n";

    ffile << rtd->filepath << "\n";

    for (size_t i=0;i<rtd->fp64_props.size();i++){
        if (i>0){
            ffile << ", ";
        }
        ffile << rtd->fp64_props[i];
    }
    ffile << "\n";

    for (size_t i=0;i<rtd->int_props.size();i++){
        if (i>0){
            ffile << ", ";
        }
        ffile << rtd->int_props[i];
    }
    ffile << "\n";

    for (size_t i=0;i<rtd->str_props.size();i++){
        if (i>0){
            ffile << ",";
        }
        ffile << rtd->str_props[i];
    }
    ffile << "\n";

    return;
}

/*----------------------------------------------------------------------------*/

void loadStandardProps(RunTimeDef* rtd, std::ifstream &fin){
    std::string line;
    std::vector<std::string> lvec;

    std::getline(fin,line); //filename
    rtd->filename = line;

    std::getline(fin,line); //filepath
    rtd->filepath = line;

    std::getline(fin,line); //fp64_props
    lvec = StringParser::stringSplitString(line, ',');
    for (size_t i=0;i<lvec.size();i++){
        rtd->fp64_props.push_back(std::stod(lvec[i]));
    }

    std::getline(fin,line); //int_props
    lvec = StringParser::stringSplitString(line, ',');
    for (size_t i=0;i<lvec.size();i++){
        rtd->int_props.push_back(std::stoi(lvec[i]));
    }

    std::getline(fin,line); //str_props
    lvec = StringParser::stringSplitString(line, ',');
    for (size_t i=0;i<lvec.size();i++){
        rtd->str_props = lvec;
    }
    return;
}

/*----------------------------------------------------------------------------*/

std::string serializerSaveState(Job* job){
    // current date/time based on current system
    time_t now = time(0);

    // convert now to tm struct for UTC
    tm *gmtm = gmtime(&now);

    //create filename/directory name
    std::stringstream s;
    s << "mpm_v2.job."  << outputName << "." << gmtm->tm_mday << "." << gmtm->tm_mon << "." << gmtm->tm_year << ".";
    s << gmtm->tm_hour << "." << gmtm->tm_min << "." << gmtm->tm_sec;

    //create directory
    std::string filefolder = outputDirectory + s.str();
    std::string savefolder = StringParser::stringMakeDirectory(s.str());
    if ((mkdir(filefolder.c_str(), ACCESSPERMS)) == -1){
        std::cerr << "Unable to create save folder: " << filefolder << " !" << std::endl;
        return "ERR";
    }
    //fix name
    filefolder = StringParser::stringMakeDirectory(filefolder);

    //create filename
    std::string filename = "mpm_v2.save_file.default_vtk.txt";

    /**********/
    //create default file
    std::ofstream defaultfile((filefolder + "default_save_file.txt"), std::ios::trunc);
    if (defaultfile.is_open()){
        //save data to file in standard way for main to open later
        defaultfile << job->serializer.filename << "\n";
        defaultfile << job->serializer.filepath << "\n";
        for (size_t i=0;i<job->serializer.fp64_props.size();i++){
            if (i>0){
                defaultfile << ", ";
            }
            defaultfile << job->serializer.fp64_props[i];
        }
        defaultfile << "\n";
        for (size_t i=0;i<job->serializer.int_props.size();i++){
            if (i>0){
                defaultfile << ", ";
            }
            defaultfile << job->serializer.int_props[i];
        }
        defaultfile << "\n";
        for (size_t i=0;i<job->serializer.str_props.size();i++){
            if (i>0){
                defaultfile << ",";
            }
            defaultfile << job->serializer.str_props[i];
        }
        defaultfile << "\n";

        defaultfile << filename << "\n";
    } else {
        std::cout << "Unable to open \"" << filefolder + "defaultsave.txt" << "\" !\n";
        return "ERR";
    }
    /**********/

    //open file
    std::ofstream ffile((filefolder+filename), std::ios::trunc);

    if (ffile.is_open()){
        ffile << "# mpm_v2 save file\n";
        ffile << savefolder << "\n"; //save name of folder
        ffile << "# job\n";
        ffile << job->DIM << "\n";
        ffile << job->dt << "\n";
        ffile << job->t << "\n";
        ffile << job->activeContacts.size() << "\n";
        for (size_t c=0;c<job->activeContacts.size();c++){
            if (c>0){
                ffile << " ";
            }
            ffile << job->activeContacts[c];
        }
        ffile << "\n";
        ffile << job->activeBodies.size() << "\n";
        for (size_t b=0;b<job->activeBodies.size();b++){
            if (b>0){
                ffile << " ";
            }
            ffile << job->activeBodies[b];
        }
        ffile << "\n";

        ffile << "# serializer\n";
        //serializer public properties already saved
        ffile << frameDirectory << "\n";
        ffile << outputDirectory << "\n";
        ffile << outputName << "\n";
        ffile << sampledFrames << "\n";
        ffile << sampleRate << "\n";
        ffile << t_last_frame << "\n";

        ffile << "# driver\n";
        saveStandardProps(&(job->driver),ffile);
        ffile << job->driver.driverSaveState(job,&(job->serializer),filefolder) << "\n"; //store filename for loading state

        ffile << "# solver\n";
        saveStandardProps(&(job->solver),ffile);
        ffile << job->solver.solverSaveState(job,&(job->serializer),filefolder) << "\n";

        ffile << "# grid\n";
        saveStandardProps(&(job->grid),ffile);
        ffile << job->grid.gridSaveState(job,&(job->serializer),filefolder) << "\n";

        for (size_t b=0;b<job->bodies.size();b++){
            ffile << "# body\n";
            ffile << job->bodies[b].bodySaveState(job,&(job->serializer),filefolder) << "\n";
            ffile << job->bodies[b].points.pointsSaveState(job,&(job->bodies[b]),&(job->serializer),filefolder) << "\n";
            ffile << job->bodies[b].nodes.nodesSaveState(job,&(job->bodies[b]),&(job->serializer),filefolder) << "\n";
            saveStandardProps(&(job->bodies[b].material),ffile);
            if (job->bodies[b].activeMaterial != 0) {
                ffile << job->bodies[b].material.materialSaveState(job, &(job->bodies[b]), &(job->serializer), filefolder) << "\n";
            } else {
                ffile << "\n";
            }
            saveStandardProps(&(job->bodies[b].boundary),ffile);
            if (job->bodies[b].activeBoundary != 0) {
                ffile << job->bodies[b].boundary.boundarySaveState(job, &(job->bodies[b]), &(job->serializer), filefolder) << "\n";
            } else {
                ffile << "\n";
            }
        }

        for (size_t c=0;c<job->contacts.size();c++){
            ffile << "# contact\n";
            saveStandardProps(&(job->contacts[c]),ffile);
            if (job->activeContacts[c] != 0) {
                ffile << job->contacts[c].contactSaveState(job, &(job->serializer), filefolder) << "\n";
            } else {
                ffile << "\n";
            }
        }

        ffile.close();
    } else {
        std::cout << "Unable to open \"" << filefolder+filename << "\" !\n";
        return "ERR";
    }

    std::cout << "Job Saved." << std::endl;

    return (filefolder+filename);
}

/*----------------------------------------------------------------------------*/

int serializerLoadState(Job* job, std::string fullpath){

    std::vector<std::string> lvec;
    lvec = StringParser::stringSplitString(fullpath,'/');
    std::string filepath = "";
    for (size_t i=0; i<(lvec.size()-1);i++){
        filepath += lvec[i];
        filepath += "/";
    } //location of file

    std::string line;
    std::stringstream ss;
    std::vector<std::string> headers = {"# job","# serializer","# driver","# solver","# grid","# body","# contact"};

    std::ifstream fin(fullpath);

    if(fin.is_open()){
        std::getline(fin,line); //header
        std::getline(fin,line); //save folder

        /********/

        std::getline(fin,line); //# job
        if (line.compare(headers[0]) != 0){
            std::cout << "Expected \"" << headers[0] << "\" in file " << fullpath << ". Got \"" << line << "\"." << std::endl;
            std::cout << "Loading failed." << std::endl;
            fin.close();
            return 0;
        }

        std::getline(fin,line); //DIM
        job->DIM = std::stod(line);
        std::getline(fin,line); //dt
        job->dt = std::stod(line);
        std::getline(fin,line); //t
        job->t = std::stod(line);
        std::getline(fin,line); //number of contacts
        job->activeContacts.resize(std::stoi(line));
        std::getline(fin,line);
        ss = std::stringstream(line); //read vector into stream
        for (size_t c=0;c<job->activeContacts.size();c++){
            ss >> job->activeContacts[c];
            job->contacts.push_back(Contact());
            job->contacts[c].id = c;
        }
        std::getline(fin,line); //number of bodies
        job->activeBodies.resize(std::stoi(line));
        std::getline(fin,line);
        ss = std::stringstream(line); //read vector into stream
        for (size_t b=0;b<job->activeBodies.size();b++){
            ss >> job->activeBodies[b];
            job->bodies.push_back(Body());
            job->bodies[b].id = b;
        }

        /********/

        std::getline(fin,line); //# serializer
        if (line.compare(headers[1]) != 0){
            std::cout << "Expected \"" << headers[1] << "\" in file " << fullpath << ". Got \"" << line << "\"." << std::endl;
            std::cout << "Loading failed." << std::endl;
            fin.close();
            return 0;
        }
        //do not need to read in serializer public properties
        std::getline(fin,line); //frameDirectory
        frameDirectory = StringParser::stringMakeDirectory(line);
        std::getline(fin,line); //outputDirectory
        outputDirectory = StringParser::stringMakeDirectory(line);
        std::getline(fin,line); //outputName
        outputName = line;
        std::getline(fin,line); //sampledFrames
        sampledFrames = std::stoi(line);
        std::getline(fin,line); //sampleRate
        sampleRate = std::stod(line);
        std::getline(fin,line); //last frame
        t_last_frame = std::stod(line);

        /********/

        std::getline(fin,line); //# driver
        if (line.compare(headers[2]) != 0){
            std::cout << "Expected \"" << headers[2] << "\" in file " << fullpath << ". Got \"" << line << "\"." << std::endl;
            std::cout << "Loading failed." << std::endl;
            fin.close();
            return 0;
        }
        loadStandardProps(&(job->driver),fin);
        if (!(job->driver.filename.empty())) {
            job->driver.driverSetPlugin(job,
                                        job->serializer.mainpath + job->driver.filepath,
                                        job->driver.filename,
                                        job->driver.fp64_props,
                                        job->driver.int_props,
                                        job->driver.str_props);
            std::getline(fin, line); //filename
            if (0 == job->driver.driverLoadState(job, &(job->serializer), filepath + line)) {
                std::cout << "Loading failed." << std::endl;
                fin.close();
                return 0;
            }
        } else {
            std::getline(fin, line); //filename
        }

        /********/

        std::getline(fin,line); //# solver
        if (line.compare(headers[3]) != 0){
            std::cout << "Expected \"" << headers[3] << "\" in file " << fullpath << ". Got \"" << line << "\"." << std::endl;
            std::cout << "Loading failed." << std::endl;
            fin.close();
            return 0;
        }
        loadStandardProps(&(job->solver),fin);
        if (!(job->solver.filename.empty())) {
            job->solver.solverSetPlugin(job,
                                        job->serializer.mainpath + job->solver.filepath,
                                        job->solver.filename,
                                        job->solver.fp64_props,
                                        job->solver.int_props,
                                        job->solver.str_props);
            std::getline(fin, line); //filename
            if (0 == job->solver.solverLoadState(job, &(job->serializer), filepath + line)) {
                std::cout << "Loading failed." << std::endl;
                fin.close();
                return 0;
            }
        } else {
            std::getline(fin, line); //filename
        }

        /********/

        std::getline(fin,line); //# grid
        if (line.compare(headers[4]) != 0){
            std::cout << "Expected \"" << headers[4] << "\" in file " << fullpath << ". Got \"" << line << "\"." << std::endl;
            std::cout << "Loading failed." << std::endl;
            fin.close();
            return 0;
        }
        loadStandardProps(&(job->grid),fin);
        if (!(job->grid.filename.empty())) {
            job->grid.gridSetPlugin(job,
                                    job->serializer.mainpath + job->grid.filepath,
                                    job->grid.filename,
                                    job->grid.fp64_props,
                                    job->grid.int_props,
                                    job->grid.str_props);
            std::getline(fin, line); //filename
            if (0 == job->grid.gridLoadState(job, &(job->serializer), filepath + line)) {
                std::cout << "Loading failed." << std::endl;
                fin.close();
                return 0;
            }
        } else {
            std::getline(fin, line); //filename
        }

        /********/

        std::getline(fin,line); //# body
        size_t body_id = 0;
        while (line.compare(headers[5]) == 0) {
            size_t id = body_id++;

            std::getline(fin,line);
            if (0 == job->bodies[id].bodyLoadState(job, &(job->serializer), filepath + line)) {
                std::cout << "Loading failed." << std::endl;
                fin.close();
                return 0;
            }

            std::getline(fin,line);
            if (0 == job->bodies[id].points.pointsLoadState(job, &(job->bodies[id]), &(job->serializer), filepath + line)) {
                std::cout << "Loading failed." << std::endl;
                fin.close();
                return 0;
            }

            std::getline(fin,line);
            if (0 == job->bodies[id].nodes.nodesLoadState(job, &(job->bodies[id]), &(job->serializer), filepath + line)) {
                std::cout << "Loading failed." << std::endl;
                fin.close();
                return 0;
            }

            loadStandardProps(&(job->bodies[id].material), fin);
            if (!(job->bodies[id].material.filename.empty()) && job->activeBodies[id] != 0 && job->bodies[id].activeMaterial != 0) {
                job->bodies[id].material.materialSetPlugin(job,
                                                           &(job->bodies[id]),
                                                           job->serializer.mainpath + job->bodies[id].material.filepath,
                                                           job->bodies[id].material.filename,
                                                           job->bodies[id].material.fp64_props,
                                                           job->bodies[id].material.int_props,
                                                           job->bodies[id].material.str_props);
                std::getline(fin, line); //filename
                if (0 == job->bodies[id].material.materialLoadState(job, &(job->bodies[id]), &(job->serializer),
                                                                    filepath + line)) {
                    std::cout << "Loading failed." << std::endl;
                    fin.close();
                    return 0;
                }
            } else {
                std::getline(fin, line); //filename
            }

            loadStandardProps(&(job->bodies[id].boundary), fin);
            if (!(job->bodies[id].boundary.filename.empty()) && job->activeBodies[id] != 0 && job->bodies[id].activeBoundary != 0) {
                job->bodies[id].boundary.boundarySetPlugin(job,
                                                           &(job->bodies[id]),
                                                           job->serializer.mainpath + job->bodies[id].boundary.filepath,
                                                           job->bodies[id].boundary.filename,
                                                           job->bodies[id].boundary.fp64_props,
                                                           job->bodies[id].boundary.int_props,
                                                           job->bodies[id].boundary.str_props);
                std::getline(fin, line); //filename
                if (0 == job->bodies[id].boundary.boundaryLoadState(job, &(job->bodies[id]), &(job->serializer),
                                                                    filepath + line)) {
                    std::cout << "Loading failed." << std::endl;
                    fin.close();
                    return 0;
                }
            } else {
                std::getline(fin, line); //filename
            }

            if (!std::getline(fin,line)){
                break;
            } //# contact
        }

        size_t contact_id = 0;
        while (line.compare(headers[6]) == 0){
            size_t id = contact_id++;

            loadStandardProps(&(job->contacts[id]),fin);
            if (!(job->contacts[id].filename.empty()) && job->activeContacts[id] != 0) {
                job->contacts[id].contactSetPlugin(job,
                                                   job->serializer.mainpath + job->contacts[id].filepath,
                                                   job->contacts[id].filename,
                                                   job->contacts[id].fp64_props,
                                                   job->contacts[id].int_props,
                                                   job->contacts[id].str_props);
                std::getline(fin, line); //filename
                if (0 == job->contacts[id].contactLoadState(job, &(job->serializer), filepath + line)) {
                    std::cout << "Loading failed." << std::endl;
                    fin.close();
                    return 0;
                }
            } else {
                std::getline(fin, line);
            }

            if (!std::getline(fin,line)){
                break;
            } //# contact
        }

        for (size_t c=0;c<job->contacts.size();c++){
            //contact
            if (job->activeContacts[c] == 0){
                continue;
            }
            job->contacts[c].contactSetFnPointers();
        }

        for (size_t b=0;b<job->bodies.size();b++){
            if (job->activeBodies[b] == 0){
                continue;
            }
            //boundary
            if(job->bodies[b].activeBoundary != 0) {
                job->bodies[b].boundary.boundarySetFnPointers();
            }
            //material
            if(job->bodies[b].activeMaterial != 0) {
                job->bodies[b].material.materialSetFnPointers();
            }
        }

        fin.close();

    } else {
        std::cout << "ERROR: Unable to open file: " << fullpath << std::endl;
        return 0;
    }

    std::cout << "Job properties (DIM = " << job->DIM << ", t = " << job->t << ", dt = " << job->dt << ")." << std::endl;
    std::cout << "Job Loaded." << std::endl;
    return 1;
}

/*----------------------------------------------------------------------------*/
