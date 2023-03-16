//
// Created by hemi-user on 2/21/23.
//

//
// Created by aaron on 12/12/19.
// finite_volume_driver.cpp
//

//
// Created by aaron on 5/10/18.
// default_driver.cpp
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <time.h>

#include "job.hpp"
#include "parser.hpp"
#include "mpm_objects.hpp"
#include "mpm_tensor.hpp"
#include "mpm_tensorarray.hpp"
#include "mpm_vector.hpp"
#include "mpm_vectorarray.hpp"

#include "fvm_objects.hpp"
#include "fvm_drivers.hpp"
#include "fvm_artificial_viscosity.hpp"
#include "registry.hpp"

/*----------------------------------------------------------------------------*/

void FVMBolideRestartDriver::init(Job* job){

    //call finite volume driver initializer
    FiniteVolumeDriver::init(job);

    if (TYPE != THERMAL){
        std::cerr << "FVMBolideRestartDriver requires THERMAL simulation. Exiting." << std::endl;
        exit(0);
    } else if (fluid_grid->object_name.compare("FVMGmsh3D") != 0){
        std::cerr << fluid_grid->object_name << " : FVMBolideRestartDriver requires FVMGmsh3D grid. Exiting." << std::endl;
        exit(0);
    }

    int GRID_DIM = 0;
    //assign grid dimension from job type
    if (job->JOB_TYPE == job->JOB_1D){
        GRID_DIM = 1;
    } else if (job->JOB_TYPE == job->JOB_2D){
        GRID_DIM = 2;
    } else if (job->JOB_TYPE == job->JOB_3D){
        GRID_DIM = 3;
    } else if (job->JOB_TYPE == job->JOB_2D_OOP){
        GRID_DIM = 2; //this is important, job->DIM =/= job->grid->GRID_DIM
    } else if (job->JOB_TYPE == job->JOB_AXISYM){
        GRID_DIM = 2; //this is important, job->DIM =/= job->grid->GRID_DIM
    } else {
        std::cerr << "Job doesn't have defined type for input " << job->JOB_TYPE << "." << std::endl;
        exit(0);
    }

    //check size of properties passed to driver object
    if (fp64_props.size() < 1+GRID_DIM+2) {
        //use default properties
        dt0 = job->dt;
        FVMBolideImpactDriver::H = 0;
        FVMBolideImpactDriver::V = 11e3;
        lambda = 0.5;
        eta = 0;
    } else if (fp64_props.size() < 1+GRID_DIM+4) {
        //assign initial height and velocity
        dt0 = job->dt;
        FVMBolideImpactDriver::H = fp64_props[1 + GRID_DIM];
        FVMBolideImpactDriver::V = fp64_props[2 + GRID_DIM];
        lambda = 0.5;
        eta = 0;
    } else {
        //assign initial height, velocity, and time-increment parameters
        dt0 = job->dt;
        FVMBolideImpactDriver::H = fp64_props[1 + GRID_DIM];
        FVMBolideImpactDriver::V = fp64_props[2 + GRID_DIM];
        lambda = fp64_props[3 + GRID_DIM];
        eta = fp64_props[4 + GRID_DIM];
    }

    if (str_props.size() < 2){
        std::cerr << "FVMBolideRestartDriver requires 2 input files in order to restart. " << str_props.size() << " given. Exiting." << std::endl;
        exit(0);
    } else {
        point_file = str_props[0];
        volume_file = str_props[1];
        std::cout << "FVMBolideRestartDriver restarting simulation using data from:\n";
        std::cout << "    Points: " << point_file << "\n";
        std::cout << "    Volumes: " << volume_file << "\n";
    }

    //loop over str-props and assign relevant flags
    std::vector<std::string> options = {"USE_ARTIFICIAL_VISCOSITY"};
    for (int i=0; i<str_props.size(); i++){
        switch (Parser::findStringID(options, str_props[i])){
            case 0:
                //USE_ARTIFICIAL_VISCOSITY
                USE_ARTIFICIAL_VISCOSITY = true;
                std::cout << "Using artificial viscosity." << std::endl;
                break;
            default:
                //do nothing
                break;
        }
    }

    // set up thermal floor
    USE_THERMAL_FLOOR = true;
    cv_floor = 717.63; T_floor = 100;

    std::cout << "Using variable time-stepping. Default: " << dt0 << "s, Mesh Quality: " << lambda << ", Viscosity: " << eta << std::endl;
    std::cout << "Bolide Impact Simulation:\n";
    std::cout << "    Inclination: 45 Degrees\n";
    std::cout << "    Initial Velocity: " << FVMBolideImpactDriver::V << " m/s\n";
    std::cout << "    Initial Height: " << FVMBolideImpactDriver::H << " m\n";
    std::cout << "    Initial Temperature: " << getAmbientTemperature(FVMBolideImpactDriver::H) << " K\n";
    std::cout << "    Initial Pressure: " << getAmbientPressure(FVMBolideImpactDriver::H) << " Pa\n";
    std::cout << "    Initial Density: " << getAmbientDensity(FVMBolideImpactDriver::H) << " kg/m^3" << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/

void FVMBolideRestartDriver::restart(Job* job) {
    // CODE TO RESET SIMULATION CONDITIONS ACCORDING TO INPUT DATA
    Points* points = job->bodies[0]->points.get();

    // 0. Open Points and Volumes Files
    std::ifstream pin(point_file);
    if (!pin.is_open()){
        std::cerr << "ERROR: Unable to open point restart file: " << point_file << ". Exiting." << std::endl;
        exit(0);
    }
    std::ifstream vin(volume_file);
    if (!vin.is_open()){
        std::cerr << "ERROR: Unable to open volume restart file: " << point_file << ". Exiting." << std::endl;
        exit(0);
    }

    // Write Update to Console
    std::cout << "Restarting Simulation." << std::endl;

    // 1. Read Point Data
    std::string line;
    std::vector<std::string> svec;
    int result = -1;
    int plen = 0;
    if (pin.is_open()){
        // Look For "POINTS" Label
        while (std::getline(pin,line)){
            //Split Line Using " "
            svec = Parser::splitString(line,' ');

            //Check for POINTS Label
            if (svec.size() > 0 && svec[0].compare("POINTS") == 0){
                //POINTS Label Found!
                //Next Value is Number of Points
                plen = std::stod(svec[1]);

                //Exit While Loop
                break;
            }
        }

        // Read in Data
        if (plen > 0){

            // Write Update to Console
            std::cout << "Reading Point Data." << std::endl;

            //resize KinematicVectors
            points->x = KinematicVectorArray(plen, job->JOB_TYPE);
            points->u = KinematicVectorArray(plen, job->JOB_TYPE);
            points->x_t = KinematicVectorArray(plen, job->JOB_TYPE);
            points->mx_t = KinematicVectorArray(plen, job->JOB_TYPE);
            points->b = KinematicVectorArray(plen, job->JOB_TYPE);

            //resize scalar vectors
            points->m.resize(plen);
            points->v.resize(plen);
            points->active.resize(plen);

            //resize tensor arrays
            points->T = MaterialTensorArray(plen);
            points->L = KinematicTensorArray(plen, job->JOB_TYPE);

            //zero out all entries to start
            points->x.setZero();
            points->u.setZero();
            points->x_t.setZero();
            points->m.setZero();
            points->v.setZero();
            points->mx_t.setZero();
            points->b.setZero();
            points->T.setZero();
            points->L.setZero();
            points->active.setZero();

            // Read in Position Data
            for (int i=0; i<plen; i++){
                // Read Line
                std::getline(pin,line);

                // Split Line
                svec = Parser::splitString(line,' ');

                // Check For Correct Line Length
                if (svec.size() != 3){
                    std::cerr << "UH OH! Point position data does not appear formatted correctly...\n";
                    std::cerr << "    [" << i << "]: " << line << std::endl;
                    std::cerr << "Exiting.\n";
                    exit(0);
                }

                // Assign Position Data
                points->x(i,0) = std::stod(svec[0]);
                points->x(i,1) = std::stod(svec[1]);
                points->x(i,2) = std::stod(svec[2]);
            }

            // Write Update to Console
            std::cout << "    Position Data Complete..." << std::endl;

            // Read Through Lines Until "POINT_DATA"
            while (std::getline(pin,line)){
                //Split Line Using " "
                svec = Parser::splitString(line,' ');

                //Check for POINTS Label
                if (svec.size() > 0 && svec[0].compare("POINT_DATA") == 0){
                    //POINT_DATA Label Found!
                    //Exit While Loop
                    break;
                }
            }

            // Write Update to Console
            std::cout << "    Found POINT_DATA.\n";

            // Keywords
            std::vector<std::string> scalar_tags = {"mass","volume","active"};
            std::vector<std::string> vector_tags = {"displacement","velocity","momentum"};
            std::vector<std::string> tensor_tags = {"cauchy_stress"};

            // Read Through Lines Until SCALARS, VECTORS, or TENSORS Found
            while (std::getline(pin,line)){
                //Split Line Using " "
                svec = Parser::splitString(line,' ');

                //Check for Data Labels
                if (svec.size() > 0 && svec[0].compare("SCALARS") == 0){
                    //Found SCALARS Label
                    //Next Value is Data Name
                    result = -1;
                    result = Parser::findStringID(scalar_tags, svec[1]);

                    //Read Next Line
                    std::getline(pin,line);

                    //Read Data Into Correct Array
                    for (int i=0; i<plen; i++) {
                        //Read Data Line
                        std::getline(pin,line);

                        //SCALARS Only Have 1 Value
                        switch (result) {
                            case 0:
                                //mass
                                points->m(i) = std::stod(line);
                                break;
                            case 1:
                                //volume
                                points->v(i) = std::stod(line);
                                break;
                            case 2:
                                //active
                                points->active(i) = std::stoi(line);
                                break;
                            default:
                                //do nothing
                                break;
                        }
                    }

                    if (result >= 0 && result < scalar_tags.size()) {
                        //Write Update to Console
                        std::cout << "    Read In " << scalar_tags[result] << " Data...\n";
                    }

                } else if (svec.size() > 0 && svec[0].compare("VECTORS") == 0){
                    //Found VECTORS Label
                    //Next Value is Data Name
                    result = -1;
                    result = Parser::findStringID(vector_tags, svec[1]);

                    //Read Data Into Correct Array
                    for (int i=0; i<plen; i++) {
                        //Read Data Line
                        std::getline(pin,line);

                        //VECTORS Have 3 Values
                        svec = Parser::splitString(line,' ');
                        for (int j=0; j<3; j++) {
                            switch (result) {
                                case 0:
                                    //displacement
                                    points->u(i,j) = std::stod(svec[j]);
                                    break;
                                case 1:
                                    //velocity
                                    points->x_t(i,j) = std::stod(svec[j]);
                                    break;
                                case 2:
                                    //momentum
                                    points->mx_t(i,j) = std::stod(svec[j]);
                                    break;
                                default:
                                    //do nothing
                                    break;
                            }
                        }
                    }

                    if (result >= 0 && result < vector_tags.size()) {
                        //Write Update to Console
                        std::cout << "    Read In " << vector_tags[result] << " Data...\n";
                    }

                } else if (svec.size() > 0 && svec[0].compare("TENSORS") == 0){
                    //Found TENSORS Label
                    //Next Value is Data Name
                    result = -1;
                    result = Parser::findStringID(tensor_tags, svec[1]);

                    //Read Data Into Correct Array
                    for (int i=0; i<plen; i++) {
                        //Read Data Line
                        std::getline(pin,line);

                        //TENSORS Have 3 Sets of 3 Values
                        svec = Parser::splitString(line,' ');
                        for (int j=0; j<3; j++) {
                            switch (result) {
                                case 0:
                                    //cauchy_stress
                                    points->T(i,j) = std::stod(svec[j]);
                                    break;
                                default:
                                    //do nothing
                                    break;
                            }
                        }
                        std::getline(pin,line);
                        svec = Parser::splitString(line,' ');
                        for (int j=0; j<3; j++) {
                            switch (result) {
                                case 0:
                                    //cauchy_stress
                                    points->T(i,3+j) = std::stod(svec[j]);
                                    break;
                                default:
                                    //do nothing
                                    break;
                            }
                        }
                        std::getline(pin,line);
                        svec = Parser::splitString(line,' ');
                        for (int j=0; j<3; j++) {
                            switch (result) {
                                case 0:
                                    //cauchy_stress
                                    points->T(i,6+j) = std::stod(svec[j]);
                                    break;
                                default:
                                    //do nothing
                                    break;
                            }
                        }

                        //TENSORS Also Has Additional Line Break
                        std::getline(pin,line);
                    }

                    if (result >= 0 && result < tensor_tags.size()) {
                        //Write Update to Console
                        std::cout << "    Read In " << tensor_tags[result] << " Data...\n";
                    }
                }
            }

        } else {
            std::cerr << "ERROR: Unable to find POINTS label in " << point_file << ". Exiting." << std::endl;
            exit(0);
        }

        //Write Update to Console
        std::cout << "    Complete.\n";

        //Close File
        pin.close();
    }

    // 2. Read Volume Data
    int vlen = 0;
    if (vin.is_open()){
        // Look For "CELL_DATA" Label
        while (std::getline(vin,line)){
            //Split Line Using " "
            svec = Parser::splitString(line,' ');

            //Check for CELL_DATA Label
            if (svec.size() > 0 && svec[0].compare("CELL_DATA") == 0){
                //CELL_DATA Label Found!
                //Next Value is Number of Cells
                vlen = std::stod(svec[1]);

                //Exit While Loop
                break;
            }
        }

        // Read in Data
        if (vlen == fluid_grid->element_count){

            // Write Update to Console
            std::cout << "Reading Volume Data." << std::endl;

            // Write Update to Console
            std::cout << "    Found CELL_DATA.\n";

            // Keywords
            std::vector<std::string> scalar_tags = {"density", "energy", "porosity"};
            std::vector<std::string> vector_tags = {"momentum"};

            // Read Through Lines Until SCALARS or VECTORS Found
            while (std::getline(vin,line)){
                //Split Line Using " "
                svec = Parser::splitString(line,' ');

                //Check for Data Labels
                if (svec.size() > 0 && svec[0].compare("SCALARS") == 0){
                    //Found SCALARS Label
                    //Next Value is Data Name
                    result = -1;
                    result = Parser::findStringID(scalar_tags, svec[1]);

                    //Read Next Line
                    std::getline(vin,line);

                    //Read Data Into Correct Array
                    for (int i=0; i<vlen; i++) {
                        //Read Data Line
                        std::getline(vin,line);

                        //SCALARS Only Have 1 Value
                        switch (result) {
                            case 0:
                                //density
                                fluid_body->rho(i) = std::stod(line);
                                break;
                            case 1:
                                //energy
                                fluid_body->rhoE(i) = std::stod(line);
                                break;
                            case 2:
                                //porosity
                                fluid_body->n_e(i) = std::stoi(line);
                                break;
                            default:
                                //do nothing
                                break;
                        }
                    }

                    if (result >= 0 && result < scalar_tags.size()) {
                        //Write Update to Console
                        std::cout << "    Read In " << scalar_tags[result] << " Data...\n";
                    }

                } else if (svec.size() > 0 && svec[0].compare("VECTORS") == 0){
                    //Found VECTORS Label
                    //Next Value is Data Name
                    result = -1;
                    result = Parser::findStringID(vector_tags, svec[1]);

                    //Read Data Into Correct Array
                    for (int i=0; i<vlen; i++) {
                        //Read Data Line
                        std::getline(vin,line);

                        //VECTORS Have 3 Values
                        svec = Parser::splitString(line,' ');
                        for (int j=0; j<3; j++) {
                            switch (result) {
                                case 0:
                                    //momentum
                                    fluid_body->p(i,j) = std::stod(svec[j]);
                                    break;
                                default:
                                    //do nothing
                                    break;
                            }
                        }
                    }

                    if (result >= 0 && result < vector_tags.size()) {
                        //Write Update to Console
                        std::cout << "    Read In " << vector_tags[result] << " Data...\n";
                    }

                }
            }

        } else {
            std::cerr << "ERROR: Unable to find CELL_DATA label in " << volume_file << ". Exiting." << std::endl;
            exit(0);
        }

        //Write Update to Console
        std::cout << "    Complete.\n";

        //Close File
        pin.close();
    }

    return;
}


/*----------------------------------------------------------------------------*/

void FVMBolideRestartDriver::run(Job* job) {
    //initialize FVM objects
    serializer->init(job, this);
    fluid_grid->init(job, this);
    solver->init(job, this);
    fluid_material->init(job, this);
    fluid_body->init(job, this);

    //attempt restart
    restart(job);

    //set counters to zero
    int stepCount = 0;
    int frameCount = 0;

    struct timespec timeStart, timeFrame, timeFinish, timeStep;
    clock_gettime(CLOCK_MONOTONIC, &timeStart);
    timeFrame = timeStart;
    timeStep = timeStart;
    //clock_t clockSim = clock();
    //clock_t clockFrame = clock();
    double tSim = 0;
    double tFrame = 0;
    double tStep = 0;

    //initialize gravity
    generateGravity(job);
    applyGravity(job);

    //bolide simulation variables
    double dV = 0;  //total change in velocity
    double KE = 0;  //total kinetic energy (in Earth Ref. Frame)
    double M = 0;   //total bolide mass
    double tmpKE, tmpt, KE0;
    double dKE = 0; //total change in kinetic energy
    double dKEdt = 0; //rate of change in kinetic energy
    double V0 = FVMBolideImpactDriver::V;  //initial velocity
    double H0 = FVMBolideImpactDriver::H;  //initial height
    double cosTHETA = std::sqrt(0.5); //cosine of inclination angle
    double YMin = FVMBolideImpactDriver::Y0;  //minimum bolide point y-position
    double YMax = FVMBolideImpactDriver::Y0;  //maximum bolide point y-position (for radius calculation)
    double YDot = 0;   //y-velocity of leading bolide point
    double YTmp = 0;
    int YMinID  = 0;   //ID of point with minimum y-position

    bool YMin_Set = false;
    for (int i=0; i<job->bodies[0]->points->x.size(); i++){
        // Only Consider Active Points
        if (job->bodies[0]->points->active(i)){
            if (!YMin_Set){
                //Initialize YMin
                YMin = job->bodies[0]->points->x(i,1);
                YMax = job->bodies[0]->points->x(i,1);
                YMin_Set = true;
            } else {
                YTmp = job->bodies[0]->points->x(i,1);
                if (YTmp < YMin){
                    YMin = YTmp;
                } else if (YTmp > YMax){
                    YMax = YTmp;
                }
            }
        }
    }
    FVMBolideImpactDriver::Y0 = YMin;

    // Compute Radius and Centroid of Bolide
    double R0 = (YMax - YMin) / 2.0;
    double YMed = (YMax + YMin) / 2.0;
    KinematicVector X0 = KinematicVector(job->JOB_TYPE);
    KinematicVector XTmp = KinematicVector(job->JOB_TYPE);
    X0.setZero(); X0(1) = YMed;

    // Set Initial Conditions
    if (true){
        for (int e=0; e<fluid_grid->element_count; e++){
            XTmp = fluid_grid->getElementCentroid(job, e);
            // Check That Element is Outside of Bolide
            if ((XTmp - X0).dot(XTmp - X0) > 2.25*R0*R0) {
                // Subtract Kinetic Energy From Total Energy
                fluid_body->rhoE(e) -= 0.5 * fluid_body->p(e).dot(fluid_body->p(e)) / fluid_body->rho(e);
                // Update Momentum
                fluid_body->p(e, 1) = fluid_body->rho(e) * FVMBolideImpactDriver::V;
                // Add Kinetic Energy To Total Energy
                fluid_body->rhoE(e) += 0.5 * fluid_body->p(e).dot(fluid_body->p(e)) / fluid_body->rho(e);
            }
        }
    }

    // Compute Initial Kinetic Energy
    KinematicVector vTMP = KinematicVector(job->JOB_TYPE);
    for (int i=0; i<job->bodies[0]->points->x.size(); i++){
        if (job->bodies[0]->points->active(i) == 1) {
            vTMP = job->bodies[0]->points->x_t[i];
            vTMP[1] -= FVMBolideImpactDriver::V;
            KE += 0.5 * job->bodies[0]->points->m(i) * vTMP.dot(vTMP);
            M += job->bodies[0]->points->m(i);
        }
    }
    KE0 = KE;
    tmpt = job->t;


    //create new output file
    //open and clear file
    std::ostringstream oss;
    oss << "Bolide_History_" << (FVMBolideImpactDriver::V/1e3) << "kms_" << (FVMBolideImpactDriver::H/1e3) << "km.csv";
    std::string output_filename = oss.str();
    std::ofstream file (output_filename,std::ios::trunc);
    if (file.is_open()){
        //success!
        //write file header
        file << "Time [s], Time Increment [s], Velocity [km/s], Altitude [km], Pressure [Pa], Temperature [K], Density [kg/m^3], Leading Edge Position [m], Kinetic Energy [J], Change in Velocity [km/s], Change in Kinetic Energy [J], Rate of Change of Kinetic Energy [J/s], Bolide Mass [kg]\n";
        file << job->t << ", ";
        file << job->dt << ", ";
        file << (FVMBolideImpactDriver::V/1e3) << ", ";
        file << (FVMBolideImpactDriver::H/1e3) << ", ";
        file << getAmbientPressure(FVMBolideImpactDriver::H) << ", ";
        file << getAmbientTemperature(FVMBolideImpactDriver::H) << ", ";
        file << getAmbientDensity(FVMBolideImpactDriver::H) << ", ";
        file << YMin << ", ";
        file << KE << ", ";
        file << (dV/1e3) << ", ";
        file << dKE << ", ";
        file << dKEdt << ", ";
        file << M << "\n";
        file.close();
    } else {
        std::cerr << "ERROR! Cannot open " << output_filename << "! Exiting." << std::endl;
        exit(0);
    }

    //element length scales
    Eigen::VectorXd l(fluid_grid->element_count);
    l.setZero();
    std::vector<int> e_faces;
    double l_min, l_e;
    for (int e=0; e<fluid_grid->element_count; e++){
        e_faces = fluid_grid->getElementFaces(e);
        for (int f=0; f<e_faces.size(); f++){
            l_e = (fluid_grid->getElementCentroid(job, e) - fluid_grid->getFaceCentroid(job, e_faces[f])).norm();
            if (f==0 || l_e < l_min){
                l_min = l_e;
            }
        }
        l(e) = l_min;
    }

    //temporary velocity magnitude and stable time step
    double c, v, dts;

    //artificial viscosity
    std::vector<ArtificialViscosityCalculator::fluxVector> fluxVectors;
    if (USE_ARTIFICIAL_VISCOSITY){
        fluxVectors.resize(fluid_grid->face_count);
    }

    //reality check
    bool failed = false;

    //run simulation until stop_time
    while (job->t <= stop_time){

        //FIRST: Identify Leading Material Point
        bool YMin_Set = false;
        YMinID = 0;
        for (int i=0; i<job->bodies[0]->points->x.size(); i++){
            // Only Consider Active Points
            if (job->bodies[0]->points->active(i)){
                if (!YMin_Set){
                    //Initialize YMin
                    YMin = job->bodies[0]->points->x(i,1);
                    YMinID = i;
                    YMin_Set = true;
                } else {
                    YTmp = job->bodies[0]->points->x(i,1);
                    if (YTmp < YMin){
                        YMin = YTmp;
                        YMinID = i;
                    }
                }
            }
        }

        //SECOND: Determine Velocity of Leading Material Point
        if (YMin_Set) {
            YDot = job->bodies[0]->points->x_t(YMinID, 1);
        }

        //THIRD: Adjust Reference Frame of Solution Fields to New Velocity
        if (YMin_Set){
            // (1) Adjust Material Point Velocities and Momenta
            for (int i=0; i<job->bodies[0]->points->x.size(); i++){
                if (job->bodies[0]->points->active(i)){
                    job->bodies[0]->points->x_t(i,1) -= YDot;
                    job->bodies[0]->points->mx_t(i) = job->bodies[0]->points->m(i) * job->bodies[0]->points->x_t(i);
                }
            }

            // (2) Adjust Fluid Momenta and Energies
            for (int e=0; e<fluid_grid->element_count; e++){
                // Subtract Kinetic Energy From Total Energy
                fluid_body->rhoE(e) -= 0.5 * fluid_body->p(e).dot(fluid_body->p(e)) / fluid_body->rho(e);
                // Update Momentum
                fluid_body->p(e,1) -= fluid_body->rho(e) * YDot;
                // Add Kinetic Energy To Total Energy
                fluid_body->rhoE(e) += 0.5 * fluid_body->p(e).dot(fluid_body->p(e)) / fluid_body->rho(e);
            }

            // (3) Adjust Velocity BC
            FVMBolideImpactDriver::V -= YDot;
            dV -= YDot;
        }


        //FOURTH: check time increment for stability
        job->dt = dt0;
        for (int e=0; e<fluid_grid->element_count; e++){
            //check advection
            v = fluid_body->p(e).norm()/fluid_body->rho(e);
            c = fluid_material->getSpeedOfSound(job,
                                                this,
                                                fluid_body->rho(e),
                                                fluid_body->p(e),
                                                fluid_body->rhoE(e),
                                                fluid_body->n_e(e));
            dts = lambda * l(e)/(c+v);
            if (dts < job->dt){
                job->dt = dts;
            }

            //check viscosity
            if (eta > 0){
                dts = (lambda * lambda * l(e) * l(e) * fluid_body->rho(e)) / eta;
                if (dts < job->dt){
                    job->dt = dts;
                }
            }

            //check thermal floor
            if (TYPE == THERMAL && USE_THERMAL_FLOOR && fluid_material->getTemperature(job,
                                                                                       this,
                                                                                       fluid_body->rho(e),
                                                                                       fluid_body->p(e),
                                                                                       fluid_body->rhoE(e),
                                                                                       fluid_body->n_e(e)) <= T_floor){

                //Highlight Error
                std::cout << std::endl;
                std::cout << "    Thermal Floor: [" << e << "]: " << fluid_material->getTemperature(job,
                                                                                                    this,
                                                                                                    fluid_body->rho(e),
                                                                                                    fluid_body->p(e),
                                                                                                    fluid_body->rhoE(e),
                                                                                                    fluid_body->n_e(e));
                std::cout << " --> " << T_floor << std::endl;

                //Assign New Temperature
                fluid_body->rhoE(e) = fluid_body->rho(e) * cv_floor * T_floor
                                      + 0.5*fluid_body->p[e].dot(fluid_body->p[e])/fluid_body->rho(e);

            }

            if (!failed && (fluid_material->getTemperature(job,
                                                           this,
                                                           fluid_body->rho(e),
                                                           fluid_body->p(e),
                                                           fluid_body->rhoE(e),
                                                           fluid_body->n_e(e)) <= 0
                            ||
                            std::isnan(fluid_material->getTemperature(job,
                                                                      this,
                                                                      fluid_body->rho(e),
                                                                      fluid_body->p(e),
                                                                      fluid_body->rhoE(e),
                                                                      fluid_body->n_e(e))))){
                failed = true;
                std::cout << std::endl << std::endl;
                std::cout << "    First failue: [" << e << "]" << std::endl;
                std::cout << "    Density: " << fluid_body->rho(e) << std::endl;
                std::cout << "    Momentum: " << fluid_body->p(e).norm() << std::endl;
                std::cout << "    Energy: " << fluid_body->rhoE(e)<< std::endl;
                std::cout << "    Porosity: " << fluid_body->n_e(e)<< std::endl;
                std::cout << "    Position: " << EIGEN_MAP_OF_KINEMATIC_VECTOR(fluid_grid->getElementCentroid(job,e)).transpose();
                std::cout << std::endl;
                std::cout << "Exiting." << std::endl;
                exit(0);
            }
        }
        //report step length to console
        if (job->dt < dt0){
            std::cout << "dt: " << job->dt << "    \r" << std::flush;
        }


        // FIFTH: run solver
        solver->step(job,this);

        // SIXTH: Write Data to Files
        if (job->serializer->writeFrame(job) == 1) {
            //call fvm serializer to write frame as well:
            serializer->writeFrame(job,this);

            //compute KE, dKE, dKEdt
            tmpKE = KE; KE = 0; M = 0;
            for (int i=0; i<job->bodies[0]->points->x.size(); i++){
                if (job->bodies[0]->points->active(i) == 1) {
                    vTMP = job->bodies[0]->points->x_t[i];
                    vTMP[1] -= FVMBolideImpactDriver::V;
                    KE += 0.5 * job->bodies[0]->points->m(i) * vTMP.dot(vTMP);
                    M += job->bodies[0]->points->m(i);
                }
            }
            dKE = KE - KE0;
            dKEdt = (KE - tmpKE) / (job->t + job->dt - tmpt);
            tmpt = job->t + job->dt;

            //write bolide simulation data
            file.open(output_filename,std::ios::app);
            if (file.is_open()){
                //success!
                file << job->t + job->dt << ", ";
                file << job->dt << ", ";
                file << (FVMBolideImpactDriver::V/1e3) << ", ";
                file << (FVMBolideImpactDriver::H/1e3) << ", ";
                file << getAmbientPressure(FVMBolideImpactDriver::H) << ", ";
                file << getAmbientTemperature(FVMBolideImpactDriver::H) << ", ";
                file << getAmbientDensity(FVMBolideImpactDriver::H) << ", ";
                file << YMin  << ", ";
                file << KE << ", ";
                file << (dV/1e3) << ", ";
                file << dKE << ", ";
                file << dKEdt << ", ";
                file << M << "\n";
                file.close();
            } else {
                std::cerr << "ERROR! Cannot open " << output_filename << "! Exiting." << std::endl;
                exit(0);
            }

            //successful frame written
            //tFrame = (double)(clock() - clockFrame)/CLOCKS_PER_SEC;
            //tSim = (double)(clock() - clockSim)/CLOCKS_PER_SEC;
            clock_gettime(CLOCK_MONOTONIC,&timeFinish);
            tFrame = (timeFinish.tv_sec - timeFrame.tv_sec) + (timeFinish.tv_nsec - timeFrame.tv_nsec)/1000000000.0;
            tSim = (timeFinish.tv_sec - timeStart.tv_sec) + (timeFinish.tv_nsec - timeStart.tv_nsec)/1000000000.0;
            timeFrame = timeFinish;
            printf("\33[2K");
            std::cout << "Frame Written [" << ++frameCount << "]. Time/Frame [" << tFrame << " s]. Elapsed Time [" << tSim << " s]." << std::flush;
        }
        std::cout << "\r";

        //SEVENTH: apply artificial viscosity
        if (USE_ARTIFICIAL_VISCOSITY){
            applyArtificialViscosityFluxes(job, fluxVectors);
        }

        //EIGHTH: update altitude
        FVMBolideImpactDriver::H -= cosTHETA * FVMBolideImpactDriver::V * job->dt;

        //NINTH: update time
        job->t += job->dt;
    }
    //tSim = (double)(clock() - clockSim)/CLOCKS_PER_SEC;
    clock_gettime(CLOCK_MONOTONIC,&timeFinish);
    tSim = (timeFinish.tv_sec - timeStart.tv_sec) + (timeFinish.tv_nsec - timeStart.tv_nsec)/1000000000.0;
    std::cout << std::endl << std::endl << "Simulation Complete. Elapsed Time [" << tSim << "s]." << std::endl;
    return;
}
