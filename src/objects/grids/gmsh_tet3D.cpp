//
// Created by aaron on 11/30/18.
// gmsh_tet3D.cpp
//

//
// Created by aaron on 5/25/18.
// gmsh_tri2D.cpp
//


#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <math.h>

#include "mpm_objects.hpp"
#include "mpm_vector.hpp"
#include "mpm_tensor.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensorarray.hpp"
#include "mpm_sparse.hpp"
#include "job.hpp"

#include "grids.hpp"

/*----------------------------------------------------------------------------*/
//
int TetrahedralGridLinear::whichSearchCell(const KinematicVector &xIN){
    assert(xIN.VECTOR_TYPE == Lx.VECTOR_TYPE && "whichSearchCell failed; invalid KinematicVector TYPE");
    bool inDomain;
    KinematicVector xTMP = xIN;

    //adjust bounding nodes to avoid overflow
    for (int i=0; i<GRID_DIM; i++){
        if (xIN[i] == Lxmax[i]){
            xTMP[i] -= hx[i]/100.0;
        }
        //rezero xtmp
        xTMP[i] -= Lxmin[i];
    }

    int elementID = floor(xTMP[0] / hx[0]); //id in x-dimension
    for (int i = 0; i < GRID_DIM; i++) {
        inDomain = (xTMP[i] < Lx[i] && xTMP[i] >= 0);
        if (!inDomain) { //if xIn is outside domain, return -1
            return -1;
        }
        if (i == 1) {
            //add number of elements in next layer of id in higher dimensions
            elementID += floor(xTMP[i] / hx[i]) * (Nx[i - 1]);
        } else if (i == 2) {
            elementID += floor(xTMP[i] / hx[i]) * (Nx[i - 1]) * (Nx[i - 2]);
        }
    }
    return elementID;
}

bool TetrahedralGridLinear::inElement(Job* job, const KinematicVector& xIN, int idIN){
    //check if point is interior to element
    KinematicVector xi = KinematicVector(job->JOB_TYPE);
    KinematicTensor tmpMat = KinematicTensor(job->JOB_TYPE);

    //xi = Ainv * (x - x0)
    xi = Ainv(idIN)*(xIN - x_n(nodeIDs(idIN,0)));

    return (xi(0) >= -1e-10 && xi(1) >= -1e-10 && xi(2) >= -1e-10 && (xi(0) + xi(1) + xi(2)) <= (1. + 1e-10));
}


/*----------------------------------------------------------------------------*/
//
void TetrahedralGridLinear::init(Job* job){
    if (GRID_DIM != 3){
        std::cerr << "TetrahedralGridLinear requires GRID_DIM = 3, got: " << GRID_DIM << std::endl;
        exit(0);
    }

    if (fp64_props.size() < 1|| str_props.size() < 1){
        std::cout << fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need 1 filename and 1 length scale defined.\n",
                __FILE__, __func__);
        exit(0);
    } else {
        //store length, number of linear nodes, and deltas
        msh_filename = str_props[0];

        //use node to element map instead of search grid
        if (int_props.size() > 0 && int_props[0] == 1){
            use_n_to_e = true;
        }

        if (fp64_props.size() == 3) {
            //assume hx, hy, hz given
            hx = KinematicVector(job->JOB_TYPE);
            for (int pos=0;pos<GRID_DIM;pos++){
                hx(pos) = fp64_props[pos];
            }
            lc = -1.; //flag for not using lc

            std::cout << "Grid properties (filename = " << msh_filename << ", hx = " << EIGEN_MAP_OF_KINEMATIC_VECTOR(hx).transpose() << ")." << std::endl;
        } else {
            lc = fp64_props[0];

            std::cout << "Grid properties (filename = " << msh_filename << ", lc = " << lc << ")." << std::endl;
        }
    }

    int len, tag;
    std::string line; //read line
    std::vector<std::string> lvec;
    std::ifstream fin(msh_filename); //file to load from
    if (fin.is_open()) {
        //if open, read lines
        std::vector<std::string> headers= {"$MeshFormat","$Nodes","$Elements"};
        while (std::getline(fin,line)) {
            switch(Parser::findStringID(headers,line)){
                case 0:
                    //mesh format
                    std::getline(fin,line);
                    lvec = Parser::splitString(line,' ');
                    //mesh file version
                    msh_version = std::stod(lvec[0]);
                    break;
                case 1:
                    //nodes
                    std::getline(fin,line);
                    
                    if (floor(msh_version) == 2) {
                        len = std::stoi(line); //number of nodes;
                        //nodeTags.resize(len);
                        x_n = KinematicVectorArray(len,job->JOB_TYPE);
                        node_count = len;
                        nodeTags = Eigen::VectorXi(len);
                        nodeTags.setConstant(-1);

                        //loop to read in all nodes
                        for (int i=0;i<len;i++){
                            std::getline(fin,line);
                            lvec = Parser::splitString(line,' ');
                            //lvec[0] gives gmsh id (1-indexed)
                            x_n(i,0) = std::stod(lvec[1]); //x-coord
                            x_n(i,1) = std::stod(lvec[2]); //y-coord
                            x_n(i,2) = std::stod(lvec[3]); //z-coord
                        }
                    } else if (floor(msh_version) == 4){
                        lvec = Parser::splitString(line,' '); //entity_blocks num_nodes
                        len = std::stoi(lvec[1]);
                        
                        //nodeTags.resize(len);
                        x_n = KinematicVectorArray(len,job->JOB_TYPE);
                        node_count = len;
                        nodeTags = Eigen::VectorXi(len);
                        nodeTags.setConstant(-1);

                        //version 4 nodes are broken into blocks
                        int k = 0;
                        int block_length;
                        while (k < len){
                            //read line
                            std::getline(fin,line);
                            lvec = Parser::splitString(line, ' ');
                            //length of block
                            block_length = std::stoi(lvec[3]);
                            //skip tags
                            for (int i=0;i<block_length;i++){
                                std::getline(fin,line);
                            }
                            //read nodes in block
                            for (int i=0;i<block_length;i++){
                                std::getline(fin,line);
                                lvec = Parser::splitString(line,' ');
                                //lvec[0] gives gmsh id (1-indexed)
                                x_n(k,0) = std::stod(lvec[0]); //x-coord
                                x_n(k,1) = std::stod(lvec[1]); //y-coord
                                x_n(k,2) = std::stod(lvec[2]); //z-coord
                                //increment counter
                                k++;
                            }
                        }
                    } else {
                        std::cerr << "Unrecognized Gmsh version: " << msh_version << ". Exiting." << std::endl;
                        exit(0);
                    }
                    break;
                case 2:
                    //elements
                    if (floor(msh_version) == 2) {
                        std::getline(fin, line);
                        len = std::stoi(line); //number of potential elements

                        nodeIDs.resize(len, npe); //element to node map

                        int i = -1;
                        while (std::getline(fin, line)) {
                            if (line.compare("$EndElements") == 0) {
                                break;
                            }

                            lvec = Parser::splitString(line, ' ');
                            //check that element type is tetrahedron
                            if (std::stoi(lvec[1]) == 4) {
                                i++; //increment i
                                len = std::stoi(lvec[2]); //num tags
                                nodeIDs(i, 0) = std::stoi(lvec[3 + len]) - 1;
                                nodeIDs(i, 1) = std::stoi(lvec[4 + len]) - 1;
                                nodeIDs(i, 2) = std::stoi(lvec[5 + len]) - 1;
                                nodeIDs(i, 3) = std::stoi(lvec[6 + len]) - 1;
                            }

                            //if element type is triangle, set node tags according to line tag
                            if (std::stoi(lvec[1]) == 2){
                                len = std::stoi(lvec[2]);           //num tags
                                tag = std::stoi(lvec[3]);           //first tag
                                nodeTags(std::stoi(lvec[3 + len]) - 1) = tag;
                                nodeTags(std::stoi(lvec[4 + len]) - 1) = tag;
                                nodeTags(std::stoi(lvec[5 + len]) - 1) = tag;
                            }
                        }

                        element_count = i + 1;
                        nodeIDs.conservativeResize(i + 1, npe);

                    } else if (floor(msh_version) == 4){
                        std::getline(fin, line);
                        lvec = Parser::splitString(line,' ');
                        len = std::stoi(lvec[1]); //number of potential entities

                        nodeIDs.resize(len, npe); //element to node map
                        int i = -1;
                        int block_length, block_type;
                        while (std::getline(fin, line)){
                            if (line.compare("$EndElements") == 0){
                                break;
                            }

                            lvec = Parser::splitString(line, ' ');
                            //block info
                            tag = std::stoi(lvec[1]);
                            block_type = std::stoi(lvec[2]);
                            block_length = std::stoi(lvec[3]);
                            if (block_type == 4) {
                                for (int k = 0; k < block_length; k++) {
                                    std::getline(fin, line);
                                    lvec = Parser::splitString(line, ' ');

                                    i++; //increment i
                                    nodeIDs(i, 0) = std::stoi(lvec[1]) - 1;
                                    nodeIDs(i, 1) = std::stoi(lvec[2]) - 1;
                                    nodeIDs(i, 2) = std::stoi(lvec[3]) - 1;
                                    nodeIDs(i, 3) = std::stoi(lvec[4]) - 1;
                                }
                            } else if (block_type == 2) {
                                //if entity is a triangle, it defines a set of physical tags for nodes
                                for (int k=0; k<block_length; k++){
                                    std::getline(fin, line);
                                    lvec = Parser::splitString(line, ' ');
                                    nodeTags(std::stoi(lvec[1]) - 1) = tag;
                                    nodeTags(std::stoi(lvec[2]) - 1) = tag;
                                    nodeTags(std::stoi(lvec[3]) - 1) = tag;
                                }
                            } else {
                                //read through block and continue
                                for (int k=0; k<block_length; k++){
                                    std::getline(fin, line);
                                }
                            }

                        }
                    element_count = i + 1;
                    nodeIDs.conservativeResize(i + 1, npe);

                    } else {
                        std::cerr << "Unrecognized Gmsh version: " << msh_version << ". Exiting." << std::endl;
                        exit(0);
                    }
                    break;
                default:
                    //do nothing
                    break;
            }
        }
    } else {
        std::cerr << "ERROR: Cannot open " << msh_filename << "! Exiting." << std::endl;
        exit(0);
    }

    //call initializing function
    hiddenInit(job);

    std::cout << "Grid Initialized." << std::endl;

    return;
}
void TetrahedralGridLinear::hiddenInit(Job* job){

    //called from loading or initializing functions

    v_e.resize(nodeIDs.rows());
    v_n.resize(x_n.size());
    v_n.setZero();

    //element volume
    KinematicVector a, b, c;
    for (int e=0; e<nodeIDs.rows(); e++){
        a = x_n(nodeIDs(e,1)) - x_n(nodeIDs(e,0));
        b = x_n(nodeIDs(e,2)) - x_n(nodeIDs(e,0));
        c = x_n(nodeIDs(e,3)) - x_n(nodeIDs(e,0));

        v_e(e) = (a.dot((b.cross(c))))/6.0;

        //nodal volume
        for (int c=0;c<nodeIDs.cols();c++) {
            v_n(nodeIDs(e, c)) += v_e(e) / nodeIDs.cols();
        }
    }

    //setup A matrices
    A = KinematicTensorArray(nodeIDs.rows(), job->JOB_TYPE);
    Ainv = KinematicTensorArray(nodeIDs.rows(), job->JOB_TYPE);

    double s;
    for (int e=0;e<nodeIDs.rows();e++){
        //A
        A(e,0,0) = x_n(nodeIDs(e,1),0) - x_n(nodeIDs(e,0),0);
        A(e,1,0) = x_n(nodeIDs(e,1),1) - x_n(nodeIDs(e,0),1);
        A(e,2,0) = x_n(nodeIDs(e,1),2) - x_n(nodeIDs(e,0),2);

        A(e,0,1) = x_n(nodeIDs(e,2),0) - x_n(nodeIDs(e,0),0);
        A(e,1,1) = x_n(nodeIDs(e,2),1) - x_n(nodeIDs(e,0),1);
        A(e,2,1) = x_n(nodeIDs(e,2),2) - x_n(nodeIDs(e,0),2);

        A(e,0,2) = x_n(nodeIDs(e,3),0) - x_n(nodeIDs(e,0),0);
        A(e,1,2) = x_n(nodeIDs(e,3),1) - x_n(nodeIDs(e,0),1);
        A(e,2,2) = x_n(nodeIDs(e,3),2) - x_n(nodeIDs(e,0),2);

        //Ainv
        Ainv(e) = A(e).inverse();
    }

    //create node to element map, or search grid
    int tmpn;
    bool new_elem;
    if (use_n_to_e){
        //setup node to element mapping
        node_to_element_map.resize(node_count);
        for (int e=0; e<element_count; e++){
            for (int n=0; n<npe; n++){
                tmpn = nodeIDs(e,n);
                node_to_element_map[tmpn].push_back(e);
            }
            std::cout << "Initializing node_to_element_map... " << e << "/" << element_count << "\r";
        }
        std::cout << std::endl;
        //set up element neighbor mapping
        element_neighbors.resize(element_count);
        //loop through elements
        for (int e=0; e<element_count; e++){
            //loop over nodes in element
            for (int n=0; n<npe; n++){
                tmpn = nodeIDs(e,n);
                //loop over element list of node
                for (int i=0; i< node_to_element_map[tmpn].size(); i++){
                    new_elem = true;
                    //check if element is new
                    if (node_to_element_map[tmpn][i] == e){
                        //do not add self to list
                        continue;
                    }
                    for (int ee = 0; ee < element_neighbors[e].size(); ee++){
                        if (node_to_element_map[tmpn][i] == element_neighbors[e][ee]){
                            new_elem = false;
                            break;
                        }
                    }
                    if (new_elem){
                        //if element is not in unique list, add it
                        element_neighbors[e].push_back(node_to_element_map[tmpn][i]);
                    }
                }
            }
        }
    }

    //setup search grid
    //store length, number of linear nodes, and deltas
    Lx = KinematicVector(job->JOB_TYPE);
    Lxmin = KinematicVector(job->JOB_TYPE);
    Lxmax = KinematicVector(job->JOB_TYPE);
    Lxmin.setZero();
    Lxmax.setZero();
    for (int i = 0; i < x_n.size(); i++) {
        for (int pos = 0; pos < GRID_DIM; pos++) {
            if (x_n(i, pos) > Lxmax(pos)) {
                Lxmax(pos) = x_n(i, pos);
            }else if (x_n(i, pos) < Lxmin(pos)) {
                Lxmin(pos) = x_n(i, pos);
            }
        }
    }
    Lx = Lxmax - Lxmin;

    if (lc != -1.) {
        Nx = Eigen::VectorXi(job->DIM);;
        for (size_t pos = 0; pos < GRID_DIM; pos++) {
            Nx(pos) = std::floor(Lx(pos) / lc);
        }

        hx = KinematicVector(job->JOB_TYPE);
        for (int pos = 0; pos < GRID_DIM; pos++) {
            hx(pos) = Lx(pos) / Nx(pos);
        }
    } else {
        //hx already set by inputs
        Nx = Eigen::VectorXi(job->DIM);;
        for (size_t pos = 0; pos < GRID_DIM; pos++) {
            Nx(pos) = std::floor(Lx(pos) / hx(pos));
        }
    }

    for (int pos = GRID_DIM; pos < hx.size(); pos++) {
        Lx(pos) = 0;
        Nx(pos) = 0;
        hx(pos) = 0;
    }

    int len = 1; //number of cells in search grid
    for (int pos = 0; pos < GRID_DIM; pos++) {
        len *= Nx(pos);
    }

    //check if this is a large domain
    if (len > LARGE_DOMAIN_CUTOFF) {
        large_domain = true;
    } else {
        large_domain = false;
    }


    Eigen::VectorXi ijk = Eigen::VectorXi(GRID_DIM);
    int sc;

    //form min/max list
    element_min_max.resize(6 * element_count);
    for (int e = 0; e < element_count; e++) {
        /*
        //find search cell
        sc = whichSearchCell(x_n(nodeIDs(e,0)));
        //convert to ijk coords
        ijk = n_to_ijk(job, sc);
        */

        //get ijk from components of x_n
        ijk[0] = floor((x_n(nodeIDs(e, 0))[0] - Lxmin[0]) / hx[0]);
        ijk[1] = floor((x_n(nodeIDs(e, 0))[1] - Lxmin[1]) / hx[1]);
        ijk[2] = floor((x_n(nodeIDs(e, 0))[2] - Lxmin[2]) / hx[2]);

        //assign min/max from first node
        for (int pos = 0; pos < ijk.rows(); pos++) {
            element_min_max[6 * e + pos] = ijk[pos];
            element_min_max[6 * e + pos + 3] = ijk[pos];
        }

        //repeat for other 3 nodes
        for (int i = 1; i < npe; i++) {
            /*
            //find search cell
            sc = whichSearchCell(x_n(nodeIDs(e,i)));
            //convert to ijk coords
            ijk = n_to_ijk(job, sc);
             */
            ijk[0] = floor((x_n(nodeIDs(e, i))[0] - Lxmin[0]) / hx[0]);
            ijk[1] = floor((x_n(nodeIDs(e, i))[1] - Lxmin[1]) / hx[1]);
            ijk[2] = floor((x_n(nodeIDs(e, i))[2] - Lxmin[2]) / hx[2]);

            for (int pos = 0; pos < ijk.rows(); pos++) {
                //if ijk < min, assign to min
                if (ijk(pos) < element_min_max[6 * e + pos]) {
                    element_min_max[6 * e + pos] = ijk[pos];
                } else if (ijk(pos) > element_min_max[6 * e + 3 + pos]) {
                    element_min_max[6 * e + 3 + pos] = ijk[pos];
                }
            }
        }
    }

    std::cout << "Lx: " << EIGEN_MAP_OF_KINEMATIC_VECTOR(Lx).transpose() << std::endl;
    std::cout << "Nx: " << Nx.transpose() << std::endl;
    std::cout << "len: " << len << std::endl;

    //if small domain, use same method as in 2D
    if (!large_domain) {
        bool cell_in_range;
        //insert elements into search cell list
        for (int i = 0; i < len; i++) {
            search_offsets.push_back(search_cells.size());
            //find ijk of search cell
            ijk = n_to_ijk(job, i);
            //fill cell
            for (int e = 0; e < element_count; e++) {
                cell_in_range = true;
                for (int pos = 0; pos < ijk.rows(); pos++) {
                    //insert element if search cell is within min/max range
                    if (ijk(pos) < element_min_max[6 * e + pos] || ijk(pos) > element_min_max[6 * e + 3 + pos]) {
                        cell_in_range = false;
                        break;
                    }
                }
                //hasn't failed, therefore add to list
                if (cell_in_range) {
                    search_cells.push_back(e);
                }
            }
            std::cout << "Initializing... " << i << "/" << len << "\r";
        }
        search_offsets.push_back(search_cells.size());
    } else {
        //use large domain version
        //reserve std::vector
        search_cells_large_domain.resize(len);

        int n = -1;

        //loop over elements
        for (int e = 0; e < element_count; e++) {
            //loop over 3D ijk and insert element into search cell bins
            for (int i = element_min_max[6 * e + 0]; i <= element_min_max[6 * e + 3]; i++) {
                ijk(0) = i;
                for (int j = element_min_max[6 * e + 1]; j <= element_min_max[6 * e + 4]; j++) {
                    ijk(1) = j;
                    for (int k = element_min_max[6 * e + 2]; k <= element_min_max[6 * e + 5]; k++) {
                        ijk(2) = k;
                        //find associated cell
                        n = ijk_to_n(job, ijk);
                        //insert element number into cell bin
                        if (n < len) {
                            search_cells_large_domain[n].push_back(e);
                        } //else n is out of bounds...
                    }
                }
            }
            std::cout << "Initializing... " << e << "/" << element_count << "\r";
        }
    }

    return;
}

/*----------------------------------------------------------------------------*/
//
std::vector<int> TetrahedralGridLinear::getEdgeList(Job* job){
    std::vector<int> result = std::vector<int>(0);

    std::vector<std::vector<int>> node_neighbors = std::vector<std::vector<int>>(node_count);

    //loop over triangles and identify edges
    int n0, n1, n2, n3;
    bool already_counted = false;
    for (int e=0; e<element_count; e++) {
        //get node ids
        n0 = nodeIDs(e,0);
        n1 = nodeIDs(e,1);
        n2 = nodeIDs(e,2);
        n3 = nodeIDs(e,3);

        //for each edge, check if already accounted for in edge list
        //if not, add edge

        //EDGE 1 (n0 --- n1)
        already_counted = false;
        for (int i = 0; i<node_neighbors[n0].size(); i++){
            if (n1 == node_neighbors[n0][i]){
                already_counted = true;
                break;
            }
        }
        if (!already_counted){
            //add ids to neighbor lists
            node_neighbors[n0].push_back(n1);
            node_neighbors[n1].push_back(n0);

            //add edge to list
            result.push_back(n0);
            result.push_back(n1);
        }

        //EDGE 2 (n1 --- n2)
        already_counted = false;
        for (int i = 0; i<node_neighbors[n1].size(); i++){
            if (n2 == node_neighbors[n1][i]){
                already_counted = true;
                break;
            }
        }
        if (!already_counted){
            //add ids to neighbor lists
            node_neighbors[n1].push_back(n2);
            node_neighbors[n2].push_back(n1);

            //add edge to list
            result.push_back(n1);
            result.push_back(n2);
        }

        //EDGE 3 (n2 --- n0)
        already_counted = false;
        for (int i = 0; i<node_neighbors[n2].size(); i++){
            if (n0 == node_neighbors[n2][i]){
                already_counted = true;
                break;
            }
        }
        if (!already_counted){
            //add ids to neighbor lists
            node_neighbors[n2].push_back(n0);
            node_neighbors[n0].push_back(n2);

            //add edge to list
            result.push_back(n2);
            result.push_back(n0);
        }

        //EDGE 4 (n0 --- n3)
        already_counted = false;
        for (int i = 0; i<node_neighbors[n0].size(); i++){
            if (n3 == node_neighbors[n0][i]){
                already_counted = true;
                break;
            }
        }
        if (!already_counted){
            //add ids to neighbor lists
            node_neighbors[n3].push_back(n0);
            node_neighbors[n0].push_back(n3);

            //add edge to list
            result.push_back(n0);
            result.push_back(n3);
        }

        //EDGE 5 (n1 --- n3)
        already_counted = false;
        for (int i = 0; i<node_neighbors[n1].size(); i++){
            if (n3 == node_neighbors[n1][i]){
                already_counted = true;
                break;
            }
        }
        if (!already_counted){
            //add ids to neighbor lists
            node_neighbors[n3].push_back(n1);
            node_neighbors[n1].push_back(n3);

            //add edge to list
            result.push_back(n1);
            result.push_back(n3);
        }

        //EDGE 6 (n2 --- n3)
        already_counted = false;
        for (int i = 0; i<node_neighbors[n2].size(); i++){
            if (n3 == node_neighbors[n2][i]){
                already_counted = true;
                break;
            }
        }
        if (!already_counted){
            //add ids to neighbor lists
            node_neighbors[n3].push_back(n2);
            node_neighbors[n2].push_back(n3);

            //add edge to list
            result.push_back(n2);
            result.push_back(n3);
        }
    }

    return result;
}


/*----------------------------------------------------------------------------*/
//
void TetrahedralGridLinear::writeHeader(Job* job, Body* body, Serializer* serializer, std::ofstream& nfile, int SPEC){
    if (SPEC != Serializer::VTK){
        std::cerr << "ERROR: Unknown file SPEC in writeHeader: " << SPEC  << "! Exiting." << std::endl;
        exit(0);
    }

    int nlen = body->nodes->x.size();

    nfile << "ASCII\n";
    nfile << "DATASET UNSTRUCTURED_GRID\n";

    nfile << "POINTS " << nlen << " double\n";
    for (int i=0;i<nlen;i++){
        //vtk files require x,y,z
        for (int pos = 0; pos < 3; pos++){
            if (pos < GRID_DIM && (body->nodes->active(i) != 0) && std::isfinite(body->nodes->x(i,pos))){
                nfile << body->nodes->x(i,pos) << " ";
            } else {
                nfile << "0 ";
            }
        }
        nfile << "\n";
    }

    nfile << "POINTS " << node_count << " double\n";
    for (int i=0;i<node_count;i++){
        //vtk requires 3D position
        nfile << x_n(i,0) << " " << x_n(i,1) << " " << x_n(i,2) << "\n";
    }
    nfile << "CELLS " << element_count << " " << 5*element_count << "\n";
    for (int e=0;e<element_count;e++){
        nfile << "4 " << nodeIDs(e,0) << " " << nodeIDs(e,1) << " " << nodeIDs(e,2) << " " << nodeIDs(e,3) << "\n";
    }

    nfile << "CELL_TYPES " << element_count << "\n";
    for (int e=0;e<element_count;e++){
        nfile << "10\n";
    }
    nfile << "POINT_DATA " << nlen << "\n";
} //write cell types


/*----------------------------------------------------------------------------*/
//
int TetrahedralGridLinear::whichElement(Job* job, KinematicVector& xIN) {
    assert(xIN.VECTOR_TYPE == Lx.VECTOR_TYPE && "whichElement failed");
    //find search cell
    int searchID = whichSearchCell(xIN);
    if (searchID < 0) {
        return -1;
    }
    //search in that cell
    int elementID;
    if (!large_domain){
        for (int cellID = search_offsets[searchID]; cellID < search_offsets[searchID + 1]; cellID++) {
            elementID = search_cells[cellID];
            if (inElement(job, xIN, elementID)) {
                return elementID;
            }
        }
    } else {
        for(int i = 0; i < search_cells_large_domain[searchID].size(); i++){
            elementID = search_cells_large_domain[searchID][i];
            if (inElement(job, xIN, elementID)){
                return elementID;
            }
        }
    }

    /*
    std::cout << "OUT OF DOMAIN!!! " << std::endl;

    if (large_domain){
        std::cout << searchID << " :";
        for (int i=0; i < search_cells_large_domain[searchID].size(); i++){
            std::cout << " " << search_cells_large_domain[searchID][i];
        }
        std::cout << std::endl;
    }
    */

    return -1;
}

int TetrahedralGridLinear::whichElement(Job* job, KinematicVector& xIN, int& elem_guess) {
    assert(xIN.VECTOR_TYPE == Lx.VECTOR_TYPE && "whichElement failed");
    if (elem_guess > 0) {
        if (inElement(job, xIN, elem_guess)) {
            return elem_guess;
        }
        if (use_n_to_e) {
            int tmpn;
            //search nearest neighbors first
            for (int e=0; e<element_neighbors[elem_guess].size(); e++){
                if (inElement(job, xIN, element_neighbors[elem_guess][e])){
                    elem_guess = element_neighbors[elem_guess][e];
                    return elem_guess;
                }
            }
        }
    }

    //find search cell
    int searchID = whichSearchCell(xIN);
    if (searchID < 0) {
        return -1;
    }
    //search in that cell
    int elementID;
    if (!large_domain){
        for (int cellID = search_offsets[searchID]; cellID < search_offsets[searchID + 1]; cellID++) {
            elementID = search_cells[cellID];
            if (inElement(job, xIN, elementID)) {
                /*
                std::cout << "Element: " << elementID << " , Guess: " << elem_guess << std::endl << "Neighbors:";
                if (elem_guess > 0) {
                    for (int ee = 0; ee < element_neighbors[elem_guess].size(); ee++) {
                        std::cout << " " << element_neighbors[elem_guess][ee];
                    }
                }
                std::cout << std::endl;
                 */
                elem_guess = elementID;
                return elementID;
            }
        }
    } else {
        for (int i = 0; i < search_cells_large_domain[searchID].size(); i++){
            elementID = search_cells_large_domain[searchID][i];
            if (inElement(job, xIN, elementID)){
                elem_guess = elementID;
                return elementID;
            }
        }
    }

    //std::cout << "OUT OF DOMAIN!" << std::endl;

    elem_guess = -1;
    return -1;
}

bool TetrahedralGridLinear::inDomain(Job* job, KinematicVector& xIN){
    assert(xIN.VECTOR_TYPE == Lx.VECTOR_TYPE && "inDomain failed");

    if (whichSearchCell(xIN) == -1 || whichElement(job,xIN) < 0){
        return false;
    } else {
        return true;
    }
}

KinematicVector TetrahedralGridLinear::nodeIDToPosition(Job* job, int idIN){
    return x_n(idIN);
}

/*----------------------------------------------------------------------------*/
//
void TetrahedralGridLinear::evaluateBasisFnValue(Job* job, KinematicVector& xIN, std::vector<int>& nID, std::vector<double>& nVAL){
    int elementID = whichElement(job,xIN);
    if (elementID < 0){
        return;
    }

    KinematicVector xi = KinematicVector(job->JOB_TYPE);

    //xi = Ainv * (x - x0)
    xi = Ainv(elementID)*(xIN - x_n(nodeIDs(elementID,0)));

    //fill shape function vals
    double tmp;
    tmp = 1 - xi(0) - xi(1) - xi(2);
    nID.push_back(nodeIDs(elementID,0));
    nVAL.push_back(tmp);

    tmp = xi(0);
    nID.push_back(nodeIDs(elementID,1));
    nVAL.push_back(tmp);

    tmp = xi(1);
    nID.push_back(nodeIDs(elementID,2));
    nVAL.push_back(tmp);

    tmp = xi(2);
    nID.push_back(nodeIDs(elementID,3));
    nVAL.push_back(tmp);

    return;
}

void TetrahedralGridLinear::evaluateBasisFnGradient(Job* job, KinematicVector& xIN, std::vector<int>& nID, KinematicVectorArray& nGRAD){
    int elementID = whichElement(job,xIN);
    if (elementID < 0){
        return;
    }

    KinematicVector xi = KinematicVector(job->JOB_TYPE);
    KinematicVector tmpVec = KinematicVector(job->JOB_TYPE);
    KinematicTensor map = A(elementID);
    KinematicTensor inv_map = Ainv(elementID);

    //xi = Ainv * (x - x0)
    xi = inv_map*(xIN - x_n(nodeIDs(elementID,0)));

    //fill shape function vals
    tmpVec(0) = -1; tmpVec(1) = -1; tmpVec(2) = -1;
    tmpVec = inv_map.transpose()*tmpVec; //map*tmpVec;
    nID.push_back(nodeIDs(elementID,0));
    nGRAD.push_back(tmpVec);

    tmpVec(0) = 1; tmpVec(1) = 0; tmpVec(2) = 0;
    tmpVec = inv_map.transpose()*tmpVec; //map*tmpVec;
    nID.push_back(nodeIDs(elementID,1));
    nGRAD.push_back(tmpVec);

    tmpVec(0) = 0; tmpVec(1) = 1; tmpVec(2) = 0;
    tmpVec = inv_map.transpose()*tmpVec; //map*tmpVec;
    nID.push_back(nodeIDs(elementID,2));
    nGRAD.push_back(tmpVec);

    tmpVec(0) = 0; tmpVec(1) = 0; tmpVec(2) = 1;
    tmpVec = inv_map.transpose()*tmpVec; //map*tmpVec;
    nID.push_back(nodeIDs(elementID,3));
    nGRAD.push_back(tmpVec);

    return;
}

void TetrahedralGridLinear::evaluateBasisFnValue(Job* job,
                                                 KinematicVector& xIN,
                                                 std::vector<int>& nID,
                                                 std::vector<double>& nVAL,
                                                 int& elem_guess){
    int elementID = whichElement(job,xIN,elem_guess);
    if (elementID < 0){
        return;
    }

    KinematicVector xi = KinematicVector(job->JOB_TYPE);

    //xi = Ainv * (x - x0)
    xi = Ainv(elementID)*(xIN - x_n(nodeIDs(elementID,0)));

    //fill shape function vals
    double tmp;
    tmp = 1 - xi(0) - xi(1) - xi(2);
    nID.push_back(nodeIDs(elementID,0));
    nVAL.push_back(tmp);

    tmp = xi(0);
    nID.push_back(nodeIDs(elementID,1));
    nVAL.push_back(tmp);

    tmp = xi(1);
    nID.push_back(nodeIDs(elementID,2));
    nVAL.push_back(tmp);

    tmp = xi(2);
    nID.push_back(nodeIDs(elementID,3));
    nVAL.push_back(tmp);

    return;
}

void TetrahedralGridLinear::evaluateBasisFnGradient(Job* job,
                                                    KinematicVector& xIN,
                                                    std::vector<int>& nID,
                                                    KinematicVectorArray& nGRAD,
                                                    int& elem_guess){
    int elementID = whichElement(job,xIN,elem_guess);
    if (elementID < 0){
        return;
    }

    KinematicVector xi = KinematicVector(job->JOB_TYPE);
    KinematicVector tmpVec = KinematicVector(job->JOB_TYPE);
    KinematicTensor map = A(elementID);
    KinematicTensor inv_map = Ainv(elementID);

    //xi = Ainv * (x - x0)
    xi = inv_map*(xIN - x_n(nodeIDs(elementID,0)));

    //fill shape function vals
    tmpVec(0) = -1; tmpVec(1) = -1; tmpVec(2) = -1;
    tmpVec = inv_map.transpose()*tmpVec; //map*tmpVec;
    nID.push_back(nodeIDs(elementID,0));
    nGRAD.push_back(tmpVec);

    tmpVec(0) = 1; tmpVec(1) = 0; tmpVec(2) = 0;
    tmpVec = inv_map.transpose()*tmpVec; //map*tmpVec;
    nID.push_back(nodeIDs(elementID,1));
    nGRAD.push_back(tmpVec);

    tmpVec(0) = 0; tmpVec(1) = 1; tmpVec(2) = 0;
    tmpVec = inv_map.transpose()*tmpVec; //map*tmpVec;
    nID.push_back(nodeIDs(elementID,2));
    nGRAD.push_back(tmpVec);

    tmpVec(0) = 0; tmpVec(1) = 0; tmpVec(2) = 1;
    tmpVec = inv_map.transpose()*tmpVec; //map*tmpVec;
    nID.push_back(nodeIDs(elementID,3));
    nGRAD.push_back(tmpVec);

    return;
}

double TetrahedralGridLinear::nodeVolume(Job* job, int idIN){
    return v_n(idIN);
}

double TetrahedralGridLinear::elementVolume(Job* job, int idIN){
    return v_e(idIN);
}

int TetrahedralGridLinear::nodeTag(Job* job, int idIN){
    return nodeTags(idIN);
}

double TetrahedralGridLinear::nodeSurfaceArea(Job *job, int idIN) {
    //need to implement this at some point
    std::cout << "WARNING: TetrahedralGridLinear::nodeSurfaceArea() isn't implemented!" << std::endl;
    return 0;
}
