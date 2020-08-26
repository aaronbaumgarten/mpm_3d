//
// Created by aaron on 4/22/18.
// main.cpp
// for generation of 2D, square, overlapping gmsh meshes with variable discretization length
//

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <array>
#include <random>

//forward declare objects
class Node2D; class Face2D; class Element2D;

//class that stored grid information
class Grid2D{
public:
    Grid2D(){
        //do nothing
    }

    //grid discretization number
    int Nx = 1;

    //functions to add new components to grid and return pointer
    int getNewNode2D();
    int getNewFace2D(int a, int b, bool bounding_face);
    int getNewElement2D(int a, int b, int c);

    //function to refine grid
    void refine();

    //lists of components
    std::vector<Node2D> nodes;
    std::vector<Face2D> faces;
    std::vector<Element2D> elements;
};

//class that stores position of node
class Node2D{
public:
    Node2D(){
        //do nothing
    }

    //position
    std::array<double,2> x = {0, 0};

    //output id
    int id = -1;
};

//class that stores face information
class Face2D{
public:
    Face2D(Grid2D* grid, int a, int b, bool bounding_face){
        //pointer to parent
        master_grid = grid;

        //face defined by two end nodes
        end_points = {a,b};

        //set flag for boundary
        on_boundary = bounding_face;
    }

    //flag for whether face has been subdivided
    bool has_children = false;

    //output id
    int id = -1;

    //flag if bounding face
    bool on_boundary = false;

    //get functions
    const std::array<int,2>& getEndPoints(){
        return end_points;
    };

    int getMidPoint(){
        if (!has_children){
            //subdivide
            refine();
        }

        return mid_point;
    }

    const std::array<int,2>& getSubSegments(){
        if (!has_children){
            //subdivide
            refine();
        }

        return sub_segments;
    };

    //function to subdivide face
    void refine(){
        //don't refine if already refined
        if (!has_children){
            //set flag to true
            has_children = true;

            //ask grid for pointer to new Node2D object
            mid_point = master_grid->getNewNode2D();

            //set node coordinates
            double d = 0.5 + 0.1*(2.0*std::rand()/RAND_MAX - 1.0); //midpoint + shift?
            //double d = 0.5;
            Node2D* A = &master_grid->nodes[end_points[0]];
            Node2D* B = &master_grid->nodes[end_points[1]];
            Node2D* C = &master_grid->nodes[mid_point];

            C->x[0] = (1.0 - d) * A->x[0] + d * B->x[0];
            C->x[1] = (1.0 - d) * A->x[1] + d * B->x[1];

            //ask grid for index of new Face2D objects
            sub_segments[0] = master_grid->getNewFace2D(end_points[0], mid_point, on_boundary);
            sub_segments[1] = master_grid->getNewFace2D(mid_point, end_points[0], on_boundary);
        }
        return;
    }

private:
    //nodes that define face
    std::array<int,2> end_points = {-1,-1};

    //mid point, only accessible through get f'n
    int mid_point = -1;

    //subdivided faces
    std::array<int,2> sub_segments = {-1,-1};

    //grid object
    Grid2D* master_grid = nullptr;
};

//class that stores face information
class Element2D{
public:
    Element2D(Grid2D* grid, int a, int b, int c){
        //pointer to parent
        master_grid = grid;

        //face defined by two end nodes
        faces = {a, b, c};
    }

    //flag for whether element has been subdivided
    bool has_children = false;

    //output id
    int id = -1;

    //get functions
    const std::array<int,3>& getFaces(){
        return faces;
    };

    const std::array<int,4>& getSubElements(){
        if (!has_children){
            //subdivide
            refine();
        }

        return sub_elements;
    };

    //function to subdivide face
    void refine(){
        //don't refine if already refined
        if (!has_children){
            //set flag to true
            has_children = true;

            /*
             * A . . D . . B
             *  .         .
             *   .       .
             *    F     E
             *     .   .
             *      . .
             *       C
             */

            int A = master_grid->faces[faces[0]].getEndPoints()[0];
            int B = master_grid->faces[faces[0]].getEndPoints()[1];
            int D = master_grid->faces[faces[0]].getMidPoint();
            int AD = master_grid->faces[faces[0]].getSubSegments()[0];
            int DB = master_grid->faces[faces[0]].getSubSegments()[1];

            int C = -1;
            int E = -1;
            int F = -1;
            int BE = -1;
            int EC = -1;
            int AF = -1;
            int FC = -1;

            for (int i=1; i<3; i++) {
                if (master_grid->faces[faces[i]].getEndPoints()[0] == A) {
                    //this is face AC
                    C = master_grid->faces[faces[i]].getEndPoints()[1];
                    F = master_grid->faces[faces[i]].getMidPoint();
                    AF = master_grid->faces[faces[i]].getSubSegments()[0];
                    FC = master_grid->faces[faces[i]].getSubSegments()[1];
                } else if (master_grid->faces[faces[i]].getEndPoints()[1] == A){
                    //this is face AC
                    C = master_grid->faces[faces[i]].getEndPoints()[0];
                    F = master_grid->faces[faces[i]].getMidPoint();
                    AF = master_grid->faces[faces[i]].getSubSegments()[1];
                    FC = master_grid->faces[faces[i]].getSubSegments()[0];
                } else if (master_grid->faces[faces[i]].getEndPoints()[0] == B) {
                    //this is face BC
                    C = master_grid->faces[faces[i]].getEndPoints()[1];
                    E = master_grid->faces[faces[i]].getMidPoint();
                    BE = master_grid->faces[faces[i]].getSubSegments()[0];
                    EC = master_grid->faces[faces[i]].getSubSegments()[1];
                } else if (master_grid->faces[faces[i]].getEndPoints()[1] == B){
                    //this is face AC
                    C = master_grid->faces[faces[i]].getEndPoints()[0];
                    E = master_grid->faces[faces[i]].getMidPoint();
                    BE = master_grid->faces[faces[i]].getSubSegments()[1];
                    EC = master_grid->faces[faces[i]].getSubSegments()[0];
                } else {
                    std::cerr << "ERROR! Element ABC missing segments AC and BC!" << std::endl;
                    exit(0);
                }
            }

            if (C == -1 || E == -1 || F == -1 || BE == -1 || EC == -1 || AF == -1 || FC == -1){
                std::cerr << "ERROR! Element ABC missing segments or nodes!" << std::endl;
                exit(0);
            }

            //create element BDE using un-oriented faces AB and BC
            int ED = master_grid->getNewFace2D(E, D, false);
            int BDE = master_grid->getNewElement2D(BE, ED, DB);

            int FE = master_grid->getNewFace2D(F, E, false);
            int CFE = master_grid->getNewElement2D(FC, FE, EC);

            int FD = master_grid->getNewFace2D(F, D, false);
            int AFD = master_grid->getNewElement2D(AF, FD, AD);

            int DEF = master_grid->getNewElement2D(ED, FD, FE);

            //assign sub_elements list
            sub_elements = {BDE, CFE, AFD, DEF};

        }
        return;
    }

private:
    //faces that define element
    std::array<int,3> faces = {-1, -1, -1};

    //subdivided faces
    std::array<int,4> sub_elements = {-1, -1, -1, -1};

    //pointer to grid object
    Grid2D* master_grid = nullptr;
};

/*----------------------------------------------------------------------------*/

//functions to add new components to grid and return pointer
int Grid2D::getNewNode2D(){
    //add new node to list
    nodes.push_back(Node2D());

    //return index of new node
    return (nodes.size() - 1);
}

int Grid2D::getNewFace2D(int a, int b, bool bounding_face){
    //add new face to list
    faces.push_back(Face2D(this, a, b, bounding_face));

    //return index of new face
    return (faces.size() - 1);
}

int Grid2D::getNewElement2D(int a, int b, int c){
    //add new element to list
    elements.push_back(Element2D(this, a, b, c));

    //return index of new element
    return (elements.size() - 1);
}

//function to refine grid
void Grid2D::refine(){
    //to refine grid, ask all elements in current list to refine themselves.
    int current_element_count = elements.size();
    for (int e=0; e<current_element_count; e++){
        elements[e].refine();
    }

    Nx *= 2;

    return;
}

/*----------------------------------------------------------------------------*/

int main(int argc, char *argv[]) {
    std::cout << "Hello, World!" << std::endl;
    std::cout << argv[0] << std::endl;

    size_t random_seed = 0;
    std::srand(random_seed);
    std::cout << "Random Seed: " << random_seed << std::endl;

    std::cout << "Exiting." << std::endl;

    return 0;
}