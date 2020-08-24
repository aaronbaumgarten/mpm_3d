//
// Created by aaron on 4/22/18.
// main.cpp
// for generation of 2D, square, overlapping gmsh meshes with variable discretization length
//

#include <iostream>
#include <fstream>
#include <stdlib.h>

//class that stores position of node
class Node2D{
public:
    Node2D(){
        //do nothing
    }

    //position
    std::array<double,2> x = {0, 0};

    //output id
    int id = 0;
};

//class that stores face information
class Face2D{
public:
    Face2D(){
        //do nothing
    }

    bool has_children = false;

private:
    //nodes that define face
    std::array<Node2D*,2> end_points = {Node2D(), Node2D()};

    //mid point, only accessible through get f'n
    Node2D mid_point = Node2D();

    //

};

int main(int argc, char *argv[]) {
    std::cout << "Hello, World!" << std::endl;
    std::cout << argv[0] << std::endl;

    std::cout << "Exiting." << std::endl;

    return 0;
}