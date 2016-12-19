#include <vector>
#include <Eigen/Sparse>
#include <iostream>
#include "unit.hpp"
#include "element.hpp"

struct Pt {
    double x[1];
    double y[1];
    double z[1];
};

struct MockBody
{
    Pt nodes[8];
    std::vector<Eigen::Triplet<double>> PhiTriplets;
    std::vector<Eigen::Triplet<double>> gradPhiXTriplets;
    std::vector<Eigen::Triplet<double>> gradPhiYTriplets;
    std::vector<Eigen::Triplet<double>> gradPhiZTriplets;
};

struct MockParticle
{
    double x[1];
    double y[1];
    double z[1];
    double corner[8][3];
    size_t id;
    double a[1];
    double v_averaging[1];
};

UTEST(shapefunction, "Testing phi and gradPhi", []() -> bool {
    bool passing = true;
    MockBody mb;
    // lay out the master element.
    mb.nodes[0].x[0] = -1;
    mb.nodes[1].x[0] = 1;
    mb.nodes[2].x[0] = -1;
    mb.nodes[3].x[0] = 1;
    mb.nodes[4].x[0] = -1;
    mb.nodes[5].x[0] = 1;
    mb.nodes[6].x[0] = -1;
    mb.nodes[7].x[0] = 1;

    mb.nodes[0].y[0] = -1;
    mb.nodes[1].y[0] = -1;
    mb.nodes[2].y[0] = 1;
    mb.nodes[3].y[0] = 1;
    mb.nodes[4].y[0] = -1;
    mb.nodes[5].y[0] = -1;
    mb.nodes[6].y[0] = 1;
    mb.nodes[7].y[0] = 1;

    mb.nodes[0].z[0] = -1;
    mb.nodes[1].z[0] = -1;
    mb.nodes[2].z[0] = -1;
    mb.nodes[3].z[0] = -1;
    mb.nodes[4].z[0] = 1;
    mb.nodes[5].z[0] = 1;
    mb.nodes[6].z[0] = 1;
    mb.nodes[7].z[0] = 1;

    MockParticle mp;

    mp.x[0] = 0;
    mp.y[0] = 0;
    mp.z[0] = 0;
    mp.corner[0][0] = 0;
    mp.corner[0][1] = 0;
    mp.corner[0][2] = 0;
    mp.a[0] = 0.5;
    mp.v_averaging[0] = 1;

    Element el;
    el.numNodes = 8;
    el.nodeID.push_back(0);
    el.nodeID.push_back(1);
    el.nodeID.push_back(2);
    el.nodeID.push_back(3);
    el.nodeID.push_back(4);
    el.nodeID.push_back(5);
    el.nodeID.push_back(6);
    el.nodeID.push_back(7);
    el.calculatePhic(&mb, &mp, 0, 0);

    passing &= (mb.PhiTriplets[0].value() == (1.0 / 8.0) / 8.0);
    passing &= (mb.PhiTriplets[1].value() == (1.0 / 8.0) / 8.0);
    passing &= (mb.PhiTriplets[2].value() == (1.0 / 8.0) / 8.0);
    passing &= (mb.PhiTriplets[3].value() == (1.0 / 8.0) / 8.0);
    passing &= (mb.PhiTriplets[4].value() == (1.0 / 8.0) / 8.0);
    passing &= (mb.PhiTriplets[5].value() == (1.0 / 8.0) / 8.0);
    passing &= (mb.PhiTriplets[6].value() == (1.0 / 8.0) / 8.0);
    passing &= (mb.PhiTriplets[7].value() == (1.0 / 8.0) / 8.0);
    // std::cout << mb.PhiTriplets[0].value() << '\n';

    /*
    // Still have to think about this...
    passing &= (mb.gradPhiZTriplets[0].value() == -(0.5 * 0.5 / 1.0) / 8.0);
    passing &= (mb.gradPhiZTriplets[1].value() == -(0.5 * 0.5 / 1.0) / 8.0);
    passing &= (mb.gradPhiZTriplets[2].value() == -(0.5 * 0.5 / 1.0) / 8.0);
    passing &= (mb.gradPhiZTriplets[3].value() == -(0.5 * 0.5 / 1.0) / 8.0);
    passing &= (mb.gradPhiZTriplets[4].value() == (0.5 * 0.5 / 1.0) / 8.0);
    passing &= (mb.gradPhiZTriplets[5].value() == (0.5 * 0.5 / 1.0) / 8.0);
    passing &= (mb.gradPhiZTriplets[6].value() == (0.5 * 0.5 / 1.0) / 8.0);
    passing &= (mb.gradPhiZTriplets[7].value() == (0.5 * 0.5 / 1.0) / 8.0);
    std::cout << mb.gradPhiZTriplets[0].value() << '\n';
    std::cout << mb.gradPhiZTriplets[1].value() << '\n';
    std::cout << mb.gradPhiZTriplets[2].value() << '\n';
    std::cout << mb.gradPhiZTriplets[3].value() << '\n';
    std::cout << mb.gradPhiZTriplets[4].value() << '\n';
    std::cout << mb.gradPhiZTriplets[5].value() << '\n';
    std::cout << mb.gradPhiZTriplets[6].value() << '\n';
    std::cout << mb.gradPhiZTriplets[7].value() << '\n';
    */
    return passing;
});
