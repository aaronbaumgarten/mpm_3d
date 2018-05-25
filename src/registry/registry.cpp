//
// Created by aaron on 4/24/18.
// registry.cpp
//

#include "registry.hpp"
#include "mpm_objects.hpp"

#include "objects/serializers/serializers.hpp"
#include "objects/drivers/drivers.hpp"
#include "objects/solvers/solvers.hpp"

#include "objects/bodies/bodies.hpp"
#include "objects/contacts/contacts.hpp"
#include "objects/grids/grids.hpp"

#include "objects/points/points.hpp"
#include "objects/nodes/nodes.hpp"
#include "objects/materials/materials.hpp"
#include "objects/boundaries/boundaries.hpp"

template<typename Base, typename Derived> std::unique_ptr<Base> createInstance() { return std::unique_ptr<Base>(new Derived); }

/*template<>
Registry<Base>::Registry() {
    object["DerivedA"] = &createInstance<Base,DerivedA>;
    object["DerivedB"] = &createInstance<Base,DerivedB>;
}*/

/*----------------------------------------------------------------------------*/
//serializers
template<>
Registry<Serializer>::Registry() {
    object["DefaultVTK"] = &createInstance<Serializer,DefaultVTK>;
}


/*----------------------------------------------------------------------------*/
//drivers
template<>
Registry<Driver>::Registry() {
    object["DefaultDriver"] = &createInstance<Driver,DefaultDriver>;
    object["ColumnCollapseDriver"] = &createInstance<Driver,ColumnCollapseDriver>;
    object["UserDefinedGravityDriver"] = &createInstance<Driver,UserDefinedGravityDriver>;
}


/*----------------------------------------------------------------------------*/
//solvers
template<>
Registry<Solver>::Registry() {
    object["ExplicitUSL"] = &createInstance<Solver,ExplicitUSL>;
}


/*----------------------------------------------------------------------------*/
//
template<>
Registry<Body>::Registry() {
    object["DefaultBody"] = &createInstance<Body,DefaultBody>;
}


/*----------------------------------------------------------------------------*/
//contacts
template<>
Registry<Contact>::Registry() {
    object["ContactHuang"] = &createInstance<Contact,ContactHuang>;
    object["SlurryMixture"] = &createInstance<Contact,SlurryMixture>;
}


/*----------------------------------------------------------------------------*/
//grids
template<>
Registry<Grid>::Registry() {
    object["CartesianLinear"] = &createInstance<Grid,CartesianLinear>;
    object["CartesianCubic"] = &createInstance<Grid,CartesianCubic>;
    object["CartesianPeriodic"] = &createInstance<Grid,CartesianPeriodic>;
}


/*----------------------------------------------------------------------------*/
//points
template<>
Registry<Points>::Registry() {
    object["DefaultPoints"] = &createInstance<Points,DefaultPoints>;
}


/*----------------------------------------------------------------------------*/
//nodes
template<>
Registry<Nodes>::Registry() {
    object["DefaultNodes"] = &createInstance<Nodes,DefaultNodes>;
}


/*----------------------------------------------------------------------------*/
//materials
template<>
Registry<Material>::Registry() {
    object["IsotropicLinearElasticity"] = &createInstance<Material,IsotropicLinearElasticity>;
}


/*----------------------------------------------------------------------------*/
//boundaries
template<>
Registry<Boundary>::Registry() {
    object["CartesianBox"] = &createInstance<Boundary,CartesianBox>;
    object["CartesianSmoothBox"] = &createInstance<Boundary,CartesianSmoothBox>;
    object["CartesianFrictionalBox"] = &createInstance<Boundary,CartesianFrictionalBox>;
    object["CartesianBoxCustom"] = &createInstance<Boundary,CartesianBoxCustom>;
}