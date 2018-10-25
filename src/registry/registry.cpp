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
    object["MinimalVTK"] = &createInstance<Serializer,MinimalVTK>;
}


/*----------------------------------------------------------------------------*/
//drivers
template<>
Registry<Driver>::Registry() {
    object["DefaultDriver"] = &createInstance<Driver,DefaultDriver>;
    object["ColumnCollapseDriver"] = &createInstance<Driver,ColumnCollapseDriver>;
    object["UserDefinedGravityDriver"] = &createInstance<Driver,UserDefinedGravityDriver>;
    object["CavityFlowDriver"] = &createInstance<Driver,CavityFlowDriver>;
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
    object["WheelBody"] = &createInstance<Body,WheelBody>;
    object["HydrostaticBody"] = &createInstance<Body,HydrostaticBody>;
    object["LaunchedBody"] = &createInstance<Body,LaunchedBody>;
}


/*----------------------------------------------------------------------------*/
//contacts
template<>
Registry<Contact>::Registry() {
    object["ContactHuang"] = &createInstance<Contact,ContactHuang>;
    object["SlurryMixture"] = &createInstance<Contact,SlurryMixture>;
    object["SlurryContact"] = &createInstance<Contact,SlurryContact>;
    object["SlurryContact_ReflectedBoundary"] = &createInstance<Contact,SlurryContact_ReflectedBoundary>;
}


/*----------------------------------------------------------------------------*/
//grids
template<>
Registry<Grid>::Registry() {
    object["CartesianLinear"] = &createInstance<Grid,CartesianLinear>;
    object["CartesianCubic"] = &createInstance<Grid,CartesianCubic>;
    object["CartesianPeriodic"] = &createInstance<Grid,CartesianPeriodic>;
    object["CartesianCustom"] = &createInstance<Grid,CartesianCustom>;
    object["CartesianCubicCustom"] = &createInstance<Grid,CartesianCubicCustom>;
    object["TriangularGridLinear"] = &createInstance<Grid,TriangularGridLinear>;
    object["Regular2DTaylorCouetteCell"] = &createInstance<Grid,Regular2DTaylorCouetteCell>;
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
    object["Sand_SachithLocal"] = &createInstance<Material,Sand_SachithLocal>;
    object["SlurryGranularPhase"] = &createInstance<Material,SlurryGranularPhase>;
    object["SlurryFluidPhase"] = &createInstance<Material,SlurryFluidPhase>;
    object["BarotropicViscousFluid"] = &createInstance<Material,BarotropicViscousFluid>;
    object["Cornstarch"] = &createInstance<Material,Cornstarch>;
    object["CustomDSTModel"] = &createInstance<Material,CustomDSTModel>;
}


/*----------------------------------------------------------------------------*/
//boundaries
template<>
Registry<Boundary>::Registry() {
    object["CartesianBox"] = &createInstance<Boundary,CartesianBox>;
    object["CartesianSmoothBox"] = &createInstance<Boundary,CartesianSmoothBox>;
    object["CartesianFrictionalBox"] = &createInstance<Boundary,CartesianFrictionalBox>;
    object["CartesianBoxCustom"] = &createInstance<Boundary,CartesianBoxCustom>;
    object["Regular2DTaylorCouetteCustom"] = &createInstance<Boundary,Regular2DTaylorCouetteCustom>;
}