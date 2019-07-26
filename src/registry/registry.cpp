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
    object["SlicePointsVTK"] = &createInstance<Serializer,SlicePointsVTK>;
}


/*----------------------------------------------------------------------------*/
//drivers
template<>
Registry<Driver>::Registry() {
    object["DefaultDriver"] = &createInstance<Driver,DefaultDriver>;
    object["ColumnCollapseDriver"] = &createInstance<Driver,ColumnCollapseDriver>;
    object["UserDefinedGravityDriver"] = &createInstance<Driver,UserDefinedGravityDriver>;
    object["CavityFlowDriver"] = &createInstance<Driver,CavityFlowDriver>;
    object["BallisticDriver"] = &createInstance<Driver,BallisticDriver>;
}


/*----------------------------------------------------------------------------*/
//solvers
template<>
Registry<Solver>::Registry() {
    object["ExplicitUSL"] = &createInstance<Solver,ExplicitUSL>;
    object["ParallelExplicitUSL"] = &createInstance<Solver,ParallelExplicitUSL>;
    object["ThreadPoolExplicitUSL"] = &createInstance<Solver,ThreadPoolExplicitUSL>;
}


/*----------------------------------------------------------------------------*/
//
template<>
Registry<Body>::Registry() {
    object["DefaultBody"] = &createInstance<Body,DefaultBody>;
    object["WheelBody"] = &createInstance<Body,WheelBody>;
    object["HydrostaticBody"] = &createInstance<Body,HydrostaticBody>;
    object["LaunchedBody"] = &createInstance<Body,LaunchedBody>;
    object["PrestressedBody"] = &createInstance<Body,PrestressedBody>;
}


/*----------------------------------------------------------------------------*/
//contacts
template<>
Registry<Contact>::Registry() {
    object["ContactHuang"] = &createInstance<Contact,ContactHuang>;
    object["SlurryMixture"] = &createInstance<Contact,SlurryMixture>;
    object["SlurryContact"] = &createInstance<Contact,SlurryContact>;
    object["SlurryContact_ReflectedBoundary"] = &createInstance<Contact,SlurryContact_ReflectedBoundary>;
    object["ContactHuang_ReflectedBoundary"] = &createInstance<Contact,ContactHuang_ReflectedBoundary>;
    object["ContactRigid_ReflectedBoundary"] = &createInstance<Contact,ContactRigid_ReflectedBoundary>;
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
    object["CartesianCubic_Offset"] = &createInstance<Grid,CartesianCubic_Offset>;
    object["TetrahedralGridLinear"] = &createInstance<Grid,TetrahedralGridLinear>;
}


/*----------------------------------------------------------------------------*/
//points
template<>
Registry<Points>::Registry() {
    object["DefaultPoints"] = &createInstance<Points,DefaultPoints>;
    object["GmshPoints"] = &createInstance<Points,GmshPoints>;
    object["CartesianPoints"] = &createInstance<Points,CartesianPoints>;
    object["ThreadPoolPoints"] = &createInstance<Points,ThreadPoolPoints>;
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
}


/*----------------------------------------------------------------------------*/
//boundaries
template<>
Registry<Boundary>::Registry() {
    object["CartesianBox"] = &createInstance<Boundary,CartesianBox>;
    object["CartesianSmoothBox"] = &createInstance<Boundary,CartesianSmoothBox>;
    object["CartesianFrictionalBox"] = &createInstance<Boundary,CartesianFrictionalBox>;
    object["CartesianBoxCustom"] = &createInstance<Boundary,CartesianBoxCustom>;
    object["AssignedVelocity"] = &createInstance<Boundary,AssignedVelocity>;
}

/*----------------------------------------------------------------------------*/
//parts
template<>
Registry<Part>::Registry() {
    object["Ball"] = &createInstance<Part,Ball>;
    object["Box"] = &createInstance<Part,Box>;
}