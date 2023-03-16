//
// Created by aaron on 4/24/18.
// registry.cpp
//

#include <fvm/fvm_grids.hpp>
#include <fvm/fvm_bodies.hpp>
#include "registry.hpp"
#include "mpm_objects.hpp"
#include "fvm_objects.hpp"

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

#include "fvm/fvm_grids.hpp"
#include "fvm/fvm_bodies.hpp"
#include "fvm/fvm_drivers.hpp"
#include "fvm/fvm_materials.hpp"
#include "fvm/fvm_solvers.hpp"
#include "fvm/fvm_serializers.hpp"

template<typename Base, typename Derived>
static std::unique_ptr<Base> createInstance() {
    return std::unique_ptr<Base>(new Derived);
}

// REGISTER allows the usage of DERIVED (with the name of DERIVED as a string)
// by placing it into the appropriate registry based on BASE. As a concrete
// example, calling
//
// REGISTER(Serializer, DefaultVTK)
//
// allows you to specify the string "DefaultVTK" as the Serializer class in the
// configuration file.
//
// Note that it is a macro, not a template, so we can use the stringization
// operator `#` to automatically generate the name from DERIVED.
#define REGISTER(BASE, DERIVED) do { \
    object[#DERIVED] = &createInstance<BASE, DERIVED>; \
} while (0)

/*template<>
Registry<Base>::Registry() {
    object["DerivedA"] = &createInstance<Base,DerivedA>;
    object["DerivedB"] = &createInstance<Base,DerivedB>;
}*/

/*----------------------------------------------------------------------------*/
//serializers
template<>
Registry<Serializer>::Registry() {
    REGISTER(Serializer, DefaultVTK);
    REGISTER(Serializer, MinimalVTK);
    REGISTER(Serializer, SlicePointsVTK);
    REGISTER(Serializer, DebugVTK);
}


/*----------------------------------------------------------------------------*/
//drivers
template<>
Registry<Driver>::Registry() {
    REGISTER(Driver, DefaultDriver);
    REGISTER(Driver, ColumnCollapseDriver);
    REGISTER(Driver, UserDefinedGravityDriver);
    REGISTER(Driver, CavityFlowDriver);
    REGISTER(Driver, BallisticDriver);
    REGISTER(Driver, FiniteVolumeDriver);
    REGISTER(Driver, FishRandomSampleDriver);
    REGISTER(Driver, FishMixtureModelRandomSampleDriver);
    REGISTER(Driver, FVMVariableStepDriver);
    REGISTER(Driver, FVMNumericalDampingDriver);
    REGISTER(Driver, ChuteFlowDriver);
    REGISTER(Driver, FVMBolideImpactDriver);
    REGISTER(Driver, FVMBolideRestartDriver);
}


/*----------------------------------------------------------------------------*/
//solvers
template<>
Registry<Solver>::Registry() {
    REGISTER(Solver, ExplicitUSL);
    REGISTER(Solver, ParallelExplicitUSL);
    REGISTER(Solver, ThreadPoolExplicitUSL);
    REGISTER(Solver, TGVErrorSolver);
    REGISTER(Solver, ShiftedTGVErrorSolver);
    REGISTER(Solver, ExplicitUSLwithVolumetricStrainSmoothing);
    REGISTER(Solver, GeneralizedVortexErrorSolver);
    REGISTER(Solver, ExplicitMixtureSolver);
}


/*----------------------------------------------------------------------------*/
//
template<>
Registry<Body>::Registry() {
    REGISTER(Body, DefaultBody);
    REGISTER(Body, WheelBody);
    REGISTER(Body, HydrostaticBody);
    REGISTER(Body, LaunchedBody);
    REGISTER(Body, PrestressedBody);
}


/*----------------------------------------------------------------------------*/
//contacts
template<>
Registry<Contact>::Registry() {
    REGISTER(Contact, ContactHuang);
    REGISTER(Contact, SlurryMixture);
    REGISTER(Contact, SlurryContact);
    REGISTER(Contact, SlurryContact_ReflectedBoundary);
    REGISTER(Contact, ContactHuang_ReflectedBoundary);
    REGISTER(Contact, ContactRigid_ReflectedBoundary);
    REGISTER(Contact, ContactLinked);
}


/*----------------------------------------------------------------------------*/
//grids
template<>
Registry<Grid>::Registry() {
    REGISTER(Grid, CartesianLinear);
    REGISTER(Grid, CartesianCubic);
    REGISTER(Grid, CartesianPeriodic);
    REGISTER(Grid, CartesianCustom);
    REGISTER(Grid, CartesianCubicCustom);
    REGISTER(Grid, TriangularGridLinear);
    REGISTER(Grid, CartesianCubic_Offset);
    REGISTER(Grid, TetrahedralGridLinear);
    REGISTER(Grid, Linear1DNonUniform);
    REGISTER(Grid, CartesianQuadratic);
    REGISTER(Grid, CartesianQuadraticCustom);
    REGISTER(Grid, CartesianUGIMP);
    REGISTER(Grid, CartesianUGIMPCustom);
}


/*----------------------------------------------------------------------------*/
//points
template<>
Registry<Points>::Registry() {
    REGISTER(Points, DefaultPoints);
    REGISTER(Points, GmshPoints);
    REGISTER(Points, CartesianPoints);
    REGISTER(Points, ThreadPoolPoints);
    REGISTER(Points, ImprovedQuadraturePoints);
}


/*----------------------------------------------------------------------------*/
//nodes
template<>
Registry<Nodes>::Registry() {
    REGISTER(Nodes, DefaultNodes);
}


/*----------------------------------------------------------------------------*/
//materials
template<>
Registry<Material>::Registry() {
    REGISTER(Material, IsotropicLinearElasticity);
    REGISTER(Material, Sand_SachithLocal);
    REGISTER(Material, SlurryGranularPhase);
    REGISTER(Material, SlurryFluidPhase);
    REGISTER(Material, BarotropicViscousFluid);
    REGISTER(Material, Cornstarch);
    REGISTER(Material, SlurryGranularPhase_wUnderCompaction);
    REGISTER(Material, Fish);
    REGISTER(Material, CompressibleNeohookeanElasticity);
    REGISTER(Material, BreakageMechanicsSand);
    REGISTER(Material, CompressibleBreakageMechanicsSand);
    REGISTER(Material, TillotsonEOSFluid);
}


/*----------------------------------------------------------------------------*/
//boundaries
template<>
Registry<Boundary>::Registry() {
    REGISTER(Boundary, CartesianBox);
    REGISTER(Boundary, CartesianSmoothBox);
    REGISTER(Boundary, CartesianFrictionalBox);
    REGISTER(Boundary, CartesianBoxCustom);
    REGISTER(Boundary, AssignedVelocity);
    REGISTER(Boundary, GeneralCustomBoundary);
}

/*----------------------------------------------------------------------------*/
//parts
template<>
Registry<Part>::Registry() {
    REGISTER(Part, Ball);
    REGISTER(Part, Box);
    REGISTER(Part, SineWave);
    REGISTER(Part, SandPile);
}

/*----------------------------------------------------------------------------*/
//finite volume objects
template<>
Registry<FiniteVolumeSolver>::Registry() {
    REGISTER(FiniteVolumeSolver, FVMDefaultSolver);
    REGISTER(FiniteVolumeSolver, FVMRungeKuttaSolver);
    REGISTER(FiniteVolumeSolver, FVMSteadyStateSolver);
    REGISTER(FiniteVolumeSolver, FVMMixtureSolver);
    REGISTER(FiniteVolumeSolver, FVMMixtureSolverRK4);
    REGISTER(FiniteVolumeSolver, FVMStaticMixtureSolverRK4);
    REGISTER(FiniteVolumeSolver, ParallelMixtureSolverRK4);
}

template<>
Registry<FiniteVolumeGrid>::Registry() {
    REGISTER(FiniteVolumeGrid, FVMCartesian);
    REGISTER(FiniteVolumeGrid, FVMGmsh2D);
    REGISTER(FiniteVolumeGrid, FVMLinear1DNonUniform);
    REGISTER(FiniteVolumeGrid, FVMGmsh3D);
}

template<>
Registry<FiniteVolumeBody>::Registry() {
    REGISTER(FiniteVolumeBody, FVMDefaultBody);
    REGISTER(FiniteVolumeBody, FVMRocketBody);
}

template<>
Registry<FiniteVolumeMaterial>::Registry() {
    REGISTER(FiniteVolumeMaterial, FVMBarotropicViscousFluid);
    REGISTER(FiniteVolumeMaterial, FVMSlurryFluidPhase);
    REGISTER(FiniteVolumeMaterial, FVMIdealGas);
    REGISTER(FiniteVolumeMaterial, FVMSlurryGasPhase);
    REGISTER(FiniteVolumeMaterial, FVMCarmanKozenyFluid);
}

template<>
Registry<FiniteVolumeSerializer>::Registry() {
    REGISTER(FiniteVolumeSerializer, FVMDefaultVTK);
}
