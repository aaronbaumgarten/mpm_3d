//
// Created by aaron on 4/24/18.
// registry.cpp
//

#include "registry.hpp"
#include "mpm_objects.hpp"

#include "objects/serializers/serializers.hpp"
#include "objects/drivers/drivers.hpp"

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
}


/*----------------------------------------------------------------------------*/
//points
template<>
Registry<Points>::Registry() {
    object["DefaulPoints"] = &createInstance<Points,DefaultPoints>;
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
}