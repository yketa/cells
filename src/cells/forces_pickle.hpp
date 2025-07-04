/*
Provide pickling ability for C++ force computation objects via pybind11.
*/

#ifndef FORCES_PICKLE_HPP
#define FORCES_PICKLE_HPP

#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "base_pickle.hpp"
#include "forces.hpp"
#include "system.hpp"

/*
 *  PerimeterForce
 *
 */

// save state
pybind11::tuple PerimeterForce::pybind11_getstate() const {
    return pybind11::make_tuple(
        // unique identifying string for the force object
        "PerimeterForce",
        // state
        parameters);
}

// load state
template<> void
VertexModel::addVertexForce<PerimeterForce, pybind11::tuple const&>(
    std::string const& name, pybind11::tuple const& t) {
    // check
    checkSize(t, 2);
    assert(t[0].cast<std::string>() == "PerimeterForce");
    // initialise force
    ParametersType const parameters = t[1].cast<ParametersType>();
    addVertexForce<PerimeterForce, double const&, double const&>(
        name, parameters.at("kP"), parameters.at("P0"));
}

/*
 *  AreaForce
 *
 */

// save state
pybind11::tuple AreaForce::pybind11_getstate() const {
    return pybind11::make_tuple(
        // unique identifying string for the force object
        "AreaForce",
        // state
        parameters);
}

// load state
template<> void
VertexModel::addVertexForce<AreaForce, pybind11::tuple const&>(
    std::string const& name, pybind11::tuple const& t) {
    // check
    checkSize(t, 2);
    assert(t[0].cast<std::string>() == "AreaForce");
    // initialise force
    ParametersType const parameters = t[1].cast<ParametersType>();
    addVertexForce<AreaForce, double const&, double const&>(
        name, parameters.at("kA"), parameters.at("A0"));
}

/*
 *  SurfaceForce
 *
 */

// save state
pybind11::tuple SurfaceForce::pybind11_getstate() const {
    return pybind11::make_tuple(
        // unique identifying string for the force object
        "SurfaceForce",
        // state
        parameters, volume);
}

// load state
template<> void
VertexModel::addVertexForce<SurfaceForce, pybind11::tuple const&>(
    std::string const& name, pybind11::tuple const& t) {
    // check
    checkSize(t, 3);
    assert(t[0].cast<std::string>() == "SurfaceForce");
    // initialise force
    ParametersType const parameters = t[1].cast<ParametersType>();
    addVertexForce<SurfaceForce,
        double const&, double const&, double const&>(
        name,
        parameters.at("Lambda"), parameters.at("V0"), parameters.at("tauV"));
    // set internal degrees of freedom state
    std::map<long int, double> const volume =
        t[2].cast<std::map<long int, double>>();
    std::shared_ptr<SurfaceForce> sf =
        std::static_pointer_cast<SurfaceForce>(
            vertexForces[name]);
    sf->setVolume(volume);
}

/*
 *  VolumeForce
 *
 */

// save state
pybind11::tuple VolumeForce::pybind11_getstate() const {
    return pybind11::make_tuple(
        // unique identifying string for the force object
        "VolumeForce",
        // state
        parameters, height);
}

// load state
template<> void
VertexModel::addVertexForce<VolumeForce, pybind11::tuple const&>(
    std::string const& name, pybind11::tuple const& t) {
    // check
    checkSize(t, 3);
    assert(t[0].cast<std::string>() == "VolumeForce");
    // initialise force
    ParametersType const parameters = t[1].cast<ParametersType>();
    addVertexForce<VolumeForce, double const&, double const&, double const&>(
        name, parameters.at("kV"), parameters.at("H0"), parameters.at("A0"));
    // set internal degrees of freedom state
    std::map<long int, double> const height =
        t[2].cast<std::map<long int, double>>();
    std::shared_ptr<VolumeForce> vf =
        std::static_pointer_cast<VolumeForce>(
            vertexForces[name]);
    vf->setHeight(height);
}

/*
 *  LinearVolumeForce
 *
 */

// save state
pybind11::tuple LinearVolumeForce::pybind11_getstate() const {
    return pybind11::make_tuple(
        // unique identifying string for the force object
        "LinearVolumeForce",
        // state
        parameters, height);
}

// load state
template<> void
VertexModel::addVertexForce<LinearVolumeForce, pybind11::tuple const&>(
    std::string const& name, pybind11::tuple const& t) {
    // check
    checkSize(t, 3);
    assert(t[0].cast<std::string>() == "LinearVolumeForce");
    // initialise force
    ParametersType const parameters = t[1].cast<ParametersType>();
    addVertexForce<LinearVolumeForce,
        double const&, double const&,
        double const&, double const&,
        double const&, double const&, double const&>(
        name,
        parameters.at("kA"), parameters.at("A0"),
        parameters.at("kP"), parameters.at("P0"),
        parameters.at("taur"), parameters.at("H0"), parameters.at("taua"));
    // set internal degrees of freedom state
    std::map<long int, double> const height =
        t[2].cast<std::map<long int, double>>();
    std::shared_ptr<LinearVolumeForce> vf =
        std::static_pointer_cast<LinearVolumeForce>(
            vertexForces[name]);
    vf->setHeight(height);
}

/*
 *  PressureForce
 *
 */

// save state
pybind11::tuple PressureForce::pybind11_getstate() const {
    return pybind11::make_tuple(
        // unique identifying string for the force object
        "PressureForce",
        // state
        parameters);
}

// load state
template<> void
VertexModel::addVertexForce<PressureForce, pybind11::tuple const&>(
    std::string const& name, pybind11::tuple const& t) {
    // check
    checkSize(t, 2);
    assert(t[0].cast<std::string>() == "PressureForce");
    // initialise force
    ParametersType const parameters = t[1].cast<ParametersType>();
    addVertexForce<PressureForce, double const&, bool const&>(
        name, parameters.at("F"), parameters.at("fixedForce"));
}

/*
 *  GrowingAreaPerimeterForce
 *
 */

// save state
pybind11::tuple GrowingAreaPerimeterForce::pybind11_getstate() const {
    return pybind11::make_tuple(
        // unique identifying string for the force object
        "GrowingAreaPerimeterForce",
        // state
        parameters, targetArea);
}

// load state
template<> void
VertexModel::addVertexForce<GrowingAreaPerimeterForce, pybind11::tuple const&>(
    std::string const& name, pybind11::tuple const& t) {
    // check
    checkSize(t, 3);
    assert(t[0].cast<std::string>() == "GrowingAreaPerimeterForce");
    // initialise force
    ParametersType const parameters = t[1].cast<ParametersType>();
    addVertexForce<GrowingAreaPerimeterForce,
        double const&, double const&,
        double const&, double const&>(
        name,
        parameters.at("kA"), parameters.at("s0"),
        parameters.at("A0"), parameters.at("tauA"));
    // set internal degrees of freedom state
    std::map<long int, double> const targetArea =
        t[2].cast<std::map<long int, double>>();
    std::shared_ptr<GrowingAreaPerimeterForce> gapf =
        std::static_pointer_cast<GrowingAreaPerimeterForce>(
            vertexForces[name]);
    gapf->setTargetArea(targetArea);
}

/*
 *  TensionForce
 *
 */

// save state
pybind11::tuple TensionForce::pybind11_getstate() const {
    return pybind11::make_tuple(
        // unique identifying string for the force object
        "TensionForce",
        // state
        parameters);
}

// load state
template<> void
VertexModel::addHalfEdgeForce<TensionForce, pybind11::tuple const&>(
    std::string const& name, pybind11::tuple const& t) {
    // check
    checkSize(t, 2);
    assert(t[0].cast<std::string>() == "TensionForce");
    // initialise force
    ParametersType const parameters = t[1].cast<ParametersType>();
    addHalfEdgeForce<TensionForce, double const&>(
        name, parameters.at("Lambda"));
}

/*
 *  BoundaryTension
 *
 */

// save state
pybind11::tuple BoundaryTension::pybind11_getstate() const {
    return pybind11::make_tuple(
        // unique identifying string for the force object
        "BoundaryTension",
        // state
        parameters);
}

// load state
template<> void
VertexModel::addVertexForce<BoundaryTension, pybind11::tuple const&>(
    std::string const& name, pybind11::tuple const& t) {
    // check
    checkSize(t, 2);
    assert(t[0].cast<std::string>() == "BoundaryTension");
    // initialise force
    ParametersType const parameters = t[1].cast<ParametersType>();
    addVertexForce<BoundaryTension, double const&>(
        name, parameters.at("gamma"));
}

/*
 *  ActiveBrownianForce
 *
 */

// save state
pybind11::tuple ActiveBrownianForce::pybind11_getstate() const {
    return pybind11::make_tuple(
        // unique identifying string for the force object
        "ActiveBrownianForce",
        // state
        parameters, theta);
}

// load state
template<> void
VertexModel::addVertexForce<ActiveBrownianForce, pybind11::tuple const&>(
    std::string const& name, pybind11::tuple const& t) {
    // check
    checkSize(t, 3);
    assert(t[0].cast<std::string>() == "ActiveBrownianForce");
    // initialise force
    ParametersType const parameters =
        t[1].cast<ParametersType>();
    addVertexForce<ActiveBrownianForce, double const&, double const&>(
        name, parameters.at("v0"), parameters.at("taup"));
    // set internal degrees of freedom state
    std::map<long int, double> const theta =
        t[2].cast<std::map<long int, double>>();
    std::shared_ptr<ActiveBrownianForce> abp =
        std::static_pointer_cast<ActiveBrownianForce>(
            vertexForces[name]);
    abp->setTheta(theta);
}

/*
 *  ActiveBrownianCellForce
 *
 */

// save state
pybind11::tuple ActiveBrownianCellForce::pybind11_getstate() const {
    return pybind11::make_tuple(
        // unique identifying string for the force object
        "ActiveBrownianCellForce",
        // state
        parameters, theta);
}

// load state
template<> void
VertexModel::addVertexForce<ActiveBrownianCellForce, pybind11::tuple const&>(
    std::string const& name, pybind11::tuple const& t) {
    // check
    checkSize(t, 3);
    assert(t[0].cast<std::string>() == "ActiveBrownianCellForce");
    // initialise force
    ParametersType const parameters =
        t[1].cast<ParametersType>();
    addVertexForce<ActiveBrownianCellForce, double const&, double const&>(
        name, parameters.at("v0"), parameters.at("taup"));
    // set internal degrees of freedom state
    std::map<long int, double> const theta =
        t[2].cast<std::map<long int, double>>();
    std::shared_ptr<ActiveBrownianCellForce> abcp =
        std::static_pointer_cast<ActiveBrownianCellForce>(
            vertexForces[name]);
    abcp->setTheta(theta);
}

/*
 *  OrnsteinUhlenbeckTension
 *
 */

// save state
pybind11::tuple OrnsteinUhlenbeckTension::pybind11_getstate() const {
    return pybind11::make_tuple(
        // unique identifying string for the force object
        "OrnsteinUhlenbeckTension",
        // state
        parameters, tension);
}

// load state
template<> void
VertexModel::addHalfEdgeForce<OrnsteinUhlenbeckTension,
    pybind11::tuple const&>(
    std::string const& name, pybind11::tuple const& t) {
    // check
    checkSize(t, 3);
    assert(t[0].cast<std::string>() == "OrnsteinUhlenbeckTension");
    // initialise force
    ParametersType const parameters =
        t[1].cast<ParametersType>();
    addHalfEdgeForce<OrnsteinUhlenbeckTension,
        double const&, double const&, double const&>(
        name,
        parameters.at("t0"), parameters.at("st0"), parameters.at("taup"));
    // set internal degrees of freedom state
    std::map<long int, double> const tension =
        t[2].cast<std::map<long int, double>>();
    std::shared_ptr<OrnsteinUhlenbeckTension> out =
        std::static_pointer_cast<OrnsteinUhlenbeckTension>(
            halfEdgeForces[name]);
    out->setTension(tension);
}

/*
 *  Model0
 *
 */

// save state
pybind11::tuple Model0::pybind11_getstate() const {
    return pybind11::make_tuple(
        // unique identifying string for the force object
        "Model0",
        // state
        parameters, perimeter, noise, tension);
}

// load state
template<> void
VertexModel::addHalfEdgeForce<Model0,
    pybind11::tuple const&>(
    std::string const& name, pybind11::tuple const& t) {
    // check
    checkSize(t, 5);
    assert(t[0].cast<std::string>() == "Model0");
    // initialise force
    ParametersType const parameters =
        t[1].cast<ParametersType>();
    addHalfEdgeForce<Model0,
        double const&, double const&, double const&, double const&>(
        name,
        parameters.at("Gamma"), parameters.at("P0"),
        parameters.at("sigma"), parameters.at("taup"));
    // set internal degrees of freedom state
    std::map<long int, double> const perimeter =
        t[2].cast<std::map<long int, double>>();
    std::map<long int, double> const noise =
        t[3].cast<std::map<long int, double>>();
    std::map<long int, double> const tension =
        t[4].cast<std::map<long int, double>>();
    std::shared_ptr<Model0> m0 =
        std::static_pointer_cast<Model0>(
            halfEdgeForces[name]);
    m0->setPerimeter(perimeter);
    m0->setNoise(noise);
    m0->setTension(tension);
}

/*
 *  Model1
 *
 */

// save state
pybind11::tuple Model1::pybind11_getstate() const {
    return pybind11::make_tuple(
        // unique identifying string for the force object
        "Model1",
        // state
        parameters, perimeter, tension);
}

// load state
template<> void
VertexModel::addHalfEdgeForce<Model1,
    pybind11::tuple const&>(
    std::string const& name, pybind11::tuple const& t) {
    // check
    checkSize(t, 4);
    assert(t[0].cast<std::string>() == "Model1");
    // initialise force
    ParametersType const parameters =
        t[1].cast<ParametersType>();
    addHalfEdgeForce<Model1,
        double const&, double const&, double const&, double const&>(
        name,
        parameters.at("Gamma"), parameters.at("P0"),
        parameters.at("sigma"), parameters.at("taup"));
    // set internal degrees of freedom state
    std::map<long int, double> const perimeter =
        t[2].cast<std::map<long int, double>>();
    std::map<long int, double> const tension =
        t[3].cast<std::map<long int, double>>();
    std::shared_ptr<Model1> m1 =
        std::static_pointer_cast<Model1>(
            halfEdgeForces[name]);
    m1->setPerimeter(perimeter);
    m1->setTension(tension);
}

/*
 *  Model2
 *
 */

// save state
pybind11::tuple Model2::pybind11_getstate() const {
    return pybind11::make_tuple(
        // unique identifying string for the force object
        "Model2",
        // state
        parameters, length, restLength, noise, tension);
}

// load state
template<> void
VertexModel::addHalfEdgeForce<Model2,
    pybind11::tuple const&>(
    std::string const& name, pybind11::tuple const& t) {
    // check
    checkSize(t, 6);
    assert(t[0].cast<std::string>() == "Model2");
    // initialise force
    ParametersType const parameters =
        t[1].cast<ParametersType>();
    addHalfEdgeForce<Model2,
        double const&, double const&, double const&, double const&>(
        name,
        parameters.at("Gamma"), parameters.at("taur"),
        parameters.at("sigma"), parameters.at("taup"));
    // set internal degrees of freedom state
    std::map<long int, double> const length =
        t[2].cast<std::map<long int, double>>();
    std::map<long int, double> const restLength =
        t[3].cast<std::map<long int, double>>();
    std::map<long int, double> const noise =
        t[4].cast<std::map<long int, double>>();
    std::map<long int, double> const tension =
        t[5].cast<std::map<long int, double>>();
    std::shared_ptr<Model2> m2 =
        std::static_pointer_cast<Model2>(
            halfEdgeForces[name]);
    m2->setLength(length);
    m2->setRestLength(restLength);
    m2->setNoise(noise);
    m2->setTension(tension);
}

/*
 *  Model3
 *
 */

// save state
pybind11::tuple Model3::pybind11_getstate() const {
    return pybind11::make_tuple(
        // unique identifying string for the force object
        "Model3",
        // state
        parameters, length, restLength, tension);
}

// load state
template<> void
VertexModel::addHalfEdgeForce<Model3,
    pybind11::tuple const&>(
    std::string const& name, pybind11::tuple const& t) {
    // check
    checkSize(t, 5);
    assert(t[0].cast<std::string>() == "Model3");
    // initialise force
    ParametersType const parameters =
        t[1].cast<ParametersType>();
    addHalfEdgeForce<Model3,
        double const&, double const&, double const&>(
        name,
        parameters.at("Gamma"),
        parameters.at("sigma"), parameters.at("taup"));
    // set internal degrees of freedom state
    std::map<long int, double> const length =
        t[2].cast<std::map<long int, double>>();
    std::map<long int, double> const restLength =
        t[3].cast<std::map<long int, double>>();
    std::map<long int, double> const tension =
        t[4].cast<std::map<long int, double>>();
    std::shared_ptr<Model3> m3 =
        std::static_pointer_cast<Model3>(
            halfEdgeForces[name]);
    m3->setLength(length);
    m3->setRestLength(restLength);
    m3->setTension(tension);
}

/*
 *  Model4
 *
 */

// save state
pybind11::tuple Model4::pybind11_getstate() const {
    return pybind11::make_tuple(
        // unique identifying string for the force object
        "Model4",
        // state
        parameters, length, restLength, tension);
}

// load state
template<> void
VertexModel::addHalfEdgeForce<Model4,
    pybind11::tuple const&>(
    std::string const& name, pybind11::tuple const& t) {
    // check
    checkSize(t, 5);
    assert(t[0].cast<std::string>() == "Model4");
    // initialise force
    ParametersType const parameters =
        t[1].cast<ParametersType>();
    addHalfEdgeForce<Model4,
        double const&, double const&, double const&, double const&>(
        name,
        parameters.at("Gamma"), parameters.at("taur"),
        parameters.at("sigma"), parameters.at("taup"));
    // set internal degrees of freedom state
    std::map<long int, double> const length =
        t[2].cast<std::map<long int, double>>();
    std::map<long int, double> const restLength =
        t[3].cast<std::map<long int, double>>();
    std::map<long int, double> const tension =
        t[4].cast<std::map<long int, double>>();
    std::shared_ptr<Model4> m4 =
        std::static_pointer_cast<Model4>(
            halfEdgeForces[name]);
    m4->setLength(length);
    m4->setRestLength(restLength);
    m4->setTension(tension);
}

/*
 *  KeratinModel
 *
 */

// save state
pybind11::tuple KeratinModel::pybind11_getstate() const {
    return pybind11::make_tuple(
        // unique identifying string for the force object
        "KeratinModel",
        // state
        parameters, keratin, targetArea);
}

// load state
template<> void
VertexModel::addVertexForce<KeratinModel,
    pybind11::tuple const&>(
    std::string const& name, pybind11::tuple const& t) {
    // check
    checkSize(t, 4);
    assert(t[0].cast<std::string>() == "KeratinModel");
    // initialise force
    ParametersType const parameters =
        t[1].cast<ParametersType>();
    addVertexForce<KeratinModel,
        double const&, double const&, double const&,
        double const&, double const&,
        double const&, double const&, double const&,
        double const&, double const&, double const&>(
        name,
        parameters.at("K"), parameters.at("A0"), parameters.at("taur"),
        parameters.at("Gamma"), parameters.at("p0"),
        parameters.at("alpha"), parameters.at("beta"), parameters.at("kth"),
        parameters.at("tau"), parameters.at("sigma"), parameters.at("ron"));
    // set internal degrees of freedom state
    std::shared_ptr<KeratinModel> k =
        std::static_pointer_cast<KeratinModel>(
            vertexForces[name]);
    std::map<long int, double> const keratin =
        t[2].cast<std::map<long int, double>>();
    k->setKeratin(keratin);
    std::map<long int, double> const targetArea =
        t[3].cast<std::map<long int, double>>();
    k->setTargetArea(targetArea);
}

/*
 *  ClassFactory<VertexForce<ForcesType>>
 *
 */

// template specialisation
template<> void
pybind11_setstate_force_class_factory<VertexForce<ForcesType>>(
    VertexModel& vm, std::map<std::string, pybind11::tuple> const& stateMap) {
    for (auto it=stateMap.begin(); it != stateMap.end(); ++it) {
        std::string const forceName = (it->second)[0].cast<std::string>();
        // templated lambda function to add force to VertexModel
        auto addVertexForce = [&vm, &it]<class Force>() {
            vm.addVertexForce<Force, pybind11::tuple const&>(
                it->first, it->second);
        };
        // --- MANUAL ASSOCIATION TO FORCES ---
        if (forceName == "PerimeterForce") {
            addVertexForce.template operator()<PerimeterForce>();
        }
        else if (forceName == "AreaForce") {
            addVertexForce.template operator()<AreaForce>();
        }
        else if (forceName == "SurfaceForce") {
            addVertexForce.template operator()<SurfaceForce>();
        }
        else if (forceName == "VolumeForce") {
            addVertexForce.template operator()<VolumeForce>();
        }
        else if (forceName == "LinearVolumeForce") {
            addVertexForce.template operator()<LinearVolumeForce>();
        }
        else if (forceName == "PressureForce") {
            addVertexForce.template operator()<PressureForce>();
        }
        else if (forceName == "BoundaryTension") {
            addVertexForce.template operator()<BoundaryTension>();
        }
        else if (forceName == "GrowingAreaPerimeterForce") {
            addVertexForce.template operator()<GrowingAreaPerimeterForce>();
        }
        else if (forceName == "ActiveBrownianForce") {
            addVertexForce.template operator()<ActiveBrownianForce>();
        }
        else if (forceName == "ActiveBrownianCellForce") {
            addVertexForce.template operator()<ActiveBrownianCellForce>();
        }
        else if (forceName == "KeratinModel") {
            addVertexForce.template operator()<KeratinModel>();
        }
        // throw error if force not recognised
        else {
            throw std::runtime_error(
                "Force object '" + forceName + "' does not exist.");
        }
    }
}

/*
 *  ClassFactory<HalfEdgeForce<ForcesType>>
 *
 */

// template specialisation
template<> void
pybind11_setstate_force_class_factory<HalfEdgeForce<ForcesType>>(
    VertexModel& vm, std::map<std::string, pybind11::tuple> const& stateMap) {
    for (auto it=stateMap.begin(); it != stateMap.end(); ++it) {
        std::string const forceName = (it->second)[0].cast<std::string>();
        // templated lambda function to add force to VertexModel
        auto addHalfEdgeForce = [&vm, &it]<class Force>() {
            vm.addHalfEdgeForce<Force, pybind11::tuple const&>(
                it->first, it->second);
        };
        // --- MANUAL ASSOCIATION TO FORCES ---
        if (forceName == "TensionForce") {
            addHalfEdgeForce.template operator()<TensionForce>();
        }
        else if (forceName == "OrnsteinUhlenbeckTension") {
            addHalfEdgeForce.template operator()<OrnsteinUhlenbeckTension>();
        }
        else if (forceName == "Model0") {
            addHalfEdgeForce.template operator()<Model0>();
        }
        else if (forceName == "Model1") {
            addHalfEdgeForce.template operator()<Model1>();
        }
        else if (forceName == "Model2") {
            addHalfEdgeForce.template operator()<Model2>();
        }
        else if (forceName == "Model3") {
            addHalfEdgeForce.template operator()<Model3>();
        }
        else if (forceName == "Model4") {
            addHalfEdgeForce.template operator()<Model4>();
        }
        // throw error if force not recognised
        else {
            throw std::runtime_error(
                "Force object '" + forceName + "' does not exist.");
        }
    }
}

#endif

