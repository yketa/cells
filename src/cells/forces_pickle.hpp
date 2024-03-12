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
 *  EdgePullForce
 *
 */

// save state
pybind11::tuple EdgePullForce::pybind11_getstate() const {
    return pybind11::make_tuple(
        // unique identifying string for the force object
        "EdgePullForce",
        // state
        parameters);
}

// load state
template<> void
VertexModel::addVertexForce<EdgePullForce, pybind11::tuple const&>(
    std::string const& name, pybind11::tuple const& t) {
    // check
    checkSize(t, 2);
    assert(t[0].cast<std::string>() == "EdgePullForce");
    // initialise force
    ParametersType const parameters = t[1].cast<ParametersType>();
    addVertexForce<EdgePullForce, double const&>(
        name, parameters.at("Fpull"));
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
        parameters, time0, keratin);
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
    double const time0 =
        t[2].cast<double>();
    addVertexForce<KeratinModel,
        double const&,
        double const&, double const&,
        double const&, double const&,
        double const&, double const&, double const&,
        double const&, double const&,
        double const&, double const&, double const&>(
        name,
        time0,
        parameters.at("K"), parameters.at("A0"),
        parameters.at("Gamma"), parameters.at("P0"),
        parameters.at("l0"), parameters.at("alpha"), parameters.at("kth"),
        parameters.at("tau"), parameters.at("sigma"),
        parameters.at("ron"), parameters.at("k0"), parameters.at("p0"));
    // set internal degrees of freedom state
    std::map<long int, double> const keratin =
        t[2].cast<std::map<long int, double>>();
    std::shared_ptr<KeratinModel> k =
        std::static_pointer_cast<KeratinModel>(
            vertexForces[name]);
    k->setKeratin(keratin);
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
        if ( forceName == "PerimeterForce" ) {
            addVertexForce.template operator()<PerimeterForce>();
        }
        else if ( forceName == "AreaForce" ) {
            addVertexForce.template operator()<AreaForce>();
        }
        else if ( forceName == "ActiveBrownianForce" ) {
            addVertexForce.template operator()<ActiveBrownianForce>();
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
        if ( forceName == "OrnsteinUhlenbeckTension" ) {
            addHalfEdgeForce.template operator()<OrnsteinUhlenbeckTension>();
        }
        else if ( forceName == "Model0" ) {
            addHalfEdgeForce.template operator()<Model0>();
        }
        else if ( forceName == "Model1" ) {
            addHalfEdgeForce.template operator()<Model1>();
        }
        else if ( forceName == "Model2" ) {
            addHalfEdgeForce.template operator()<Model2>();
        }
        else if ( forceName == "Model3" ) {
            addHalfEdgeForce.template operator()<Model3>();
        }
        else if ( forceName == "Model4" ) {
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

