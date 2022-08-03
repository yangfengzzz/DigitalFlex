//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.force/fluid_model.h"

#include <iostream>

#include "vox.force/drag_solver.h"
#include "vox.force/elasticity_solver.h"
#include "vox.force/surface_tension_solver.h"
#include "vox.force/viscosity_solver.h"
#include "vox.force/vorticity_solver.h"
//#include "EmitterSystem.h"
//#include "NeighborhoodSearch.h"
#include "vox.force/simulation.h"
//#include "SPHKernels.h"
//#include "TimeManager.h"
//#include "TimeStep.h"
//#include "Utilities/Logger.h"

namespace vox::flex {

int FluidModel::NUM_PARTICLES = -1;
int FluidModel::NUM_REUSED_PARTICLES = -1;
int FluidModel::DENSITY0 = -1;
int FluidModel::DRAG_METHOD = -1;
int FluidModel::SURFACE_TENSION_METHOD = -1;
int FluidModel::VISCOSITY_METHOD = -1;
int FluidModel::VORTICITY_METHOD = -1;
int FluidModel::ELASTICITY_METHOD = -1;

FluidModel::FluidModel()
    : m_masses(),
      m_a(),
      m_v0(),
      m_x0(),
      m_x(),
      m_v(),
      m_density(),
      m_particleId(),
      m_objectId(),
      m_objectId0(),
      m_particleState() {
    m_density0 = 1000.0;
    m_pointSetIndex = 0;

    m_emitterSystem = new EmitterSystem(this);
    m_viscosity = nullptr;
    m_viscosityMethod = 0;
    m_surfaceTension = nullptr;
    m_surfaceTensionMethod = 0;
    m_vorticityMethod = 0;
    m_vorticity = nullptr;
    m_dragMethod = 0;
    m_drag = nullptr;
    m_dragMethodChanged = nullptr;
    m_surfaceTensionMethodChanged = nullptr;
    m_viscosityMethodChanged = nullptr;
    m_vorticityMethodChanged = nullptr;
    m_elasticityMethod = 0;
    m_elasticity = nullptr;
    m_elasticityMethodChanged = nullptr;

    addField({"id", FieldType::UInt, [&](const unsigned int i) -> unsigned int* { return &getParticleId(i); }, true});
    addField({"state", FieldType::UInt,
              [&](const unsigned int i) -> unsigned int* { return (unsigned int*)(&m_particleState[i]); }, true});
    addField({"object_id", FieldType::UInt, [&](const unsigned int i) -> unsigned int* { return &getObjectId(i); },
              true});
    addField({"position", FieldType::Vector3, [&](const unsigned int i) -> double* { return &getPosition(i)[0]; },
              true});
    addField({"position0", FieldType::Vector3, [&](const unsigned int i) -> double* { return &getPosition0(i)[0]; }});
    addField({"velocity", FieldType::Vector3, [&](const unsigned int i) -> double* { return &getVelocity(i)[0]; },
              true});
    addField({"density", FieldType::Scalar, [&](const unsigned int i) -> double* { return &getDensity(i); }, false});
}

FluidModel::~FluidModel() {
    removeFieldByName("id");
    removeFieldByName("state");
    removeFieldByName("object_id");
    removeFieldByName("position");
    removeFieldByName("position0");
    removeFieldByName("velocity");
    removeFieldByName("density");

    delete m_emitterSystem;
    delete m_surfaceTension;
    delete m_drag;
    delete m_vorticity;
    delete m_viscosity;
    delete m_elasticity;

    releaseFluidParticles();
}

void FluidModel::init() {
    initParameters();

    setViscosityMethod(1);
}

void FluidModel::deferredInit() {
    if (m_surfaceTension) m_surfaceTension->deferredInit();
    if (m_viscosity) m_viscosity->deferredInit();
    if (m_vorticity) m_vorticity->deferredInit();
    if (m_drag) m_drag->deferredInit();
    if (m_elasticity) m_elasticity->deferredInit();
}

void FluidModel::initParameters() {}

void FluidModel::reset() {
    setNumActiveParticles(m_numActiveParticles0);
    const unsigned int nPoints = numActiveParticles();

    struct Comparator {
        Comparator(FluidModel* _this) : m_this(_this){};
        bool operator()(unsigned int a, unsigned int b) { return m_this->getParticleId(a) < m_this->getParticleId(b); }
        FluidModel* m_this;
    };

    // Fluid
    for (unsigned int i = 0; i < nPoints; i++) {
        const Vector3D& x0 = getPosition0(i);
        getPosition(i) = x0;
        getVelocity(i) = getVelocity0(i);
        getAcceleration(i).setZero();
        m_objectId[i] = m_objectId0[i];
        m_density[i] = 0.0;
        m_particleId[i] = i;
        m_particleState[i] = ParticleState::Active;
    }
    // emitted particles
    for (unsigned int i = nPoints; i < (unsigned int)m_particleId.size(); i++) {
        m_particleId[i] = i;
    }

    NeighborhoodSearch* neighborhoodSearch = Simulation::getCurrent()->getNeighborhoodSearch();
    if (neighborhoodSearch->point_set(m_pointSetIndex).n_points() != nPoints)
        neighborhoodSearch->resize_point_set(m_pointSetIndex, &getPosition(0)[0], nPoints);

    if (m_surfaceTension) m_surfaceTension->reset();
    if (m_viscosity) m_viscosity->reset();
    if (m_vorticity) m_vorticity->reset();
    if (m_drag) m_drag->reset();
    if (m_elasticity) m_elasticity->reset();

    m_emitterSystem->reset();
}

void FluidModel::initMasses() {
    const double particleRadius = Simulation::getCurrent()->getParticleRadius();
    const int nParticles = (int)numParticles();
    const double diam = static_cast<double>(2.0) * particleRadius;
    if (Simulation::getCurrent()->is2DSimulation())
        m_V = static_cast<double>(0.8) * diam * diam;
    else
        m_V = static_cast<double>(0.8) * diam * diam * diam;

#pragma omp parallel default(shared)
    {
#pragma omp for schedule(static)
        for (int i = 0; i < nParticles; i++) {
            setMass(i, m_V * m_density0);
            // each particle represents a cube with a side length of r
            // mass is slightly reduced to prevent pressure at the beginning of the simulation
        }
    }
}

void FluidModel::resizeFluidParticles(const unsigned int newSize) {
    m_x0.resize(newSize);
    m_x.resize(newSize);
    m_v.resize(newSize);
    m_v0.resize(newSize);
    m_a.resize(newSize);
    m_masses.resize(newSize);
    m_density.resize(newSize);
    m_particleId.resize(newSize);
    m_objectId.resize(newSize);
    m_objectId0.resize(newSize);
    m_particleState.resize(newSize, ParticleState::Active);
}

void FluidModel::releaseFluidParticles() {
    m_x0.clear();
    m_x.clear();
    m_v.clear();
    m_v0.clear();
    m_a.clear();
    m_masses.clear();
    m_density.clear();
    m_particleId.clear();
    m_objectId.clear();
    m_objectId0.clear();
    m_particleState.clear();
}

void FluidModel::initModel(const std::string& id,
                           const unsigned int nFluidParticles,
                           Vector3D* fluidParticles,
                           Vector3D* fluidVelocities,
                           unsigned int* fluidObjectIds,
                           const unsigned int nMaxEmitterParticles) {
    m_id = id;
    init();
    releaseFluidParticles();
    resizeFluidParticles(nFluidParticles + nMaxEmitterParticles);

// copy fluid positions
#pragma omp parallel default(shared)
    {
#pragma omp for schedule(static)
        for (int i = 0; i < (int)nFluidParticles; i++) {
            getPosition0(i) = fluidParticles[i];
            getPosition(i) = fluidParticles[i];
            getVelocity0(i) = fluidVelocities[i];
            getVelocity(i) = fluidVelocities[i];
            getAcceleration(i).setZero();
            m_density[i] = 0.0;
            m_particleId[i] = i;
            m_objectId[i] = fluidObjectIds[i];
            m_objectId0[i] = fluidObjectIds[i];
            if (m_particleState[i] != ParticleState::Fixed) m_particleState[i] = ParticleState::Active;
        }
    }
    // set IDs for emitted particles
    for (unsigned int i = nFluidParticles; i < (nFluidParticles + nMaxEmitterParticles); i++) {
        m_particleId[i] = i;
        m_objectId[i] = 0;
    }

    // initialize masses
    initMasses();

    // Fluids
    NeighborhoodSearch* neighborhoodSearch = Simulation::getCurrent()->getNeighborhoodSearch();
    m_pointSetIndex = neighborhoodSearch->add_point_set(&getPosition(0)[0], nFluidParticles, true, true, true, this);

    m_numActiveParticles0 = nFluidParticles;
    m_numActiveParticles = m_numActiveParticles0;
}

void FluidModel::performNeighborhoodSearchSort() {
    const unsigned int numPart = numActiveParticles();
    if (numPart == 0) return;

    NeighborhoodSearch* neighborhoodSearch = Simulation::getCurrent()->getNeighborhoodSearch();

    auto const& d = neighborhoodSearch->point_set(m_pointSetIndex);
    d.sort_field(&m_x[0]);
    d.sort_field(&m_v[0]);
    d.sort_field(&m_a[0]);
    d.sort_field(&m_masses[0]);
    d.sort_field(&m_density[0]);
    d.sort_field(&m_particleId[0]);
    d.sort_field(&m_objectId[0]);
    d.sort_field(&m_particleState[0]);

    if (m_viscosity) m_viscosity->performNeighborhoodSearchSort();
    if (m_surfaceTension) m_surfaceTension->performNeighborhoodSearchSort();
    if (m_vorticity) m_vorticity->performNeighborhoodSearchSort();
    if (m_drag) m_drag->performNeighborhoodSearchSort();
    if (m_elasticity) m_elasticity->performNeighborhoodSearchSort();
}

void FluidModel::setDensity0(const double v) {
    m_density0 = v;
    initMasses();
}

const FieldDescription& FluidModel::getField(const std::string& name) {
    unsigned int index = 0;
    for (auto i = 0; i < m_fields.size(); i++) {
        if (m_fields[i].name == name) {
            index = i;
            break;
        }
    }
    return m_fields[index];
}

void FluidModel::setNumActiveParticles(const unsigned int num) { m_numActiveParticles = num; }

unsigned int FluidModel::numActiveParticles() const { return m_numActiveParticles; }

void FluidModel::setDragMethodChangedCallback(std::function<void()> const& callBackFct) {
    m_dragMethodChanged = callBackFct;
}

void FluidModel::setSurfaceMethodChangedCallback(std::function<void()> const& callBackFct) {
    m_surfaceTensionMethodChanged = callBackFct;
}

void FluidModel::setViscosityMethodChangedCallback(std::function<void()> const& callBackFct) {
    m_viscosityMethodChanged = callBackFct;
}

void FluidModel::setVorticityMethodChangedCallback(std::function<void()> const& callBackFct) {
    m_vorticityMethodChanged = callBackFct;
}

void FluidModel::setElasticityMethodChangedCallback(std::function<void()> const& callBackFct) {
    m_elasticityMethodChanged = callBackFct;
}

void FluidModel::computeSurfaceTension() {
    if (m_surfaceTension) m_surfaceTension->step();
}

void FluidModel::computeViscosity() {
    if (m_viscosity) m_viscosity->step();
}

void FluidModel::computeVorticity() {
    if (m_vorticity) m_vorticity->step();
}

void FluidModel::computeDragForce() {
    if (m_drag) m_drag->step();
}

void FluidModel::computeElasticity() {
    if (m_elasticity) m_elasticity->step();
}

void FluidModel::emittedParticles(const unsigned int startIndex) {
    if (m_viscosity) m_viscosity->emittedParticles(startIndex);
    if (m_surfaceTension) m_surfaceTension->emittedParticles(startIndex);
    if (m_vorticity) m_vorticity->emittedParticles(startIndex);
    if (m_drag) m_drag->emittedParticles(startIndex);
    if (m_elasticity) m_elasticity->emittedParticles(startIndex);
}

void FluidModel::setSurfaceTensionMethod(const std::string& val) {
    Simulation* sim = Simulation::getCurrent();
    std::vector<Simulation::NonPressureForceMethod>& methods = sim->getSurfaceTensionMethods();
    for (size_t i = 0; i < methods.size(); i++) {
        if (methods[i].m_name == val) {
            setSurfaceTensionMethod(static_cast<unsigned int>(i));
            break;
        }
    }
}

void FluidModel::setSurfaceTensionMethod(const unsigned int val) {
    Simulation* sim = Simulation::getCurrent();
    std::vector<Simulation::NonPressureForceMethod>& stMethods = sim->getSurfaceTensionMethods();

    unsigned int stm = val;
    if (stm >= stMethods.size()) stm = 0;

    if (stm == m_surfaceTensionMethod) return;

    delete m_surfaceTension;
    m_surfaceTension = nullptr;

    m_surfaceTensionMethod = stm;
    int method_id = getValue<int>(FluidModel::SURFACE_TENSION_METHOD);
    for (unsigned int i = 0; i < stMethods.size(); i++) {
        if (stMethods[i].m_id == method_id)
            m_surfaceTension = static_cast<SurfaceTensionBase*>(stMethods[i].m_creator(this));
    }

    if (m_surfaceTension != nullptr) m_surfaceTension->init();

    if (m_surfaceTensionMethodChanged != nullptr) m_surfaceTensionMethodChanged();
}

void FluidModel::setViscosityMethod(const std::string& val) {
    Simulation* sim = Simulation::getCurrent();
    std::vector<Simulation::NonPressureForceMethod>& methods = sim->getViscosityMethods();
    for (size_t i = 0; i < methods.size(); i++) {
        if (methods[i].m_name == val) {
            setViscosityMethod(static_cast<unsigned int>(i));
            break;
        }
    }
}

void FluidModel::setViscosityMethod(const unsigned int val) {
    Simulation* sim = Simulation::getCurrent();
    std::vector<Simulation::NonPressureForceMethod>& viscoMethods = sim->getViscosityMethods();

    unsigned int vm = val;
    if (vm >= viscoMethods.size()) vm = 0;

    if (vm == m_viscosityMethod) return;

    delete m_viscosity;
    m_viscosity = nullptr;

    m_viscosityMethod = vm;

    for (unsigned int i = 0; i < viscoMethods.size(); i++) {
        if (viscoMethods[i].m_id == m_viscosityMethod) {
            m_viscosity = static_cast<ViscosityBase*>(viscoMethods[i].m_creator(this));
            break;
        }
    }

    if (m_viscosity != nullptr) m_viscosity->init();

    if (m_viscosityMethodChanged != nullptr) m_viscosityMethodChanged();
}

void FluidModel::setVorticityMethod(const std::string& val) {
    Simulation* sim = Simulation::getCurrent();
    std::vector<Simulation::NonPressureForceMethod>& methods = sim->getVorticityMethods();
    for (size_t i = 0; i < methods.size(); i++) {
        if (methods[i].m_name == val) {
            setVorticityMethod(static_cast<unsigned int>(i));
            break;
        }
    }
}

void FluidModel::setVorticityMethod(const unsigned int val) {
    Simulation* sim = Simulation::getCurrent();
    std::vector<Simulation::NonPressureForceMethod>& vorticityMethods = sim->getVorticityMethods();

    unsigned int vm = val;
    if (vm >= vorticityMethods.size()) vm = 0;

    if (vm == m_vorticityMethod) return;

    delete m_vorticity;
    m_vorticity = nullptr;

    m_vorticityMethod = vm;

    int method_id = getValue<int>(FluidModel::VORTICITY_METHOD);
    for (unsigned int i = 0; i < vorticityMethods.size(); i++) {
        if (vorticityMethods[i].m_id == method_id)
            m_vorticity = static_cast<VorticityBase*>(vorticityMethods[i].m_creator(this));
    }

    if (m_vorticity != nullptr) m_vorticity->init();

    if (m_vorticityMethodChanged != nullptr) m_vorticityMethodChanged();
}

void FluidModel::setDragMethod(const std::string& val) {
    Simulation* sim = Simulation::getCurrent();
    std::vector<Simulation::NonPressureForceMethod>& methods = sim->getDragMethods();
    for (size_t i = 0; i < methods.size(); i++) {
        if (methods[i].m_name == val) {
            setDragMethod(static_cast<unsigned int>(i));
            break;
        }
    }
}

void FluidModel::setDragMethod(const unsigned int val) {
    Simulation* sim = Simulation::getCurrent();
    std::vector<Simulation::NonPressureForceMethod>& dragMethods = sim->getDragMethods();

    unsigned int dm = val;
    if (dm >= dragMethods.size()) dm = 0;

    if (dm == m_dragMethod) return;

    delete m_drag;
    m_drag = nullptr;

    m_dragMethod = dm;

    int method_id = getValue<int>(FluidModel::DRAG_METHOD);
    for (unsigned int i = 0; i < dragMethods.size(); i++) {
        if (dragMethods[i].m_id == method_id) m_drag = static_cast<DragBase*>(dragMethods[i].m_creator(this));
    }

    if (m_drag != nullptr) m_drag->init();

    if (m_dragMethodChanged != nullptr) m_dragMethodChanged();
}

void FluidModel::setElasticityMethod(const std::string& val) {
    Simulation* sim = Simulation::getCurrent();
    std::vector<Simulation::NonPressureForceMethod>& methods = sim->getElasticityMethods();
    for (size_t i = 0; i < methods.size(); i++) {
        if (methods[i].m_name == val) {
            setElasticityMethod(static_cast<unsigned int>(i));
            break;
        }
    }
}

void FluidModel::setElasticityMethod(const unsigned int val) {
    Simulation* sim = Simulation::getCurrent();
    std::vector<Simulation::NonPressureForceMethod>& elasticityMethods = sim->getElasticityMethods();

    unsigned int em = val;
    if (em >= elasticityMethods.size()) em = 0;

    if (em == m_elasticityMethod) return;

    delete m_elasticity;
    m_elasticity = nullptr;

    m_elasticityMethod = em;

    int method_id = getValue<int>(FluidModel::ELASTICITY_METHOD);
    for (unsigned int i = 0; i < elasticityMethods.size(); i++) {
        if (elasticityMethods[i].m_id == method_id)
            m_elasticity = static_cast<ElasticityBase*>(elasticityMethods[i].m_creator(this));
    }

    if (m_elasticity != nullptr) m_elasticity->init();

    if (m_elasticityMethodChanged != nullptr) m_elasticityMethodChanged();
}

void FluidModel::addField(const FieldDescription& field) {
    m_fields.push_back(field);
    std::sort(m_fields.begin(), m_fields.end(),
              [](FieldDescription& i, FieldDescription& j) -> bool { return (i.name < j.name); });
}

void FluidModel::removeFieldByName(const std::string& fieldName) {
    for (auto it = m_fields.begin(); it != m_fields.end(); it++) {
        if (it->name == fieldName) {
            m_fields.erase(it);
            break;
        }
    }
}

void FluidModel::saveState(BinaryFileWriter& binWriter) {
    binWriter.write(m_numActiveParticles);
    binWriter.writeBuffer((char*)m_particleState.data(), m_numActiveParticles * sizeof(ParticleState));

    if (m_surfaceTension) m_surfaceTension->saveState(binWriter);
    if (m_viscosity) m_viscosity->saveState(binWriter);
    if (m_vorticity) m_vorticity->saveState(binWriter);
    if (m_drag) m_drag->saveState(binWriter);
    if (m_elasticity) m_elasticity->saveState(binWriter);
    m_emitterSystem->saveState(binWriter);
}

void FluidModel::loadState(BinaryFileReader& binReader) {
    binReader.read(m_numActiveParticles);
    NeighborhoodSearch* neighborhoodSearch = Simulation::getCurrent()->getNeighborhoodSearch();
    neighborhoodSearch->update_point_sets();
    neighborhoodSearch->resize_point_set(m_pointSetIndex, &getPosition(0)[0], m_numActiveParticles);

    binReader.readBuffer((char*)m_particleState.data(), m_numActiveParticles * sizeof(ParticleState));

    if (m_surfaceTension) m_surfaceTension->loadState(binReader);
    if (m_viscosity) m_viscosity->loadState(binReader);
    if (m_vorticity) m_vorticity->loadState(binReader);
    if (m_drag) m_drag->loadState(binReader);
    if (m_elasticity) m_elasticity->loadState(binReader);

    m_emitterSystem->loadState(binReader);
}

}  // namespace vox::flex