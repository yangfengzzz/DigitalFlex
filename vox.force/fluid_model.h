//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <utility>
#include <vector>

#include "vox.force/rigid_body_object.h"
//#include "SPHKernels.h"
#ifdef USE_AVX
#include "SPlisHSPlasH/Utilities/AVX_math.h"
#endif
#include "vox.force/binary_file_reader_writer.h"

#ifdef USE_PERFORMANCE_OPTIMIZATION
// compute the value xj (empty in the optimized version)
#define compute_xj(fm_neighbor, pid)

// compute the value Vj (empty in the optimized version)
#define compute_Vj(fm_neighbor)

// compute the value Vj * gradW
#define compute_Vj_gradW() \
    const Vector3f8 &V_gradW = model->get_precomputed_V_gradW()[model->get_precomputed_indices()[i] + idx];

// compute the value Vj * gradW
#define compute_Vj_gradW_samephase() \
    const Vector3f8 &V_gradW = model->get_precomputed_V_gradW()[model->get_precomputed_indices_same_phase()[i] + j / 8];
#else
// compute the value xj
#define compute_xj(fm_neighbor, pid) \
    const Vector3f8 xj_avx =         \
            convertVec_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &fm_neighbor->getPosition(0), count);

// compute the value Vj
#define compute_Vj(fm_neighbor) const Scalarf8 Vj_avx = convert_zero(fm_neighbor->getVolume(0), count);

// compute the value Vj * gradW assuming that xj and Vj are already available
#define compute_Vj_gradW() const Vector3f8 &V_gradW = CubicKernel_AVX::gradW(xi_avx - xj_avx) * Vj_avx;

// compute the value Vj * gradW assuming that xj and Vj are already available
#define compute_Vj_gradW_samephase() const Vector3f8 &V_gradW = CubicKernel_AVX::gradW(xi_avx - xj_avx) * Vj_avx;
#endif

namespace vox::flex {
class TimeStep;
class ViscosityBase;
class SurfaceTensionBase;
class VorticityBase;
class DragBase;
class ElasticityBase;
class EmitterSystem;

enum FieldType { Scalar = 0, Vector3, Vector6, Matrix3, Matrix6, UInt };
struct FieldDescription {
    std::string name;
    FieldType type;
    // getFct(particleIndex)
    std::function<void *(const unsigned int)> getFct;
    bool storeData;

    FieldDescription(std::string n,
                     const FieldType &t,
                     std::function<void *(const unsigned int)> fct,
                     const bool s = false)
        : name(std::move(n)), type(t), getFct(std::move(fct)), storeData(s) {}
};

enum class ParticleState { Active = 0, AnimatedByEmitter, Fixed };

/** \brief The fluid model stores the particle and simulation information
 */
class FluidModel {
public:
    static int NUM_PARTICLES;
    static int NUM_REUSED_PARTICLES;
    static int DENSITY0;

    static int DRAG_METHOD;
    static int SURFACE_TENSION_METHOD;
    static int VISCOSITY_METHOD;
    static int VORTICITY_METHOD;
    static int ELASTICITY_METHOD;

    FluidModel();
    FluidModel(const FluidModel &) = delete;
    FluidModel &operator=(const FluidModel &) = delete;
    virtual ~FluidModel();

    void init();
    /** This function is called after the simulation scene is loaded and all
     * parameters are initialized. While reading a scene file several parameters
     * can change. The deferred init function should initialize all values which
     * depend on these parameters.
     */
    void deferredInit();

    [[nodiscard]] std::string getId() const { return m_id; }

protected:
    std::string m_id;
    EmitterSystem *m_emitterSystem;

    // Mass
    // If the mass is zero, the particle is static
    std::vector<double> m_masses;
    std::vector<Vector3D> m_a;
    std::vector<Vector3D> m_v0;
    std::vector<Vector3D> m_x0;
    std::vector<Vector3D> m_x;
    std::vector<Vector3D> m_v;
    std::vector<double> m_density;
    std::vector<unsigned int> m_particleId;
    std::vector<unsigned int> m_objectId;
    std::vector<unsigned int> m_objectId0;
    std::vector<ParticleState> m_particleState;
    double m_V{};

#ifdef USE_PERFORMANCE_OPTIMIZATION
    std::vector<Vector3f8, Eigen::aligned_allocator<Vector3f8>> m_precomp_V_gradW;
    std::vector<unsigned int> m_precompIndices;
    std::vector<unsigned int> m_precompIndicesSamePhase;
#endif

    unsigned int m_surfaceTensionMethod;
    SurfaceTensionBase *m_surfaceTension;
    unsigned int m_viscosityMethod;
    ViscosityBase *m_viscosity;
    unsigned int m_vorticityMethod;
    VorticityBase *m_vorticity;
    unsigned int m_dragMethod;
    DragBase *m_drag;
    unsigned int m_elasticityMethod;
    ElasticityBase *m_elasticity;
    std::vector<FieldDescription> m_fields;

    std::function<void()> m_dragMethodChanged;
    std::function<void()> m_surfaceTensionMethodChanged;
    std::function<void()> m_viscosityMethodChanged;
    std::function<void()> m_vorticityMethodChanged;
    std::function<void()> m_elasticityMethodChanged;

    double m_density0;
    unsigned int m_pointSetIndex;

    unsigned int m_numActiveParticles{};
    unsigned int m_numActiveParticles0{};

    virtual void initParameters();

    void initMasses();

    /** Resize the arrays containing the particle data.
     */
    virtual void resizeFluidParticles(unsigned int newSize);

    /** Release the arrays containing the particle data.
     */
    virtual void releaseFluidParticles();

public:
    [[nodiscard]] inline double getDensity0() const { return m_density0; }
    void setDensity0(double v);

    [[nodiscard]] unsigned int getPointSetIndex() const { return m_pointSetIndex; }

    void addField(const FieldDescription &field);
    const std::vector<FieldDescription> &getFields() { return m_fields; }
    const FieldDescription &getField(const unsigned int i) { return m_fields[i]; }
    const FieldDescription &getField(const std::string &name);
    unsigned int numberOfFields() { return static_cast<unsigned int>(m_fields.size()); }
    void removeFieldByName(const std::string &fieldName);

    void setNumActiveParticles(unsigned int num);
    [[nodiscard]] unsigned int numberOfParticles() const { return static_cast<unsigned int>(m_x.size()); }

    EmitterSystem *getEmitterSystem() { return m_emitterSystem; }

    virtual void reset();

    void performNeighborhoodSearchSort();

    void initModel(const std::string &id,
                   const unsigned int nFluidParticles,
                   Vector3D *fluidParticles,
                   Vector3D *fluidVelocities,
                   const unsigned int *fluidObjectIds,
                   const unsigned int nMaxEmitterParticles);

    const unsigned int numParticles() const { return static_cast<unsigned int>(m_masses.size()); }
    unsigned int numActiveParticles() const;

    unsigned int getNumActiveParticles0() const { return m_numActiveParticles0; }
    void setNumActiveParticles0(unsigned int val) { m_numActiveParticles0 = val; }

    void emittedParticles(const unsigned int startIndex);

    unsigned int getSurfaceTensionMethod() const { return m_surfaceTensionMethod; }
    void setSurfaceTensionMethod(const std::string &val);
    void setSurfaceTensionMethod(const unsigned int val);
    unsigned int getViscosityMethod() const { return m_viscosityMethod; }
    void setViscosityMethod(const std::string &val);
    void setViscosityMethod(const unsigned int val);
    unsigned int getVorticityMethod() const { return m_vorticityMethod; }
    void setVorticityMethod(const std::string &val);
    void setVorticityMethod(const unsigned int val);
    unsigned int getDragMethod() const { return m_dragMethod; }
    void setDragMethod(const std::string &val);
    void setDragMethod(const unsigned int val);
    unsigned int getElasticityMethod() const { return m_elasticityMethod; }
    void setElasticityMethod(const std::string &val);
    void setElasticityMethod(const unsigned int val);

    SurfaceTensionBase *getSurfaceTensionBase() { return m_surfaceTension; }
    ViscosityBase *getViscosityBase() { return m_viscosity; }
    VorticityBase *getVorticityBase() { return m_vorticity; }
    DragBase *getDragBase() { return m_drag; }
    ElasticityBase *getElasticityBase() { return m_elasticity; }

    void setDragMethodChangedCallback(std::function<void()> const &callBackFct);
    void setSurfaceMethodChangedCallback(std::function<void()> const &callBackFct);
    void setViscosityMethodChangedCallback(std::function<void()> const &callBackFct);
    void setVorticityMethodChangedCallback(std::function<void()> const &callBackFct);
    void setElasticityMethodChangedCallback(std::function<void()> const &callBackFct);

    void computeSurfaceTension();
    void computeViscosity();
    void computeVorticity();
    void computeDragForce();
    void computeElasticity();

    void saveState(BinaryFileWriter &binWriter);
    void loadState(BinaryFileReader &binReader);

#ifdef USE_PERFORMANCE_OPTIMIZATION
    inline std::vector<Vector3f8, Eigen::aligned_allocator<Vector3f8>> &get_precomputed_V_gradW() {
        return m_precomp_V_gradW;
    }
    inline std::vector<unsigned int> &get_precomputed_indices() { return m_precompIndices; }
    inline std::vector<unsigned int> &get_precomputed_indices_same_phase() { return m_precompIndicesSamePhase; }
#endif

    inline Vector3D &getPosition0(const unsigned int i) { return m_x0[i]; }

    inline const Vector3D &getPosition0(const unsigned int i) const { return m_x0[i]; }

    inline void setPosition0(const unsigned int i, const Vector3D &pos) { m_x0[i] = pos; }

    inline Vector3D &getPosition(const unsigned int i) { return m_x[i]; }

    inline const Vector3D &getPosition(const unsigned int i) const { return m_x[i]; }

    inline void setPosition(const unsigned int i, const Vector3D &pos) { m_x[i] = pos; }

    inline Vector3D &getVelocity(const unsigned int i) { return m_v[i]; }

    inline const Vector3D &getVelocity(const unsigned int i) const { return m_v[i]; }

    inline void setVelocity(const unsigned int i, const Vector3D &vel) { m_v[i] = vel; }

    inline Vector3D &getVelocity0(const unsigned int i) { return m_v0[i]; }

    inline const Vector3D &getVelocity0(const unsigned int i) const { return m_v0[i]; }

    inline void setVelocity0(const unsigned int i, const Vector3D &vel) { m_v0[i] = vel; }

    inline Vector3D &getAcceleration(const unsigned int i) { return m_a[i]; }

    inline const Vector3D &getAcceleration(const unsigned int i) const { return m_a[i]; }

    inline void setAcceleration(const unsigned int i, const Vector3D &accel) { m_a[i] = accel; }

    inline const double getMass(const unsigned int i) const { return m_masses[i]; }

    inline double &getMass(const unsigned int i) { return m_masses[i]; }

    inline void setMass(const unsigned int i, const double mass) { m_masses[i] = mass; }

    inline const double &getDensity(const unsigned int i) const { return m_density[i]; }

    inline double &getDensity(const unsigned int i) { return m_density[i]; }

    inline void setDensity(const unsigned int i, const double &val) { m_density[i] = val; }

    inline unsigned int &getParticleId(const unsigned int i) { return m_particleId[i]; }

    inline const unsigned int &getParticleId(const unsigned int i) const { return m_particleId[i]; }

    inline unsigned int &getObjectId(const unsigned int i) { return m_objectId[i]; }

    inline const unsigned int &getObjectId(const unsigned int i) const { return m_objectId[i]; }

    inline void setObjectId(const unsigned int i, const unsigned int val) { m_objectId[i] = val; }

    inline const ParticleState &getParticleState(const unsigned int i) const { return m_particleState[i]; }

    inline ParticleState &getParticleState(const unsigned int i) { return m_particleState[i]; }

    inline void setParticleState(const unsigned int i, const ParticleState &val) { m_particleState[i] = val; }

    inline const double getVolume(const unsigned int i) const { return m_V; }

    inline double &getVolume(const unsigned int i) { return m_V; }
};

}  // namespace vox::flex