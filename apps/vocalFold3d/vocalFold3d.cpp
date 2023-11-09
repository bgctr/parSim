#include "olb3D.h"
#include "olb3D.hh"

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "spline.h"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::particles;
using namespace olb::particles::dynamics;

using T = double;

using DESCRIPTOR = D3Q19<OMEGA, POROSITY, VELOCITY_NUMERATOR, VELOCITY_DENOMINATOR>;

struct DYNBEHAVIOUR_CUSTOM : public PARTICLE_DESCRIPTOR<1, ACTIVE, DYNAMICS_ID, SCALAR>, public DYNBEHAVIOUR {};

//Define particleType
using PARTICLETYPE = PARTICLE_DESCRIPTOR<
  DESCRIPTOR::d, GENERAL_TMP<DESCRIPTOR::d>,
  MOBILITY_VERLET<DESCRIPTOR::d>, SURFACE_RESOLVED<DESCRIPTOR::d>,
  FORCING_RESOLVED<DESCRIPTOR::d>, PHYSPROPERTIES_RESOLVED<DESCRIPTOR::d>,
  DYNBEHAVIOUR_CUSTOM
>;

template <unsigned ID>
struct particles_of_dynamics_id {
  template<typename T, typename PARTICLETYPE>
  static bool value(Particle<T, PARTICLETYPE>& particle) {
    using namespace descriptors;
    return particle.template getField<DYNBEHAVIOUR, DYNAMICS_ID>() == ID;
  }
  static constexpr bool dynamic = true;
};

using BulkDynamics = PorousParticleBGKdynamics<T, DESCRIPTOR>;

const int N = 81;
const T maxPhysT = 0.5;

const T lengthX = 0.1;
const T lengthY = 0.01;
const T lengthZ = 0.015;

const T upperFoldY = 0.01 + 0.01 / N;
const T lowerFoldY = 0.00 - 0.01 / N;

std::vector<T> X_lower = {0,10,20,30,40};
std::vector<T> Y_lower = {0.00 - 0.01 / N, 0.00 - 0.01 / N, 0.00 - 0.01 / N, 0.00 - 0.01 / N, 0.00 - 0.01 / N};
std::vector<T> X_upper = {41, 51, 61, 71, 81};
std::vector<T> Y_upper = {0.01 + 0.01 / N, 0.01 + 0.01 / N, 0.01 + 0.01 / N, 0.01 + 0.01 / N, 0.01 + 0.01 / N};

template<typename T, typename PARTICLETYPE>
class LowerRestrictedVerletParticleDynamics : public ParticleDynamics<T, PARTICLETYPE> {
public:
  LowerRestrictedVerletParticleDynamics() {
    this->getName() = "RestrictedVerletParticleDynamics";
  }
  void process(Particle<T, PARTICLETYPE>& particle, T timeStepSize) override {
    using namespace descriptors;
    doWhenMeetingCondition<T, PARTICLETYPE, particles_of_dynamics_id<0>>(particle, [&]() {
      OstreamManager clout(std::cout, "ParticleId0");
	  T particleId = particle.getId();
	  //int rank;
      //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	  
	  //clout << "Rank Pos 1 = " << rank << std::endl;
	  //std::vector<T> X = {0, 10, 20, 30, 40};
	  //std::vector<T> Y(5);
	  auto acceleration = particles::access::getAcceleration(particle);
      auto position = particles::access::getPosition(particle);
	  if (particleId == 10 || particleId == 20 || particleId == 30){
		acceleration[0] = 0;
        T factor = particle.template getField<DYNBEHAVIOUR,SCALAR>();
        T displacement = - (position[1] - lowerFoldY) / 0.0089;
        acceleration[1] += util::max(0,util::log(100*displacement+1));
        acceleration[1] -= factor*util::max(0,-displacement);
        acceleration[2] = 0;
        particles::dynamics::velocityVerletTranslation(
          particle, timeStepSize, timeStepSize*timeStepSize, acceleration);
		if (particleId == 10){
			Y_lower[1] = position[1];
		}else if (particleId == 20){
			Y_lower[2] = position[1];
		}else{
			Y_lower[3] = position[1];
		}
	  }
	  /*
	  if (particleId == 0){
		  Y_lower[0] = position[1];
	  }
	  if (particleId == 40){
		  Y_lower[4] = position[1];
	  }
	  
	
	  tk::spline s(X_lower,Y_lower, tk::spline::cspline);
      for (T inter = 0; inter <= 40; inter++) { // Loop von 0 bis 40, inklusive
        T x = inter;
        T y = s(x);
        Vector<T,3> newPosition = {position[0], y, position[2]};
        particle.template setField<GENERAL, POSITION>(newPosition);
        //clout << "y = " << y << std::endl;		
      }
	  
	  //clout << "Rank Pos 3 = " << rank << std::endl;
	*/
	  
  });

}
};

template<typename T, typename PARTICLETYPE>
class UpperRestrictedVerletParticleDynamics : public ParticleDynamics<T,PARTICLETYPE> {
public:
  UpperRestrictedVerletParticleDynamics( ) {
    this->getName() = "RestrictedVerletParticleDynamics";
  }
  void process(Particle<T,PARTICLETYPE>& particle, T timeStepSize) override {
    using namespace descriptors;
    doWhenMeetingCondition<T,PARTICLETYPE,particles_of_dynamics_id<1>>(particle, [&](){
      OstreamManager clout( std::cout,"ParticleID = 1" );
		T particleId = particle.getId();
	  //int rank;
      //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	  
	  //clout << "Rank Pos 1 = " << rank << std::endl;
	  //std::vector<T> X = {0, 10, 20, 30, 40};
	  //std::vector<T> Y(5);
	  auto acceleration = particles::access::getAcceleration(particle);
      auto position = particles::access::getPosition(particle);
	  if (particleId == 51 || particleId == 61 || particleId == 71){
		acceleration[0] = 0;
        T factor = particle.template getField<DYNBEHAVIOUR,SCALAR>();
        T displacement = - (position[1] - lowerFoldY) / 0.0089;
        acceleration[1] += util::max(0,util::log(100*displacement+1));
        acceleration[1] -= factor*util::max(0,-displacement);
        acceleration[2] = 0;
        particles::dynamics::velocityVerletTranslation(
          particle, timeStepSize, timeStepSize*timeStepSize, acceleration);
		if (particleId == 51){
			Y_upper[1] = position[1];
		}else if (particleId == 61){
			Y_upper[2] = position[1];
		}else{
			Y_upper[3] = position[1];
		}
	  }
	  
	  /*
	  if (particleId == 41){
		  Y_upper[0] = position[1];
	  }
	  if (particleId == 81){
		  Y_upper[4] = position[1];
	  }
	 
	  tk::spline s(X_upper,Y_upper, tk::spline::cspline);
      for (T inter = 41; inter <= 81; inter++) { // Loop von 0 bis 40, inklusive
        T x = inter;
        T y = s(x);
        Vector<T,3> newPosition = {position[0], y, position[2]};
        particle.template setField<GENERAL, POSITION>(newPosition);
        //clout << "y = " << y << std::endl;		
      }
	  */
	  //clout << "Rank Pos 3 = " << rank << std::endl;

	  
    });
  }
};

/*

template<typename T, typename PARTICLETYPE>
class UpperRestrictedVerletParticleDynamics : public ParticleDynamics<T, PARTICLETYPE> {
public:
  UpperRestrictedVerletParticleDynamics() {
    this->getName() = "RestrictedVerletParticleDynamics";
  }
  void process(Particle<T, PARTICLETYPE>& particle, T timeStepSize) override {
    using namespace descriptors;
    doWhenMeetingCondition<T, PARTICLETYPE, particles_of_dynamics_id<1>>(particle, [&]() {
      OstreamManager clout(std::cout, "ParticleId1");

      T particleId = particle.getId();

      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   
      if (particleId == 0) {
        auto position = particles::access::getPosition(particle);
        slicePosInterpolsUpper[0] = position[1];
		
      }
      else if (particleId == 1) {
        auto position = particles::access::getPosition(particle);
        slicePosInterpolsUpper[1] = position[1];
		
      }

      if (particleId == 2 || particleId == 3 || particleId == 4 || particleId > 4) {
        auto acceleration = particles::access::getAcceleration(particle);
        auto position = particles::access::getPosition(particle);
	    
		//clout << particleId << std::endl;
		//clout << position << std::endl;
		
        acceleration[0] = 0;
        T factor = particle.template getField<DYNBEHAVIOUR, SCALAR>();
        T displacement = -(position[1] - upperFoldY) / 0.0089;
        acceleration[1] += util::max(0, util::log(100 * displacement + 1));
        acceleration[1] -= factor * util::max(0, -displacement);
        acceleration[2] = 0;
        particles::dynamics::velocityVerletTranslation(particle, timeStepSize, timeStepSize * timeStepSize, acceleration);

        if (particleId == 2) {
          slicePosInterpolsUpper[2] = position[1];
        }
        else if (particleId == 3) {
          slicePosInterpolsUpper[3] = position[1];
        }
        else if (particleId == 4) {
          slicePosInterpolsUpper[4] = position[1];
         
          if (rank == 0) {

            std::vector<T> X_upper = { 1, 2, 3, 4, 5 };
            T num_points_upper = 45;
            std::vector<T> Y_upper = { slicePosInterpolsUpper[0], slicePosInterpolsUpper[2], slicePosInterpolsUpper[3], slicePosInterpolsUpper[4], slicePosInterpolsUpper[1] };

            tk::spline s; //(X,Y,tk::spline::cspline_hermite)
            s.set_points(X_upper, Y_upper);
            T dx_upper = (X_upper.back() - X_upper.front()) / (num_points_upper - 1);

            for (T i = 0; i < num_points_upper; ++i) {
              T x_upper = X_upper.front() + i * dx_upper;
              interpolated_Y_upper[i] = s(x_upper);
            }
            interpolated_Y_upper.erase(interpolated_Y_upper.begin() + 44);
            interpolated_Y_upper.erase(interpolated_Y_upper.begin() + 33);
            interpolated_Y_upper.erase(interpolated_Y_upper.begin() + 22);
            interpolated_Y_upper.erase(interpolated_Y_upper.begin() + 11);
            interpolated_Y_upper.erase(interpolated_Y_upper.begin());
          }
        }
        else {
          auto position = particles::access::getPosition(particle);
          T interpolPosUpper = interpolated_Y_upper[particleId - 5];
          Vector<T,3> newPositionUpper = {position[0], interpolPosUpper, position[2]};
          particle.template setField<GENERAL, POSITION>(newPositionUpper);
        }
      }
  });

}
};

*/

void setParticlePosition(ParticleSystem<T, PARTICLETYPE>& particleSystem) {
  // Hier können Sie auf das particleSystem zugreifen und Änderungen vornehmen
  // Beispiel:
  OstreamManager clout(std::cout, "setParticlePosition");
  //auto& particles = particleSystem.getParticles();
  //auto position = particles::access::getPosition(particle);
  //int i = 0;
  //clout << position << std::endl;
  //auto particle = particleSystem.get(particleSystem.size());
  //clout << particle << std::endl;
	
	  for (std::size_t particleID = 0; particleID < particleSystem.size(); ++particleID) {
    Particle<T, PARTICLETYPE>& particle = particleSystem.get(particleID);
    Vector<T,3> currentPosition = particles::access::getPosition(particle);
    // Hier können Sie die Position aktualisieren, z.B.:
    Vector<T,3> modifiedPosition = {currentPosition[0], currentPosition[1], currentPosition[2]}; // Neue Position setzen
    particle.template setField<GENERAL,POSITION>(modifiedPosition);
  }
	
  }
}



void prepareGeometry(UnitConverter<T, DESCRIPTOR> const& converter,
  SuperGeometry<T, 3>& superGeometry)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  Vector<T, 3> extend(lengthX, lengthY, lengthZ);
  Vector<T, 3> origin;

  superGeometry.rename(0, 2);
  superGeometry.rename(2, 1, { 1,1,1 });

  // Set material number for inflow
  extend[0] = 2. * converter.getPhysDeltaX();
  origin[0] = -converter.getPhysDeltaX();
  IndicatorCuboid3D<T> inflow(extend, origin);
  superGeometry.rename(2, 3, 1, inflow);

  // Set material number for outflow
  origin[0] = lengthX - converter.getPhysDeltaX();
  IndicatorCuboid3D<T> outflow(extend, origin);
  superGeometry.rename(2, 4, 1, outflow);

  superGeometry.clean();
  superGeometry.checkForErrors();
  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice(SuperLattice<T, DESCRIPTOR>& sLattice,
  const UnitConverter<T, DESCRIPTOR>& converter,
  SuperGeometry<T, 3>& superGeometry)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=0 -->do nothing
  sLattice.defineDynamics<NoDynamics>(superGeometry, 0);

  sLattice.defineDynamics<BulkDynamics>(superGeometry, 1);
  sLattice.defineDynamics<BounceBack>(superGeometry, 2);

  setInterpolatedPressureBoundary(sLattice, converter.getLatticeRelaxationFrequency(), superGeometry, 3);
  setInterpolatedVelocityBoundary(sLattice, converter.getLatticeRelaxationFrequency(), superGeometry, 4);

  sLattice.setParameter<descriptors::OMEGA>(omega);

  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues(SuperLattice<T, DESCRIPTOR>& sLattice,
  const UnitConverter<T, DESCRIPTOR>& converter, std::size_t iT,
  SuperGeometry<T, 3>& superGeometry)
{
  OstreamManager clout(std::cout, "setBoundaryValues");

  // No of time steps for smooth start-up
  const std::size_t iTmaxStart = converter.getLatticeTime(0.0005);
  const std::size_t iTupdate = 10;

  T frac[1] = {};
  if (iT < iTmaxStart) {
    PolynomialStartScale<T, T> scale(iTmaxStart, T(1));
    T iTvec[1] = { T(iT) };
    scale(frac, iTvec);
  }
  else {
    frac[0] = 1;
  }

  if (iT % iTupdate == 0 && iT <= iTmaxStart) {
    const T targetVelocityX = converter.getCharLatticeVelocity() * frac[0];
    std::vector<T> targetU{ -targetVelocityX, 0., 0. };
    RectanglePoiseuille3D<T> poiseuilleU(superGeometry, 4, targetU, 0.5 * converter.getPhysDeltaX(), 0.5 * converter.getPhysDeltaX(), 0.5 * converter.getPhysDeltaX());
    sLattice.defineU(superGeometry, 4, poiseuilleU);
  }

  if (iT % converter.getLatticeTime(0.0001) == 0) {
    for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
      auto& porosity = sLattice.getBlock(iC).getField<POROSITY2>();
      for (std::size_t i = 0; i < sLattice.getBlock(iC).getNcells(); ++i) {
        porosity[0][i] = 0;
      }
    }
  }
}

void getResults(SuperLattice<T, DESCRIPTOR>& sLattice,
  UnitConverter<T, DESCRIPTOR> const& converter, int iT,
  SuperGeometry<T, 3>& superGeometry, util::Timer<T>& timer,
  ParticleSystem<T, PARTICLETYPE>& particleSystem)
{
  OstreamManager clout(std::cout, "getResults");

  SuperVTMwriter3D<T> vtkWriter("vocalFold3d");
  SuperLatticePhysVelocity3D velocity(sLattice, converter);
  SuperEuklidNorm3D<T> normVel(velocity);
  BlockReduction3D2D<T> planeReductionVel( normVel, {0, 0, 1} );
  BlockGifWriter<T> gifWriter;
  gifWriter.write( planeReductionVel, iT, "vel" );
  
  
  //SuperLatticePhysPressure3D pressure(sLattice, converter);
  SuperLatticeField3D<T, DESCRIPTOR, POROSITY2> porosity(sLattice);
  
  vtkWriter.addFunctor(velocity);
  // vtkWriter.addFunctor(pressure);
  // vtkWriter.addFunctor(porosity);
  // SuperEuklidNorm3D<T> normPoro(porosity);
  // BlockReduction3D2D<T> planeReductionPoro( normPoro, {0, 0, 1} );
  // BlockGifWriter<T> gifWriter;
  // write ppm image to file system
  // gifWriter.write( planeReductionPoro, iT, "poro" );

  if (iT == 0) {
    SuperLatticeGeometry3D geometry(sLattice, superGeometry);
    SuperLatticeCuboid3D cuboid(sLattice);
    SuperLatticeRank3D rank(sLattice);
    vtkWriter.write(geometry);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);
    vtkWriter.createMasterFile();
  }

	// clout << converter.getLatticeTime(0.005) << std::endl;

  if (iT % converter.getLatticeTime(0.00005) == 0) {
    vtkWriter.write(iT);
  }

  //if (iT % converter.getLatticeTime(0.0001) == 0) {
  timer.update(iT);
  timer.printStep();
  sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
  /*{
    auto particle = particleSystem.get(particleSystem.size() / 4);
    io::printResolvedParticleInfo(particle);
  }*/
  /*{
    auto particle = particleSystem.get(3 * (particleSystem.size() / 4));
    io::printResolvedParticleInfo(particle);
  }*/
  }
//}

int main(int argc, char* argv[])
{
  // === 1st Step: Initialization ===
  olbInit(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");
  OstreamManager clout(std::cout, "main");

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
    int{ N },      // resolution: number of voxels per charPhysL
    (T)0.51,   // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)0.01,   // charPhysLength: reference length of simulation geometry
    (T)5,      // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)0.0002, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)1.0     // physDensity: physical density in __kg / m^3__
  );
  converter.print();

  // === 2nd Step: Prepare Geometry ===
  Vector<T, 3> origin;
  Vector<T, 3> extend(lengthX, lengthY, lengthZ);
  IndicatorCuboid3D<T> cuboid(extend, origin);
  CuboidGeometry3D<T> cuboidGeometry(cuboid, converter.getPhysDeltaX(), singleton::mpi().getSize());
  HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

  SuperGeometry<T, 3> superGeometry(cuboidGeometry, loadBalancer, 2);
  prepareGeometry(converter, superGeometry);

  SuperLattice<T, DESCRIPTOR> sLattice(superGeometry);
  prepareLattice(sLattice, converter, superGeometry);

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  std::shared_ptr<ParticleDynamics<T,PARTICLETYPE>> lowerParticleDynamics(new LowerRestrictedVerletParticleDynamics<T, PARTICLETYPE>());
  std::shared_ptr<ParticleDynamics<T,PARTICLETYPE>> upperParticleDynamics(new UpperRestrictedVerletParticleDynamics<T, PARTICLETYPE>());
  particleSystem.addDynamics(lowerParticleDynamics);
  particleSystem.addDynamics(upperParticleDynamics);

  ParticleManager<T, DESCRIPTOR, PARTICLETYPE> particleManager(particleSystem, superGeometry, sLattice, converter, Vector<T, 3>(0.));

  std::shared_ptr<IndicatorF3D<T>> sliceBottomI(new STLreader<T>("stimmlippe.stl", converter.getPhysDeltaX(), 1));
  for (int iZ=0; iZ < util::ceil((lengthZ+2*converter.getPhysDeltaX()) / (0.0005-2*converter.getPhysDeltaX())); ++iZ) {
    Vector<T,3> slicePos{lengthX-0.03,
                         lowerFoldY,
                         -1*converter.getPhysDeltaX() + iZ*(0.0005-converter.getPhysDeltaX())};
    if (slicePos[2] < lengthZ+converter.getPhysDeltaX()) {
      creators::addResolvedArbitraryShape3D(particleSystem, slicePos, converter.getPhysDeltaX(), sliceBottomI, T{0.5}*converter.getPhysDeltaX(), T{1});
      particleSystem.get(particleSystem.size()-1).template setField<DYNBEHAVIOUR,DYNAMICS_ID>(0); // lower
      particleSystem.get(particleSystem.size()-1).template setField<DYNBEHAVIOUR,SCALAR>(
        0.0001*util::abs((slicePos[2] - 0.5*lengthZ) / (0.5*lengthZ)));
    }
  }
  std::shared_ptr<IndicatorF3D<T>> sliceTopI(new STLreader<T>("stimmlippe_top.stl", converter.getPhysDeltaX(), 1));
  for (int iZ=0; iZ < util::ceil((lengthZ+2*converter.getPhysDeltaX()) / (0.0005-2*converter.getPhysDeltaX())); ++iZ) {
    Vector<T,3> slicePos{lengthX-0.03,
                         upperFoldY,
                         -1*converter.getPhysDeltaX() + iZ*(0.0005-converter.getPhysDeltaX())};
    if (slicePos[2] < lengthZ+converter.getPhysDeltaX()) {
      creators::addResolvedArbitraryShape3D(particleSystem, slicePos, converter.getPhysDeltaX(), sliceTopI, T{0.5}*converter.getPhysDeltaX(), T{1});
      particleSystem.get(particleSystem.size()-1).template setField<DYNBEHAVIOUR,DYNAMICS_ID>(1); // upper
      particleSystem.get(particleSystem.size()-1).template setField<DYNBEHAVIOUR,SCALAR>(
        0.0001*util::abs((slicePos[2] - 0.5*lengthZ) / (0.5*lengthZ)));
    }
  }
 
  // === 4th Step: Main Loop with Timer ===
  util::Timer<T> timer(converter.getLatticeTime(maxPhysT), superGeometry.getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT = 0; iT < converter.getLatticeTime(maxPhysT); ++iT) {
    setBoundaryValues(sLattice, converter, iT, superGeometry);

    particleManager.execute<
      couple_lattice_to_particles<T, DESCRIPTOR, PARTICLETYPE>,
      process_dynamics<T, PARTICLETYPE>,
      couple_particles_to_lattice<T, DESCRIPTOR, PARTICLETYPE>
    >();

	// Berechnung des spline
	
	// interpolation
	
	setParticlePosition(particleSystem);

    sLattice.collideAndStream();

    getResults(sLattice, converter, iT, superGeometry, timer, particleSystem);
  }

  timer.stop();
  timer.printSummary();
}
