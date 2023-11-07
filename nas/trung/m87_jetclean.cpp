//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cluster_setup.cpp
//! \brief Problem generator for HD galaxy cluster based on M87. By Trung Ha (UNT)
//!

/*
Virgo cluster setup, using similar initial condition as Guo, Stone, Kim, & Quataert, 2023

This simplified Problem Generator includes: perturbation, gravity, radiative cooling, cooling-based AMR only, Gaussian open jet profile with accretion-dependent power only, inner/outer spherical boundaries, cold accretion

For Virgo:
haloMass = 3.e13 M_sun
haloRadius = 6.e-2 Mpc
stellarMass = 3e11 M_sun
stellarRadius = 2e-3 Mpc
scaledEntropy = 7 keV cm^2
etaPower = 1.1
scaledRadius = 2e-3 Mpc
*/

// C++ headers
#include <algorithm> // min, max
#include <cmath>     // log
#include <cstring>   // strcmp()
#include <iostream>
#include <vector>
#include <random>
#include <cstdint>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../defs.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../hydro/hydro_diffusion/hydro_diffusion.hpp" // diffusion
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"
#include "../utils/utils.hpp"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// User-defined functions, explanations available at each function implementation

void innerX1Boundary(MeshBlock *pmb, Coordinates *pco,
                     AthenaArray<Real> &prim, FaceField &b,
                     Real time, Real dt,
                     int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void outerX1Boundary(MeshBlock *pmb, Coordinates *pco,
                     AthenaArray<Real> &prim, FaceField &b,
                     Real time, Real dt,
                     int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void innerX2Boundary(MeshBlock *pmb, Coordinates *pco,
                     AthenaArray<Real> &prim, FaceField &b,
                     Real time, Real dt,
                     int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void outerX2Boundary(MeshBlock *pmb, Coordinates *pco,
                     AthenaArray<Real> &prim, FaceField &b,
                     Real time, Real dt,
                     int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void innerX3Boundary(MeshBlock *pmb, Coordinates *pco,
                     AthenaArray<Real> &prim, FaceField &b,
                     Real time, Real dt,
                     int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void outerX3Boundary(MeshBlock *pmb, Coordinates *pco,
                     AthenaArray<Real> &prim, FaceField &b,
                     Real time, Real dt,
                     int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void radialBoundaries(MeshBlock *pmb, const AthenaArray<Real> &prim,
                      AthenaArray<Real> &cons, int k, int j, int i,
                      Real z, Real y, Real x, Real r);

int RefinementCondition(MeshBlock *pmb);

void allSourceFunctions(MeshBlock *pmb, const Real time, const Real dt,
                        const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
                        const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                        AthenaArray<Real> &cons_scalar);

void gravitySourceFunction(MeshBlock *pmb, const Real dt,
                           const AthenaArray<Real> &prim, AthenaArray<Real> &cons,
                           int k, int j, int i, Real z, Real y, Real x, Real r);
static Real massFuncOfR(const Real distanceFromSMBH);
static Real gravitationalAcceleration(const Real distanceFromSMBH);

static Real coreEntropy(const Real dimensionlessX);
static Real dN_dX(const Real dimensionlessX, const Real numberDensity);
static Real rungeKutta4(Real x0, Real y0, Real xn, int n, Real (*differentialFunction)(Real, Real));
static Real numberDensityFromEntropy(const Real dimensionlessX);
static Real pressureFromEntropyAndNumberDensity(const Real entropy, const Real numberDensity);
static void tabulatedNumberDensitiesAndPressures();

static void crossProduct(std::vector<Real> r, std::vector<Real> p, Real *x, Real *y, Real *z);
static Real interpolate1D(Real inputX, std::vector<Real> tabulatedX, std::vector<Real> tabulatedY);
void enforceFloors(MeshBlock *pmb, AthenaArray<Real> &cons, int k, int j, int i);
static Real coolingTime(AthenaArray<Real> &w, int k, int j, int i);
static Real soundCrossingTime(AthenaArray<Real> &w, int k, int j, int i, Real cellWidth);

void coolingSourceFunction(const Real dt, const AthenaArray<Real> &prim,
                           AthenaArray<Real> &cons, int k, int j, int i);
static Real emissivityFromTemperature(Real temperature);
static Real coolingFunction(Real numberDensity, Real emissivity);

void jetFeedbackSourceFunction(MeshBlock *pmb, AthenaArray<Real> &cons, const AthenaArray<Real> &prim,
                               int k, int j, int i, const Real dt,
                               Real z, Real y, Real x, Real polarDistance);

Real coldAccretionRate(MeshBlock *pmb, int iout);
Real hotAccretionRate(MeshBlock *pmb, int iout);
Real jetPower(MeshBlock *pmb, int iout);
Real angularMomentumX(MeshBlock *pmb, int iout);
Real angularMomentumY(MeshBlock *pmb, int iout);
Real angularMomentumZ(MeshBlock *pmb, int iout);

// Global variables put in unnamed namespace to avoid linkage issues
namespace
{
    // Global constants
    static int numberOfRefinementLevels;
    static Real mu = 0.5924489; // average mass per particle, for n_H = rho / mu * m_H equation, modified on 04/22/2023 to be consistent with yt implementation assuming fully ionized gas

    // Global variables
    Real gammaAdiabatic, gammaMinus1;
    Real innerRadius; // inner radius, used for smoothing to avoid NaN at the center
    Real outerRadius; // outer radius, used to force regions outside to hydrostatic equilibrium
    Real numberDensityFloorCGS, temperatureFloor, pressureFloorCGS, numberDensityVacuumSinkCGS;
    Real numberDensityFloorAstronomical, pressureFloorAstronomical, numberDensityVacuumSinkAstronomical, pressureVacuumSinkAstronomical, densityVacuumSinkAstronomical, densityFloorAstronomical;
    Real haloMass, haloRadius;                  // dark-matter-halo-related parameters
    Real stellarMass, stellarRadius;            // BCG-related parameters
    Real scaledEntropy, etaPower, scaledRadius; // entropy-related parameters
    Real SMBHMass;
    Real scaledEntropyCGS, scaledEntropyAstronomical;
    Real initialNumberDensityCGS, initialNumberDensityAstronomical; // initial condition at x=1, used to solve hydrostatic equilibrium equation numerically
    Real densityAtOuterBoundary, pressureAtOuterBoundary;
    Real jetLaunchingHeightScale, jetLaunchingHeight, jetFeedbackEfficiency, jetKineticFraction, jetConversionEfficiency, jetLaunchingWidthScale, jetLaunchingWidth;
    Real jetVelocityKms, jetVelocityAstronomical, jetPrecessionAngleTheta, jetPrecessionPeriod;
    Real coolingStartTime, jetStartTime;
    Real accretionUpdateFrequency, accretionMaxTemperature; // update accretion rate (aka jet power) every n Myr, set maximum temperature of gas that is counted toward accretion
    Real amrTimeOn;

    Real simulationBoxWidth, numberOfMeshBlocks, numberOfZones, smallestCellWidth;
    bool outputJetSpecial = false;
    bool outputAMRSpecial = false;
    bool restartFlag;

    int turbulenceType, seed;
    Real densityPerturbation, velocityPerturbationKms, velocityPerturbationAstronomical, velocityPerturbationWavelength;

    std::string configuration, additionalComment;

    // Create a table of radial distances, number densities, and pressure, to then interpolate to find density and pressure at a given r
    const int numberOfTabulatedEntries = 1000; // number of points in the tabulated data table, 1000 should be fine
    std::vector<Real> dimensionlessXArray(numberOfTabulatedEntries);
    std::vector<Real> entropyFuncOfXArray(numberOfTabulatedEntries);
    std::vector<Real> numberDensityFuncOfXArray(numberOfTabulatedEntries);
    std::vector<Real> densityFuncOfXArray(numberOfTabulatedEntries);
    std::vector<Real> pressureFuncOfXArray(numberOfTabulatedEntries);

    // Declare variables to store cooling curve from Schure+09, Schneider & Robertson+18, and Koyama & Inutsuka+02 (solar metallicity)
    std::vector<Real> logTemperatureArray{3., 3.03, 3.06, 3.09, 3.121, 3.151, 3.181, 3.211, 3.241,
                                          3.271, 3.302, 3.332, 3.362, 3.392, 3.422, 3.452, 3.482, 3.513,
                                          3.543, 3.573, 3.603, 3.633, 3.663, 3.693, 3.724, 3.754, 3.784,
                                          3.814, 3.844, 3.874, 3.905, 3.935, 3.965, 3.995, 4.025, 4.055,
                                          4.085, 4.116, 4.146, 4.176, 4.206, 4.236, 4.266, 4.296, 4.327,
                                          4.357, 4.387, 4.417, 4.447, 4.477, 4.508, 4.538, 4.568, 4.598,
                                          4.628, 4.658, 4.688, 4.719, 4.749, 4.779, 4.809, 4.839, 4.869,
                                          4.899, 4.93, 4.96, 4.99, 5.02, 5.05, 5.08, 5.111, 5.141,
                                          5.171, 5.201, 5.231, 5.261, 5.291, 5.322, 5.352, 5.382, 5.412,
                                          5.442, 5.472, 5.503, 5.533, 5.563, 5.593, 5.623, 5.653, 5.683,
                                          5.714, 5.744, 5.774, 5.804, 5.834, 5.864, 5.894, 5.925, 5.955,
                                          5.985, 6.015, 6.045, 6.075, 6.106, 6.136, 6.166, 6.196, 6.226,
                                          6.256, 6.286, 6.317, 6.347, 6.377, 6.407, 6.437, 6.467, 6.497,
                                          6.528, 6.558, 6.588, 6.618, 6.648, 6.678, 6.709, 6.739, 6.769,
                                          6.799, 6.829, 6.859, 6.889, 6.92, 6.95, 6.98, 7.01, 7.04,
                                          7.07, 7.101, 7.131, 7.161, 7.191, 7.221, 7.251, 7.281, 7.312,
                                          7.342, 7.372, 7.402, 7.432, 7.462, 7.492, 7.523, 7.553, 7.583,
                                          7.613, 7.643, 7.673, 7.704, 7.734, 7.764, 7.794, 7.824, 7.854,
                                          7.884, 7.915, 7.945, 7.975, 8.005, 8.035, 8.065, 8.095, 8.126,
                                          8.156, 8.186, 8.216, 8.246, 8.276, 8.307, 8.337, 8.367, 8.397,
                                          8.427, 8.457, 8.487, 8.518, 8.548, 8.578, 8.608, 8.638, 8.668,
                                          8.698, 8.729, 8.759, 8.789, 8.819, 8.849, 8.879, 8.91, 8.94,
                                          8.97, 9.};
    std::vector<Real> logEmissivityHydroArray{-26.0928, -26.075, -26.0575, -26.0401, -26.0228, -26.0057,
                                              -25.9887, -25.9719, -25.9552, -25.9386, -25.922, -25.9056,
                                              -25.8893, -25.8731, -25.8569, -25.8408, -25.8248, -25.8088,
                                              -25.7929, -25.777, -25.761, -25.7446, -25.7267, -25.7045,
                                              -25.6703, -25.6076, -25.489, -25.2909, -25.0182, -24.7013,
                                              -24.3696, -24.0409, -23.7232, -23.4198, -23.1318, -22.8592,
                                              -22.6017, -22.3587, -22.1298, -21.9141, -21.5767, -21.4745,
                                              -21.473, -21.5081, -21.5615, -21.6154, -21.647, -21.6523,
                                              -21.6292, -21.5736, -21.5057, -21.4322, -21.3547, -21.2746,
                                              -21.1955, -21.1193, -21.0472, -20.9809, -20.9209, -20.8653,
                                              -20.815, -20.7712, -20.735, -20.706, -20.6852, -20.674,
                                              -20.6796, -20.6934, -20.7097, -20.7229, -20.7213, -20.713,
                                              -20.7014, -20.6894, -20.6819, -20.6771, -20.6738, -20.6711,
                                              -20.674, -20.6935, -20.7373, -20.812, -20.9338, -21.068,
                                              -21.1941, -21.2991, -21.362, -21.3976, -21.4183, -21.4312,
                                              -21.4498, -21.4845, -21.5293, -21.5796, -21.6218, -21.6492,
                                              -21.667, -21.678, -21.687, -21.7003, -21.716, -21.7328,
                                              -21.7469, -21.7565, -21.7644, -21.7727, -21.7859, -21.8116,
                                              -21.85, -21.9014, -21.9664, -22.0362, -22.1037, -22.1665,
                                              -22.2214, -22.2645, -22.2992, -22.3261, -22.3448, -22.3549,
                                              -22.3603, -22.3615, -22.3591, -22.3534, -22.3469, -22.3403,
                                              -22.3344, -22.332, -22.3328, -22.3369, -22.3444, -22.3557,
                                              -22.3687, -22.3837, -22.4008, -22.4221, -22.4461, -22.4724,
                                              -22.5002, -22.5272, -22.5514, -22.5725, -22.5901, -22.6024,
                                              -22.6114, -22.6175, -22.6208, -22.6212, -22.6197, -22.6166,
                                              -22.6121, -22.606, -22.5988, -22.5908, -22.5821, -22.5728,
                                              -22.5628, -22.5526, -22.542, -22.5309, -22.5198, -22.5087,
                                              -22.4976, -22.4864, -22.4752, -22.4638, -22.4524, -22.441,
                                              -22.4294, -22.4179, -22.4065, -22.3949, -22.3813, -22.3678,
                                              -22.3542, -22.3406, -22.3271, -22.3135, -22.2999, -22.2864,
                                              -22.2728, -22.2592, -22.2457, -22.2321, -22.2185, -22.2049,
                                              -22.1914, -22.1778, -22.1642, -22.1507, -22.1371, -22.1235,
                                              -22.11, -22.0964, -22.0828, -22.0693, -22.0557, -22.0421,
                                              -22.0286, -22.015};

    // Code conversion, assuming codeTemperature = 1
    static Real codeLength;
    static Real codeMass;
    static Real codeVolume;
    static Real codeTemperature;
    static Real codeArea;
    static Real codeTime;
    static Real codeVelocity;
    static Real codeMomentum;
    static Real codeMomentumDensity;
    static Real codeDensity;
    static Real codeNumberDensity;
    static Real codeAcceleration;
    static Real codeForce;
    static Real codePressure;
    static Real codeEnergy;
    static Real codePower;
    static Real codeEnergyDensity;
    static Real codeGravitationalConst;
    static Real codeBoltzmannConst;

    // Some units to cgs
    const Real kpcCGS = 3.08567758096e+21;
    const Real MpcCGS = kpcCGS * 1.e3;
    const Real kmCGS = 1.e5;
    const Real gravitationalConstCGS = 6.67384e-08;
    const Real protonMassCGS = 1.67373522381e-24;
    const Real solarMassCGS = 1.988e33;
    const Real boltzmannConstCGS = 1.3806488e-16;
    const Real yearCGS = 31557600.0;
    const Real MyrCGS = yearCGS * 1.e6;
    const Real eV_CGS = 1.60217733e-12;
    const Real hydrogenMassCGS = 1.6735575e-24; // hydrogen atom mass in gram
    const Real cCGS = 2.9979e10;

    // cgs to astronomical units (Msun, Mpc, Myr, K)
    const Real boltzmannConstAstronomical = boltzmannConstCGS / (solarMassCGS * pow(MpcCGS, 2) * pow(MyrCGS, -2));
    const Real eV_Astronomical = eV_CGS / (solarMassCGS * SQR(MpcCGS) * pow(MyrCGS, -2.));
    const Real hydrogenMassAstronomical = hydrogenMassCGS / solarMassCGS;                                                              // hydrogen atom mass in M_sun (8.418e-64)
    const Real gravitationalConstAstronomical = gravitationalConstCGS / (pow(MpcCGS, 3.) * pow(MyrCGS, -2.) * pow(solarMassCGS, -1.)); // gravitational constant in Mpc^3 Myr^-2 M_sun^-1

    // Additional unit conversion
    const Real kms_Astronomical = kmCGS / (MpcCGS / MyrCGS); // convert velocity from km s^-1 to Mpc Myr^-1

} // namespace

//----------------------------------------------------------------------------------------
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.

void Mesh::InitUserMeshData(ParameterInput *pin)
{
    // Read Problem Parameters
    gammaAdiabatic = pin->GetReal("hydro", "gamma"); // adiabatic gamma factor
    gammaMinus1 = gammaAdiabatic - 1;

    simulationBoxWidth = mesh_size.x1max - mesh_size.x1min;
    numberOfMeshBlocks = pin->GetInteger("meshblock", "nx1");
    numberOfZones = pin->GetInteger("mesh", "nx1");
    numberOfRefinementLevels = pin->GetInteger("mesh", "numlevel");
    smallestCellWidth = simulationBoxWidth / numberOfZones / pow(2., numberOfRefinementLevels - 1);

    codeLength = pin->GetReal("conversion", "box_width") / simulationBoxWidth; // set unit code length based on the width of the cubic box in Mpc
    codeMass = pin->GetReal("conversion", "base_mass");                        // unit code mass in M_sun
    codeTime = pin->GetReal("conversion", "base_time");                        // unit code time in Myr
    codeTemperature = pin->GetOrAddReal("conversion", "base_temp", 1.);        // unit code temperature in Kelvin, default 1

    numberDensityFloorCGS = pin->GetOrAddReal("problem", "n_floor", 1.e-3); // enforce number density floor in cm^-3
    numberDensityFloorAstronomical = numberDensityFloorCGS / pow(MpcCGS, -3);
    densityFloorAstronomical = numberDensityFloorAstronomical * mu * hydrogenMassAstronomical;
    numberDensityVacuumSinkCGS = pin->GetOrAddReal("problem", "n_sink", 5.e-3); // enforce number density inside central sink region in cm^-3
    numberDensityVacuumSinkAstronomical = numberDensityVacuumSinkCGS / pow(MpcCGS, -3);
    densityVacuumSinkAstronomical = numberDensityVacuumSinkAstronomical * mu * hydrogenMassAstronomical;
    temperatureFloor = pin->GetOrAddReal("problem", "t_floor", 2.e5);                                           // enforce temperature floor in K
    pressureFloorAstronomical = numberDensityFloorAstronomical * boltzmannConstAstronomical * temperatureFloor; // enforce pressure floor, inferred from P = nkT, in astronomical units
    pressureVacuumSinkAstronomical = numberDensityVacuumSinkAstronomical * boltzmannConstAstronomical * temperatureFloor;

    haloMass = pin->GetOrAddReal("problem", "M_DM", 3.e13);    // dark matter halo mass in M_sun
    haloRadius = pin->GetOrAddReal("problem", "r0_DM", 6.e-2); // halo scale radius in Mpc

    stellarMass = pin->GetOrAddReal("problem", "M_star", 3.e11); // BCG stellar mass in M_sun
    stellarRadius = pin->GetOrAddReal("problem", "r0_s", 2.e-3); // BCG stellar radius in Mpc
    SMBHMass = pin->GetOrAddReal("problem", "M_SMBH", 6.5e9);    // SMBH mass in M_sun

    scaledEntropy = pin->GetOrAddReal("problem", "K0", 7.); // central entropy in keV cm^2
    scaledEntropyCGS = scaledEntropy * eV_CGS * 1.e3;
    scaledEntropyAstronomical = scaledEntropyCGS / (solarMassCGS * pow(MpcCGS, 4) * pow(MyrCGS, -2));
    etaPower = pin->GetOrAddReal("problem", "eta", 1.1);               // power slope
    scaledRadius = pin->GetOrAddReal("problem", "r0_e", 2.e-3);        // scaled radius of the entropy profile in Mpc
    initialNumberDensityCGS = pin->GetOrAddReal("problem", "n0", 0.1); // initial condition of number density in cm^-3
    initialNumberDensityAstronomical = initialNumberDensityCGS / pow(MpcCGS, -3);

    innerRadius = pin->GetOrAddReal("boundary", "r_in", 3e-5);  // inner (smoothing) radius in Mpc
    outerRadius = pin->GetOrAddReal("boundary", "r_out", 0.18); // by default, set outer radius to 180 kpc, or slightly inside the box

    accretionMaxTemperature = pin->GetOrAddReal("accretion", "max_temp", 1.e6); // set a maximum temperature that is considered to be accreted by the SMBH

    coolingStartTime = pin->GetOrAddReal("cooling", "t_cool", 5.); // only start cooling after this much time has passed, in Myr

    jetStartTime = pin->GetOrAddReal("jet", "t_jet", 5.);            // only start jet feedback after this much time has passed, in Myr.
    jetKineticFraction = pin->GetOrAddReal("jet", "f_kin", 0.5);     // which fraction of the jet power is injected as kinetic vs thermal energy, from Martizzi+2019
    jetLaunchingHeightScale = pin->GetOrAddReal("jet", "h_jet", 6.); // the scale height from which the jet launching layer is set, default to six times innerRadius
    jetLaunchingHeight = innerRadius * (jetLaunchingHeightScale);
    jetLaunchingWidthScale = pin->GetOrAddReal("jet", "w_jet", 0.5); // the scale width of the jet launching platform, wrt the innerRadius
    jetLaunchingWidth = innerRadius * jetLaunchingWidthScale;
    jetVelocityKms = pin->GetOrAddReal("jet", "v_jet", 9.e3); // initial velocity of the jet in km/s
    jetVelocityAstronomical = jetVelocityKms * kms_Astronomical;
    jetConversionEfficiency = 0.5 * SQR(jetVelocityAstronomical * (MpcCGS / MyrCGS) / cCGS) / jetKineticFraction; // the conversion efficiency of accretion rate to jet power, epsilon = v^2 / 2fc^2

    accretionUpdateFrequency = pin->GetOrAddReal("jet", "m_upd", 5.); // update accretion rate (aka jet power) every n Myr
    if (jetStartTime < accretionUpdateFrequency)
    {
        jetStartTime = accretionUpdateFrequency; // if non-static jet, start jet only after galaxy has accumulated mass for n years
        outputJetSpecial = true;
    }

    jetPrecessionAngleTheta = pin->GetOrAddReal("jet", "a_precess", 15.); // jet precession angle in degrees
    jetPrecessionAngleTheta *= PI / 180.;                                 // convert to radians
    jetPrecessionPeriod = pin->GetOrAddReal("jet", "t_precess", 10.);     // frequency of changing precession angle in Myr

    seed = 1;                                                    // some random seeding for random generator, doesn't matter what the actual value is
    turbulenceType = pin->GetInteger("turbulence", "turb_type"); // 0 for density perturbation, 1 for velocity perturbation
    if (turbulenceType == 0)
    {
        densityPerturbation = pin->GetOrAddReal("turbulence", "turb_amp", 5.e-2); // deviation of density from hydrostatic equilibrium solution as a fraction
    }
    else
    {
        velocityPerturbationKms = pin->GetOrAddReal("turbulence", "turb_vel", 100.); // dispersion of velocity perturbation in km s^-1
        velocityPerturbationAstronomical = velocityPerturbationKms * kms_Astronomical;
        // velocityPerturbationWavelength = pin->GetOrAddReal("problem", "lambda_pert", 8.e-3); // wavelength of velocity perturbation
    }
    restartFlag = (pin->GetInteger("job", "res_flag") == 1 ? true : false); // check whether job is new or restart

    configuration = pin->GetString("comment", "configure");
    additionalComment = pin->GetOrAddString("comment", "comment", ""); // add any necessary comments about the simulation setup

    if (adaptive)
    {
        amrTimeOn = pin->GetOrAddReal("amr", "time_on", 40.); // time to turn on condition-specific AMR, might save some time
    }

    // Conversions to code units-----------------------------------------------------------------
    codeVolume = pow(codeLength, 3);
    codeArea = pow(codeLength, 2);
    codeVelocity = codeLength / codeTime;
    codeMomentum = codeVelocity * codeMass;
    codeDensity = codeMass / codeVolume;
    codeNumberDensity = 1. / codeVolume;
    codeMomentumDensity = codeVelocity * codeDensity; // momentum density, Athena++ defines momentum per cell as the base IM1, IM2, IM3
    codeAcceleration = codeVelocity / codeTime;
    codeForce = codeMass * codeAcceleration;
    codePressure = codeForce / codeArea;
    codeEnergy = codeForce * codeLength;
    codePower = codeEnergy / codeTime;
    codeEnergyDensity = codeDensity * codeLength * codeAcceleration;                                            // total energy density, Athena++ defines energy per cell as the base IEN
    codeGravitationalConst = gravitationalConstAstronomical / (pow(codeLength, 3.) / SQR(codeTime) / codeMass); // convert gravitational constant to code unit
    codeBoltzmannConst = boltzmannConstAstronomical / (codeEnergy / codeTemperature);                           // convert Boltzmann constant to code unit

    // Initial conditions------------------------------------------------------------------------
    tabulatedNumberDensitiesAndPressures(); // generate a tabulated dataset of initial densities and pressures as a function of x
    densityFuncOfXArray = numberDensityFuncOfXArray;
    for (auto &it : densityFuncOfXArray)
        it *= mu * hydrogenMassAstronomical;

    // Boundary Conditions-----------------------------------------------------------------------
    densityAtOuterBoundary = interpolate1D(outerRadius / scaledRadius, dimensionlessXArray, densityFuncOfXArray); // the density and pressure outside the outerRadius is fixed
    pressureAtOuterBoundary = interpolate1D(outerRadius / scaledRadius, dimensionlessXArray, pressureFuncOfXArray);
    if (mesh_bcs[BoundaryFace::inner_x1] == GetBoundaryFlag("user"))
    {
        EnrollUserBoundaryFunction(BoundaryFace::inner_x1, innerX1Boundary);
    }
    if (mesh_bcs[BoundaryFace::outer_x1] == GetBoundaryFlag("user"))
    {
        EnrollUserBoundaryFunction(BoundaryFace::outer_x1, outerX1Boundary);
    }
    if (mesh_bcs[BoundaryFace::inner_x2] == GetBoundaryFlag("user"))
    {
        EnrollUserBoundaryFunction(BoundaryFace::inner_x2, innerX2Boundary);
    }
    if (mesh_bcs[BoundaryFace::outer_x2] == GetBoundaryFlag("user"))
    {
        EnrollUserBoundaryFunction(BoundaryFace::outer_x2, outerX2Boundary);
    }
    if (mesh_bcs[BoundaryFace::inner_x3] == GetBoundaryFlag("user"))
    {
        EnrollUserBoundaryFunction(BoundaryFace::inner_x3, innerX3Boundary);
    }
    if (mesh_bcs[BoundaryFace::outer_x3] == GetBoundaryFlag("user"))
    {
        EnrollUserBoundaryFunction(BoundaryFace::outer_x3, outerX3Boundary);
    }

    // Enroll Adaptive Mesh Refinement conditions------------------------------------------------
    if (adaptive)
    {
        EnrollUserRefinementCondition(RefinementCondition);
    }

    // Enroll source function--------------------------------------------------------------------
    EnrollUserExplicitSourceFunction(allSourceFunctions); // all source functions (gravity, radial boundary, cooling, jet) are inside this function

    // Print data to output file-----------------------------------------------------------------
    if (Globals::my_rank == 0)
    {
        std::cout.precision(3);

        if (restartFlag == true)
        {
            std::cout << "----RESTART JOB----"
                      << "\n";
            std::cout << "\n";
        }

        std::cout << "----SIMULATION PARAMETERS----"
                  << "\n";
        std::cout << "Simulation box width = " << simulationBoxWidth << "\n";
        std::cout << "Number of zones = " << numberOfZones << "\n";
        std::cout << "Number of cell per mesh block dimension = " << numberOfMeshBlocks << "\n";
        std::cout << "Time between hst dumps = " << pin->GetReal("output1", "dt") * codeTime << " Myr \n";
        std::cout << "Time between hdf5 dumps = " << pin->GetReal("output2", "dt") * codeTime << " Myr \n";
        std::cout << "Simulation start time = " << time * codeTime << " Myr \n";
        std::cout << "Final simulation time = " << pin->GetReal("time", "tlim") * codeTime << " Myr \n";
        std::cout << "CFL number = " << pin->GetReal("time", "cfl_number") << "\n";
        std::cout << "Time integrator = " << pin->GetString("time", "integrator") << "\n";
        std::cout << "Number of SMR levels = " << numberOfRefinementLevels - 1 << "\n";
        std::cout << "\n";

        if (adaptive)
        {
            std::cout << "----ADAPTIVE MESH REFINEMENT----"
                      << "\n";
            // std::cout << "AMR is ON, up to an additional " << maxAdditionalAMRLevels << " refinement levels \n";
            std::cout << "AMR based on cooling / sound crossing condition \n";
            if (amrTimeOn > 0.)
            {
                std::cout << "AMR is delayed until " << amrTimeOn << " Myr \n";
                outputAMRSpecial = true;
            }
            std::cout << "\n";
        }

        std::cout << "----UNIT CONVERSIONS----"
                  << "\n";
        std::cout << "Code length = " << codeLength << " Mpc \n";
        std::cout << "Code mass = " << codeMass << " Msun \n";
        std::cout << "Code time = " << codeTime << " Myr \n";
        std::cout << "Code temperature = " << codeTemperature << " Kelvin \n";
        std::cout << "\n";

        std::cout << "---PHYSICAL PARAMETERS----"
                  << "\n";
        std::cout << "Scaled entropy = " << scaledEntropy << " keV cm^2 \n";
        std::cout << "Scaled radius = " << scaledRadius * 1.e3 << " kpc \n";
        std::cout << "\n";

        std::cout << "Halo mass = " << haloMass << " Msun \n";
        std::cout << "Stellar mass = " << stellarMass << " Msun \n";
        std::cout << "SMBH mass = " << SMBHMass << " Msun \n";
        std::cout << "\n";

        std::cout << "Inner radius = " << innerRadius * 1.e6 << " pc \n";
        std::cout << "Dark matter radius = " << haloRadius * 1.e3 << " kpc \n";
        std::cout << "Stellar radius = " << stellarRadius * 1.e3 << " kpc \n";
        std::cout << "\n";

        std::cout << "----PHYSICS----"
                  << "\n";
        std::cout << "Stellar potential is ON"
                  << "\n";
        std::cout << "Inner vacuum sink region is ON"
                  << "\n";
        std::cout << "Outer radial boundary condition is ON"
                  << "\n";
        std::cout << "Cooling is ON"
                  << "\n";

        std::cout << "Jet feedback is ON"
                  << "\n";
        if (turbulenceType == 0)
        {
            std::cout << "Turbulent density perturbation is ON"
                      << "\n";
        }
        else
        {
            std::cout << "Turbulent velocity perturbation is ON"
                      << "\n";
        }
        std::cout << "Hot accretion is OFF"
                  << "\n";
        std::cout << "\n";

        std::cout << "----BOUNDARY CONDITIONS----"
                  << "\n";
        std::cout << "Outer box boundary condition along X is: " << pin->GetString("mesh", "ix1_bc") << "\n";
        std::cout << "Outer box boundary condition along Y is: " << pin->GetString("mesh", "ix2_bc") << "\n";
        std::cout << "Outer box boundary condition along Z is: " << pin->GetString("mesh", "ix3_bc") << "\n";

        std::cout << "Outer radius = " << outerRadius * 1.e3 << " kpc \n";
        std::cout << "Number density at boundary = " << densityAtOuterBoundary / (mu * hydrogenMassAstronomical) * pow(MpcCGS, -3) << " cm^-3 \n";
        std::cout << "Pressure at boundary = " << pressureAtOuterBoundary * (solarMassCGS / MpcCGS / SQR(MyrCGS)) / 10 << " Pa \n";

        if (outerRadius > codeLength)
        {
            std::cout << "WARNING: outer radius is outside the box, might return nonphysical values"
                      << "\n";
        }
        std::cout << "\n";

        std::cout << "----TURBULENT PERTURBATION----"
                  << "\n";
        if (turbulenceType == 0)
        {
            std::cout << "Amplitude as a fraction = " << densityPerturbation << "\n";
        }
        else
        {
            std::cout << "Velocity perturbation dispersion = " << velocityPerturbationKms << " km s^-1 \n";
            // std::cout << "Velocity perturbation wavelength = " << velocityPerturbationWavelength * 1.e3 << " kpc \n";
        }
        std::cout << "\n";

        std::cout << "----ACCRETION----"
                  << "\n";
        std::cout << "Maximum accretion temperature = " << accretionMaxTemperature << " K \n";
        std::cout << "\n";

        std::cout << "----RADIATIVE COOLING----"
                  << "\n";
        std::cout << "Cooling start time = " << coolingStartTime << " Myr \n";
        std::cout << "\n";

        std::cout << "----AGN JET FEEDBACK----"
                  << "\n";
        std::cout << "Jet start time = " << jetStartTime << " Myr \n";
        if (outputJetSpecial)
        {
            std::cout << "Jet start time delay enforced. Set t_jet >= m_upd \n";
            outputJetSpecial = false;
        }
        std::cout << "Jet launch height = " << jetLaunchingHeight * 1.e6 << " pc \n";
        std::cout << "Jet base radius = " << jetLaunchingWidth * 1.e6 << " pc \n";
        std::cout << "Jet velocity = " << jetVelocityKms << " km s^-1 \n";
        std::cout << "Jet conversion efficiency = " << jetConversionEfficiency << "\n";
        std::cout << "Jet kinetic fraction = " << jetKineticFraction << "\n";
        std::cout << "M_dot updates frequency = " << accretionUpdateFrequency << " Myr \n";
        std::cout << "Jet profile is Gaussian \n";
        std::cout << "Jet precession is ON"
                  << "\n";
        std::cout << "Precession polar angle = " << jetPrecessionAngleTheta * 180. / PI << " degrees \n";
        std::cout << "Period of precession = " << jetPrecessionPeriod << " Myr \n";
        std::cout << "\n";

        std::cout << "----FLOOR VALUES----"
                  << "\n";
        std::cout << "Number density floor = " << numberDensityFloorCGS << " cm^-3 \n";
        std::cout << "Number density inside sink region = " << numberDensityVacuumSinkCGS << " cm^-3 \n";
        std::cout << "Temperature floor = " << temperatureFloor << " K \n";
        // std::cout << "Pressure floor = " << pressureFloorAstronomical * (solarMassCGS / MpcCGS / SQR(MyrCGS)) / 10 << " Pa \n";
        std::cout << "\n";

        std::cout << "----DIAGNOSTIC PARAMETERS----"
                  << "\n";
        std::cout << "Width of the smallest cell = " << smallestCellWidth * codeLength * 1.e6 << " pc \n";
        // std::cout << "Array of X starts at " << dimensionlessXArray[0] * scaledRadius * 1.e6 << " pc \n";
        // std::cout << "Array of X ends at " << dimensionlessXArray[numberOfTabulatedEntries - 1] * scaledRadius * 1.e3 << " kpc \n";
        std::cout << "Expected mass of the inner region = " << pow(innerRadius, 3.) * densityVacuumSinkAstronomical * PI * 4. / 3. << " Msun \n";
        // std::cout << "Root refinement level = " << root_level << "\n";
        std::cout << "\n";

        std::cout << "----COMMENTS----"
                  << "\n";
        std::cout << "To configure, do: python configure.py " << configuration << "\n";
        if (!additionalComment.empty())
        {
            std::cout << additionalComment << "\n";
        }
        std::cout << "\n";

        std::cout << "----SIMULATION STARTS----"
                  << "\n";
    }

    // Declare global variables that need frequent updates
    AllocateRealUserMeshDataField(8);
    ruser_mesh_data[0].NewAthenaArray(1);
    ruser_mesh_data[0](0) = 0.; // base mass inside inner region
    ruser_mesh_data[1].NewAthenaArray(2);
    ruser_mesh_data[1](0) = 0.;           // cold gas
    ruser_mesh_data[1](1) = 0.;           // hot gas
    ruser_mesh_data[2].NewAthenaArray(2); // mass accreted in the last n time step (corresponding to the accretionUpdateFrequency), used to load onto jet
    ruser_mesh_data[2](0) = 0.;
    ruser_mesh_data[2](1) = 0.;           // hot mass accreted in the last n time steps
    ruser_mesh_data[3].NewAthenaArray(2); // total mass accreted
    ruser_mesh_data[3](0) = 0.;
    ruser_mesh_data[3](1) = 0.;           // total hot mass accreted
    ruser_mesh_data[4].NewAthenaArray(1); // a mesh data object that keeps track of the last time jet power was updated
    ruser_mesh_data[4](0) = time;
    ruser_mesh_data[5].NewAthenaArray(2); // keep track of current <M_dot> in code units
    ruser_mesh_data[5](0) = 0.;
    ruser_mesh_data[5](1) = 0.;           // keep track of current <M_dot_hot> in code units
    ruser_mesh_data[6].NewAthenaArray(1); // keep track of current jet power in code units
    ruser_mesh_data[6](0) = 0.;
    ruser_mesh_data[7].NewAthenaArray(3); // keep track of the cell height at jet base and precession angle, for consistency (also other jet values in the future, modify as appropriate)
    ruser_mesh_data[7](0) = 0.;           // jet base height
    ruser_mesh_data[7](1) = 0.;           // current precession angle, THIS SEEMS TO BECOME 279 EVERY TIME RUN IS RESTARTED, SEED NEEDS TO BE SAVED ACROSS RUNS, OOPS
    ruser_mesh_data[7](2) = time;         // last time precession angle was updated

    // Declare seed for consistency
    AllocateIntUserMeshDataField(1);
    iuser_mesh_data[0].NewAthenaArray(1);
    iuser_mesh_data[0](0) = seed;

    // Declare history output variables
    AllocateUserHistoryOutput(3);
    EnrollUserHistoryOutput(0, coldAccretionRate, "mdot_c", UserHistoryOperation::max);
    EnrollUserHistoryOutput(1, hotAccretionRate, "mdot_h", UserHistoryOperation::max);
    EnrollUserHistoryOutput(2, jetPower, "p_jet", UserHistoryOperation::max);
    // EnrollUserHistoryOutput(3, angularMomentumX, "L_x", UserHistoryOperation::sum);
    // EnrollUserHistoryOutput(4, angularMomentumY, "L_y", UserHistoryOperation::sum);
    // EnrollUserHistoryOutput(5, angularMomentumZ, "L_z", UserHistoryOperation::sum);

    return;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the Virgo cluster

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
    // Prepare index bounds
    int il = is - NGHOST;
    int iu = ie + NGHOST;
    int jl = js;
    int ju = je;
    if (block_size.nx2 > 1)
    {
        jl -= (NGHOST);
        ju += (NGHOST);
    }
    int kl = ks;
    int ku = ke;
    if (block_size.nx3 > 1)
    {
        kl -= (NGHOST);
        ku += (NGHOST);
    }

    // Read problem parameters
    for (int k = kl; k <= ku; ++k)
    {
        Real z = pcoord->x3v(k);
        for (int j = jl; j <= ju; ++j)
        {
            Real y = pcoord->x2v(j);
            for (int i = il; i <= iu; ++i)
            {
                Real x = pcoord->x1v(i);
                Real r = sqrt(SQR(x) + SQR(y) + SQR(z)); // radial distance in code unit

                Real rAstronomical = r * codeLength;                // radial distance in Mpc
                Real dimensionlessX = rAstronomical / scaledRadius; // a dimensionless x parameter, equals to r/r_0

                Real interpolatedDensity = interpolate1D(dimensionlessX, dimensionlessXArray, densityFuncOfXArray);
                Real interpolatedPressure = interpolate1D(dimensionlessX, dimensionlessXArray, pressureFuncOfXArray);

                if (r <= innerRadius / codeLength) // inner region
                {
                    phydro->u(IDN, k, j, i) = densityVacuumSinkAstronomical / codeDensity;
                    phydro->u(IEN, k, j, i) = pressureVacuumSinkAstronomical / codePressure / gammaMinus1;
                    phydro->u(IM1, k, j, i) = 0.;
                    phydro->u(IM2, k, j, i) = 0.;
                    phydro->u(IM3, k, j, i) = 0.;
                }
                else
                {
                    phydro->u(IDN, k, j, i) = interpolatedDensity / codeDensity; // u(IDN, k, j, i) is the conserved quantity, w(IDN, k, j, i) is the primitive quantity
                    phydro->u(IEN, k, j, i) = interpolatedPressure / codePressure / gammaMinus1;
                    phydro->u(IM1, k, j, i) = 0.0;
                    phydro->u(IM2, k, j, i) = 0.0; // initial motion of the gas field
                    phydro->u(IM3, k, j, i) = 0.0;

                    // Initial turbulence field
                    if (turbulenceType == 0)
                    {
                        // generate random density perturbation
                        auto randomGenerator = std::mt19937(++(pmy_mesh->iuser_mesh_data[0](0)));
                        auto normalDistribution = std::normal_distribution<Real>{0.0, densityPerturbation / 2};
                        Real randomDensityPerturbation = normalDistribution(randomGenerator);

                        phydro->u(IDN, k, j, i) *= 1 + randomDensityPerturbation;
                    }
                    else
                    { // THIS IS BROKEN AGAIN, COMES OUT TO ZERO PERTURBATION
                        // generate random velocity perturbation with given standard deviation and smoothing wavelength
                        auto randomGenerator = std::mt19937(++(pmy_mesh->iuser_mesh_data[0](0)));
                        auto normalDistribution = std::normal_distribution<Real>{0.0, velocityPerturbationAstronomical / 2};
                        Real randomVelocityPerturbationAstronomical = normalDistribution(randomGenerator); // normal distribution, centered around 0, with a 1 sigma of velocity perturbation

                        // Added 04/03/2023 to try a built-in random number generator instead
                        int64_t iseed = -1 - ++(pmy_mesh->iuser_mesh_data[0](0));

                        Real randomPhiRad = 2. * PI * ran2(&iseed);
                        Real randomThetaRad = PI * ran2(&iseed - 1);

                        if (pmy_mesh->iuser_mesh_data[0](0) % 1000000 == 0 && Globals::my_rank == 0)
                        {
                            std::cout << "seed = " << pmy_mesh->iuser_mesh_data[0](0) << ", velocity perturbation = " << randomVelocityPerturbationAstronomical / kms_Astronomical << " km s^-1, phi = " << randomPhiRad << ", theta = " << randomThetaRad << "\n";
                        }

                        Real velocityPerturbationX = randomVelocityPerturbationAstronomical * sin(randomThetaRad) * cos(randomPhiRad);
                        Real velocityPerturbationY = randomVelocityPerturbationAstronomical * sin(randomThetaRad) * sin(randomPhiRad);
                        Real velocityPerturbationZ = randomVelocityPerturbationAstronomical * cos(randomThetaRad);

                        phydro->u(IM1, k, j, i) += phydro->u(IDN, k, j, i) * (velocityPerturbationX / codeVelocity);
                        phydro->u(IM2, k, j, i) += phydro->u(IDN, k, j, i) * (velocityPerturbationY / codeVelocity);
                        phydro->u(IM3, k, j, i) += phydro->u(IDN, k, j, i) * (velocityPerturbationZ / codeVelocity);
                    }
                }
            }
        }
    }
    return;
}

//----------------------------------------------------------------------------------------
//  \brief User-defined functions

//----------------------------------------------------------------------------------------
//! \fn int RefinementCondition(MeshBlock *pmb)
//  \brief Specify adaptive mesh refinement condition, added 07/08/2023

int RefinementCondition(MeshBlock *pmb)
{
    // Mimicking SMR setup but for radial spacing
    int currentRefinementLevel, requiredRefinementLevel;
    int requiredSMRLevel = numberOfRefinementLevels; // at first, assuming the refinement needed is highest
    int requiredAMRLevel = 0;
    int maximumAllowedAMRLevel = numberOfRefinementLevels - 5; // highest level of AMR allowed w/r to the total number of refinement levels; reduce this number to speed up simulation
    int maximumAllowedDifferenceSMR_AMRLevel = 2;              // highest allowable difference in AMR vs. the base SMR at this meshblock
    int refinementCountFlag = 0;
    Real cellLengthCode = 0.;
    Real rootCellLengthCode = simulationBoxWidth / numberOfZones; // the width of the biggest cell in the domain
    Real maximumRadialDistanceCode = simulationBoxWidth / 2.;

    AthenaArray<Real> &w = pmb->phydro->w;

    // Average location of meshblock
    Real meshblockX1 = (pmb->block_size.x1max + pmb->block_size.x1min) / 2.;
    Real meshblockX2 = (pmb->block_size.x2max + pmb->block_size.x2min) / 2.;
    Real meshblockX3 = (pmb->block_size.x3max + pmb->block_size.x3min) / 2.;
    // Real meshblockR = sqrt(SQR(meshblockX1) + SQR(meshblockX2) + SQR(meshblockX3));
    Real maximumMeshBlockX = fmax(fmax(fabs(meshblockX1), fabs(meshblockX2)), fabs(meshblockX3));
    requiredSMRLevel = fmin(requiredSMRLevel, int(log2(maximumRadialDistanceCode / maximumMeshBlockX))); // every time we go twice in radial distance, refinement level goes down by one

    // Added 08/02/2023: a flag implementation to exit the nested for loop in case certain conditions are met
    int nestedLoopFlag = 0;
    for (int k = pmb->ks - 1; k <= pmb->ke + 1; k++)
    {
        for (int j = pmb->js - 1; j <= pmb->je + 1; j++)
        {
            for (int i = pmb->is - 1; i <= pmb->ie + 1; i++)
            {
                if (refinementCountFlag == 0)
                {
                    cellLengthCode = pmb->pcoord->GetEdge1Length(k, j, i);
                    currentRefinementLevel = round(log2(rootCellLengthCode / cellLengthCode));
                    refinementCountFlag = 1;
                }

                if (requiredSMRLevel >= maximumAllowedAMRLevel)
                {
                    nestedLoopFlag = 1;
                    break; // if the SMR level is already high enough, ignore AMR
                }

                if (pmb->pmy_mesh->time < amrTimeOn / codeTime)
                {
                    nestedLoopFlag = 1;
                    break; // before AMR is turned on, only one step of the loop to record current refinement level is needed
                }
                else
                {
                    Real coolingBySoundCrossingRatio = coolingTime(w, k, j, i) / soundCrossingTime(w, k, j, i, cellLengthCode);

                    if (coolingBySoundCrossingRatio <= 6.)
                    { // check condition for AMR
                        requiredAMRLevel = fmin(currentRefinementLevel + 1, maximumAllowedAMRLevel);
                        nestedLoopFlag = 1;
                        break; // if the condition has been met, skip through the nested for loop to save computational time
                    }
                }
            }

            if (nestedLoopFlag == 1)
            {
                break;
            }
        }

        if (nestedLoopFlag == 1)
        {
            break;
        }
    }

    // Combine SMR and AMR condition
    if (requiredAMRLevel > requiredSMRLevel + maximumAllowedDifferenceSMR_AMRLevel)
    {
        requiredRefinementLevel = requiredSMRLevel + maximumAllowedDifferenceSMR_AMRLevel; // if the required AMR is more than two levels above SMR at that meshblock, only refine by 2 at most (avoid runaway refinement that exceeds memory allocation)
    }
    else
    {
        requiredRefinementLevel = fmax(requiredAMRLevel, requiredSMRLevel); // take the maximum between the two requirements otherwise
    }

    if (requiredRefinementLevel > currentRefinementLevel)
        return 1;
    if (requiredRefinementLevel < currentRefinementLevel)
        return -1;
    return 0;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh:UserWorkInLoop()
//  \brief Analysis of accretion rate after each cycle

void Mesh::UserWorkInLoop()
{
    Real accretionPerTimeStepCode;      // the difference in the base mass before and after reset is the accreted mass
    Real innerMassCode = 0.;            // base mass of the inner region
    Real updatedInnerColdMassCode = 0.; // cold mass after accretion
    Real updatedInnerHotMassCode = 0.;  // hot mass after accretion

    Real localInnerMass, localUpdatedColdMass, localUpdatedHotMass; // value per MPI rank, to be summed over all ranks below
    // only need localInnerMass, since there should be no hot gas after update

    if (ncycle == 0)
    {
        ruser_mesh_data[1](0) = ruser_mesh_data[0](0); // do this to save a repeated calculation on each cell at cycle 0
    }

    localInnerMass = ruser_mesh_data[0](0);
    localUpdatedColdMass = ruser_mesh_data[1](0);
    localUpdatedHotMass = ruser_mesh_data[1](1);

#ifdef MPI_PARALLEL
    MPI_Allreduce(&localInnerMass, &innerMassCode, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD); // sync data from local MPI rank to the rest by summing all values

    MPI_Allreduce(&localUpdatedColdMass, &updatedInnerColdMassCode, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&localUpdatedHotMass, &updatedInnerHotMassCode, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif

    if (ncycle == 0)
    {
        if (Globals::my_rank == 0)
        {
            std::cout << "Total mass of the inner region = " << innerMassCode * codeMass << " Msun \n";
            std::cout << "\n";
        }
    }
    else
    {
        accretionPerTimeStepCode = fmax(updatedInnerColdMassCode - innerMassCode, 0.); // in the case where no cold gas is accreted, the actual updated cold mass might be less than the expected cold mass, so set it to zero to prevent negative jet power

        ruser_mesh_data[2](0) += accretionPerTimeStepCode; // tally up accretion until next accretion update
        ruser_mesh_data[3](0) += accretionPerTimeStepCode; // tally up accretion for the entire simulation

        ruser_mesh_data[2](1) += updatedInnerHotMassCode; // tally up hot accretion for analysis purpose
        ruser_mesh_data[3](1) += updatedInnerHotMassCode; // tally up hot accretion for the entire simulation

        if (time >= ruser_mesh_data[4](0) + accretionUpdateFrequency / codeTime)
        {                                                                                          // update the amount of mass loaded onto jet per accretionUpdateFrequency
            ruser_mesh_data[5](0) = ruser_mesh_data[2](0) / (accretionUpdateFrequency / codeTime); // <M_dot> = M_total / time, in code units
            ruser_mesh_data[5](1) = ruser_mesh_data[2](1) / (accretionUpdateFrequency / codeTime); // <M_dot_hot>

            Real jetPowerCGS = jetConversionEfficiency * (ruser_mesh_data[5](0) * codeMass / codeTime) * (solarMassCGS / MyrCGS) * pow(cCGS, 2); // P_jet = eff * <M_dot> * c^2
            ruser_mesh_data[6](0) = jetPowerCGS / (solarMassCGS * pow(MpcCGS, 2) * pow(MyrCGS, -3)) / codePower;                                 // update jet power in code units, save for restart

            if (Globals::my_rank == 0)
            {
                std::cout << "Time = " << time * codeTime << " Myr, <M_dot> = " << ruser_mesh_data[5](0) * codeMass / codeTime / 1.e6 << " Msun yr^-1\n";
                std::cout << "Jet power = " << jetPowerCGS << " erg s^-1 \n";
            }

            ruser_mesh_data[2](0) = 0.;                                   // reset mass accreted per unit update frequency
            ruser_mesh_data[2](1) = 0.;                                   // reset hot mass accreted per unit update frequency
            ruser_mesh_data[4](0) += accretionUpdateFrequency / codeTime; // update time tracker
        }

        if (time > ruser_mesh_data[7](2) + jetPrecessionPeriod / codeTime)
        { // precession
            int64_t iseed = -1 - ++(iuser_mesh_data[0](0));
            ruser_mesh_data[7](1) = ran2(&iseed) * 2. * PI; // generate a random angle between 0 and 2 pi

            if (Globals::my_rank == 0)
            {
                std::cout << "Jet azimuthal angle = " << ruser_mesh_data[7](1) * 180. / PI << " degrees \n";
            }

            ruser_mesh_data[7](2) += jetPrecessionPeriod / codeTime; // update time tracker
        }

        if (ncycle % ncycle_out == 0 && Globals::my_rank == 0)
        {
            std::cout.precision(4);
            // std::cout << "Mass of the inner region = " << innerMassCode * codeMass << " Msun \n";
            // every ncycle, output the instantaneous accretion rate
            std::cout << "Hot accretion rate = " << updatedInnerHotMassCode * codeMass / (dt * codeTime * 1.e6) << " Msun yr^-1 \n";
            std::cout << "Cold accretion rate = " << (updatedInnerColdMassCode - innerMassCode) * codeMass / (dt * codeTime * 1.e6) << " Msun yr^-1 \n";
            // std::cout << "\n";

            if (time > amrTimeOn / codeTime && outputAMRSpecial)
            {
                std::cout << "AMR is turned on from now \n";
                outputAMRSpecial = false;
            }
        }
    }
    ruser_mesh_data[1](0) = 0.;
    ruser_mesh_data[1](1) = 0.; // reset mass accreted for next time step
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh:UserWorkAfterLoop(ParameterInput *pin)
//  \brief Print out the total mass accreted at the end of the simulation
void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{
#ifdef MPI_PARALLEL
    MPI_Allreduce(MPI_IN_PLACE, &(ruser_mesh_data[7](0)), 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD); // sync data from local MPI rank to the rest by taking the maximum value
    MPI_Allreduce(MPI_IN_PLACE, &(ruser_mesh_data[7](1)), 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
#endif

    if (Globals::my_rank == 0)
    {
        std::cout.precision(4);
        std::cout << "\n";
        std::cout << "Current jet power = " << ruser_mesh_data[6](0) * codePower * (solarMassCGS * pow(MpcCGS, 2) * pow(MyrCGS, -3)) << " erg s^-1 \n";
        std::cout << "Current <M_dot> = " << ruser_mesh_data[5](0) * codeMass / codeTime / 1.e6 << " Msun yr^-1 \n";
        std::cout << "Jet base cell height = " << ruser_mesh_data[7](0) * codeLength * 1.e6 << " pc \n";
        std::cout << "Jet precession azimuthal angle = " << ruser_mesh_data[7](1) * 180. / PI << " degrees \n";
        std::cout << "Total accreted mass = " << ruser_mesh_data[3](0) * codeMass << " Msun \n";
    }
}

// BOUNDARY CONDITIONS------------------------------------------------------------------------------------

/*
Force initial hydrostatic equilibrium BCs
*/

void innerX1Boundary(MeshBlock *pmb, Coordinates *pco,
                     AthenaArray<Real> &prim, FaceField &b,
                     Real time, Real dt,
                     int il, int iu, int jl, int ju, int kl, int ku, int ngh)
{
    for (int k = kl; k <= ku; ++k)
    {
        for (int j = jl; j <= ju; ++j)
        {
            for (int i = 1; i <= ngh; ++i)
            {
                prim(IDN, k, j, il - i) = densityAtOuterBoundary / codeDensity;
                prim(IPR, k, j, il - i) = pressureAtOuterBoundary / codePressure;
                prim(IVX, k, j, il - i) = 0.;
                prim(IVY, k, j, il - i) = 0.;
                prim(IVZ, k, j, il - i) = 0.;
            }
        }
    }
}

void outerX1Boundary(MeshBlock *pmb, Coordinates *pco,
                     AthenaArray<Real> &prim, FaceField &b,
                     Real time, Real dt,
                     int il, int iu, int jl, int ju, int kl, int ku, int ngh)
{
    for (int k = kl; k <= ku; ++k)
    {
        for (int j = jl; j <= ju; ++j)
        {
            for (int i = 1; i <= ngh; ++i)
            {
                prim(IDN, k, j, iu + i) = densityAtOuterBoundary / codeDensity;
                prim(IPR, k, j, iu + i) = pressureAtOuterBoundary / codePressure;
                prim(IVX, k, j, iu + i) = 0.;
                prim(IVY, k, j, iu + i) = 0.;
                prim(IVZ, k, j, iu + i) = 0.;
            }
        }
    }
}

void innerX2Boundary(MeshBlock *pmb, Coordinates *pco,
                     AthenaArray<Real> &prim, FaceField &b,
                     Real time, Real dt,
                     int il, int iu, int jl, int ju, int kl, int ku, int ngh)
{
    for (int k = kl; k <= ku; ++k)
    {
        for (int j = 1; j <= ngh; ++j)
        {
            for (int i = il; i <= iu; ++i)
            {
                prim(IDN, k, jl - j, i) = densityAtOuterBoundary / codeDensity;
                prim(IPR, k, jl - j, i) = pressureAtOuterBoundary / codePressure;
                prim(IVX, k, jl - j, i) = 0.;
                prim(IVY, k, jl - j, i) = 0.;
                prim(IVZ, k, jl - j, i) = 0.;
            }
        }
    }
}

void outerX2Boundary(MeshBlock *pmb, Coordinates *pco,
                     AthenaArray<Real> &prim, FaceField &b,
                     Real time, Real dt,
                     int il, int iu, int jl, int ju, int kl, int ku, int ngh)
{
    for (int k = kl; k <= ku; ++k)
    {
        for (int j = 1; j <= ngh; ++j)
        {
            for (int i = il; i <= iu; ++i)
            {
                prim(IDN, k, ju + j, i) = densityAtOuterBoundary / codeDensity;
                prim(IPR, k, ju + j, i) = pressureAtOuterBoundary / codePressure;
                prim(IVX, k, ju + j, i) = 0.;
                prim(IVY, k, ju + j, i) = 0.;
                prim(IVZ, k, ju + j, i) = 0.;
            }
        }
    }
}

void innerX3Boundary(MeshBlock *pmb, Coordinates *pco,
                     AthenaArray<Real> &prim, FaceField &b,
                     Real time, Real dt,
                     int il, int iu, int jl, int ju, int kl, int ku, int ngh)
{
    for (int k = 1; k <= ngh; ++k)
    {
        for (int j = jl; j <= iu; ++j)
        {
            for (int i = il; i <= iu; ++i)
            {
                prim(IDN, kl - k, j, i) = densityAtOuterBoundary / codeDensity;   // force density = density | r = box_width / 2
                prim(IPR, kl - k, j, i) = pressureAtOuterBoundary / codePressure; // force pressure = pressure | r = box_width / 2
                prim(IVX, kl - k, j, i) = 0.;                                     // force zero velocity
                prim(IVY, kl - k, j, i) = 0.;
                prim(IVZ, kl - k, j, i) = 0.;
            }
        }
    }
}

void outerX3Boundary(MeshBlock *pmb, Coordinates *pco,
                     AthenaArray<Real> &prim, FaceField &b,
                     Real time, Real dt,
                     int il, int iu, int jl, int ju, int kl, int ku, int ngh)
{
    for (int k = 1; k <= ngh; ++k)
    {
        for (int j = jl; j <= iu; ++j)
        {
            for (int i = il; i <= iu; ++i)
            {
                prim(IDN, ku + k, j, i) = densityAtOuterBoundary / codeDensity;   // force density = density | r = box_width / 2
                prim(IPR, ku + k, j, i) = pressureAtOuterBoundary / codePressure; // force pressure = pressure | r = box_width / 2
                prim(IVX, ku + k, j, i) = 0.;                                     // force zero velocity
                prim(IVY, ku + k, j, i) = 0.;
                prim(IVZ, ku + k, j, i) = 0.;
            }
        }
    }
}

/*
Radial boundaries
*/

void radialBoundaries(MeshBlock *pmb, const AthenaArray<Real> &prim,
                      AthenaArray<Real> &cons, int k, int j, int i,
                      Real z, Real y, Real x, Real r)
{
    Real primDensity = prim(IDN, k, j, i);
    // Real valueW = pmb->phydro->w(IDN, k, j, i);
    Real valueW1 = pmb->phydro->w1(IDN, k, j, i);

    if (r <= innerRadius / codeLength)
    {
        if (valueW1 == primDensity) // check to see if the integrator has moved to the second sub-step of the rk2. Only record accretion here to avoid data being counted twice
        {                           // BAD FIX, CONSIDER WRITING MESHBLOCK::USERWORKINLOOP OR MEMORY ADDRESS POINTER TO AVOID CASES WHERE W = W1

            // Added 06/23/2023: calculate current temperature of cell to determine if accretion is hot or cold
            Real conservedNumberDensityCode = cons(IDN, k, j, i) / (mu * hydrogenMassAstronomical / codeMass);
            Real kineticEnergyDensityCode = 0.5 * (SQR(cons(IM1, k, j, i)) + SQR(cons(IM2, k, j, i)) + SQR(cons(IM3, k, j, i))) / cons(IDN, k, j, i); // unit: mass / length / time^2 (checked)
            Real thermalEnergyDensityCode = cons(IEN, k, j, i) - kineticEnergyDensityCode;
            Real conservedTemperatureCode = gammaMinus1 * thermalEnergyDensityCode / (conservedNumberDensityCode * codeBoltzmannConst);

            Real cellVolumeCode = pmb->pcoord->GetCellVolume(k, j, i);

            if (pmb->pmy_mesh->ncycle == 0)
            { // at the start, record the base mass; at every time step after, record the modified mass
                pmb->pmy_mesh->ruser_mesh_data[0](0) += primDensity * cellVolumeCode;
            }

            if (conservedTemperatureCode <= accretionMaxTemperature / codeTemperature)
            {
                pmb->pmy_mesh->ruser_mesh_data[1](0) += cons(IDN, k, j, i) * cellVolumeCode; // add up the mass contribution inside the inner region before density is reset to floor
            }
            else
            {
                pmb->pmy_mesh->ruser_mesh_data[1](1) += cons(IDN, k, j, i) * cellVolumeCode; // add up the mass contribution inside the inner region before density is reset to floor
            }
        }

        cons(IDN, k, j, i) = densityVacuumSinkAstronomical / codeDensity;                 // force sink density
        cons(IEN, k, j, i) = pressureVacuumSinkAstronomical / codePressure / gammaMinus1; // force floor pressure based on floor temperature
        cons(IM1, k, j, i) = 0.;                                                          // force zero velocity
        cons(IM2, k, j, i) = 0.;
        cons(IM3, k, j, i) = 0.;
    }
    else if (r > outerRadius / codeLength) // force outer region to be one value
    {
        cons(IDN, k, j, i) = densityAtOuterBoundary / codeDensity;                 // force density = density | r = outer radius
        cons(IEN, k, j, i) = pressureAtOuterBoundary / codePressure / gammaMinus1; // force pressure = pressure | r = outer radius
        cons(IM1, k, j, i) = 0.;                                                   // force zero velocity
        cons(IM2, k, j, i) = 0.;
        cons(IM3, k, j, i) = 0.;
    }
}

// END BOUNDARY CONDITIONS------------------------------------------------------------------------

// SOURCE FUNCTION--------------------------------------------------------------------------------
void allSourceFunctions(MeshBlock *pmb, const Real time, const Real dt,
                        const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
                        const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                        AthenaArray<Real> &cons_scalar)
{
    for (int k = pmb->ks; k <= pmb->ke; ++k)
    {
        Real z = pmb->pcoord->x3v(k);
        for (int j = pmb->js; j <= pmb->je; ++j)
        {
            Real y = pmb->pcoord->x2v(j);
            for (int i = pmb->is; i <= pmb->ie; ++i)
            {
                Real x = pmb->pcoord->x1v(i);
                Real r = sqrt(SQR(x) + SQR(y) + SQR(z));    // radial distance in code unit
                Real polarDistance = sqrt(SQR(x) + SQR(y)); // polar distance in code unit

                // GRAVITY
                gravitySourceFunction(pmb, dt, prim, cons, k, j, i, z, y, x, r);
                // END GRAVITY

                // COOLING
                if (time > coolingStartTime / codeTime && r > innerRadius / codeLength && r < outerRadius / codeLength) // apply cooling only for region between the inner and outer radius
                {
                    coolingSourceFunction(dt, prim, cons, k, j, i);
                }
                // END COOLING

                // ENFORCE TEMPERATURE AND DENSITY FLOOR
                enforceFloors(pmb, cons, k, j, i); // 06/17/2023: moving this up here, since only cooling would result in a temperature crash
                // END ENFORCE TEMPERATURE AND DENSITY FLOOR

                // JET FEEDBACK
                if (time > jetStartTime / codeTime && polarDistance <= jetLaunchingWidth / codeLength && r < 2 * jetLaunchingHeight / codeLength) // apply jet feedback only to the regions close to the jet platform
                {
                    jetFeedbackSourceFunction(pmb, cons, prim, k, j, i, dt, z, y, x, polarDistance);
                }
                // END JET FEEDBACK

                // RADIAL BOUNDARIES
                radialBoundaries(pmb, prim, cons, k, j, i, z, y, x, r);
                // END RADIAL BOUNDARIES
            }
        }
    }
}

// END SOURCE FUNCTION------------------------------------------------------------------------------------

// GRAVITY------------------------------------------------------------------------------------------------

/*
Gravity function, part of the source function
*/

void gravitySourceFunction(MeshBlock *pmb, const Real dt, const AthenaArray<Real> &prim,
                           AthenaArray<Real> &cons, int k, int j, int i, Real z, Real y, Real x, Real r)
{
    Real rDoubleDot = gravitationalAcceleration(r * codeLength) / codeAcceleration; // radial acceleration in code units

    Real deltaU = dt * rDoubleDot * x / r; // sin(theta) * cos(phi) = x/r
    Real deltaV = dt * rDoubleDot * y / r; // sin(theta) * sin(phi) = y/r
    Real deltaW = dt * rDoubleDot * z / r; // cos(theta) = z/r

    cons(IM1, k, j, i) += prim(IDN, k, j, i) * deltaU;
    cons(IM2, k, j, i) += prim(IDN, k, j, i) * deltaV;
    cons(IM3, k, j, i) += prim(IDN, k, j, i) * deltaW;
    if (NON_BAROTROPIC_EOS)
        cons(IEN, k, j, i) += dt * prim(IDN, k, j, i) * rDoubleDot / r * (prim(IVX, k, j, i) * x + prim(IVY, k, j, i) * y + prim(IVZ, k, j, i) * z);
}

/*
Enclosed mass profile of the three main components: SMBH, stars, dark matter
*/

static Real massFuncOfR(const Real distanceFromSMBH)
{
    Real haloMassFuncOfR = haloMass * (log(1 + distanceFromSMBH / haloRadius) - (distanceFromSMBH / (haloRadius + distanceFromSMBH)));
    Real stellarMassFuncOfR = stellarMass * (log(1. + distanceFromSMBH / stellarRadius) - (distanceFromSMBH / (stellarRadius + distanceFromSMBH)));
    return haloMassFuncOfR + SMBHMass + stellarMassFuncOfR;
}

/*
Gravitational Acceleration resulting from the enclosed mass
*/

static Real gravitationalAcceleration(const Real distanceFromSMBH)
{
    Real gravitationalAccelerationFuncOfR = -gravitationalConstAstronomical * massFuncOfR(distanceFromSMBH) / SQR(distanceFromSMBH);

    if (distanceFromSMBH <= innerRadius) // softening by a potential factor, Eq 13 of Guo+23
    {
        Real r_a = 0.5 * innerRadius;
        Real factorU = exp(-SQR(distanceFromSMBH / r_a));
        gravitationalAccelerationFuncOfR *= pow(distanceFromSMBH, 3) * (1 - factorU) / pow(SQR(distanceFromSMBH) + SQR(r_a) * factorU, 1.5);
    }
    return gravitationalAccelerationFuncOfR;
}

// END GRAVITY--------------------------------------------------------------------------------------------

// INITIAL CONDITIONS-------------------------------------------------------------------------------------

/*
Flat core entropy profile
*/

static Real coreEntropy(const Real dimensionlessX)
{
    return scaledEntropyAstronomical / 2 * (1 + pow(dimensionlessX, etaPower));
}

/*
Density solver, including a 4th-order Runge-Kutta scheme
*/

// dN/dX function, solved analytically
static Real dN_dX(const Real dimensionlessX, const Real numberDensity)
{
    return (2 * pow(numberDensity, 2 - gammaAdiabatic) * scaledRadius * gravitationalAcceleration(dimensionlessX * scaledRadius) * mu * hydrogenMassAstronomical / scaledEntropyAstronomical - numberDensity * pow(dimensionlessX, etaPower - 1) * etaPower) / ((1 + pow(dimensionlessX, etaPower)) * gammaAdiabatic);
}

// KR4 differentiator
static Real rungeKutta4(Real x0, Real y0, Real xn, int n, Real (*differentialFunction)(Real, Real))
{ // (x0,y0) are initial conditions, xn is the final value of x, n is the number of iterations, and differentialFunction is the dy/dx expression
    Real k1, k2, k3, k4;
    Real h = (xn - x0) / n;

    for (int i = 0; i < n; i++)
    {
        k1 = h * differentialFunction(x0, y0);
        k2 = h * differentialFunction(x0 + 0.5 * h, y0 + 0.5 * k1);
        k3 = h * differentialFunction(x0 + 0.5 * h, y0 + 0.5 * k2);
        k4 = h * differentialFunction(x0 + h, y0 + k3);

        y0 = y0 + (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
        x0 = x0 + h;
    }

    return y0;
}

// Find number density at given radial distance, using the KR4 approximation
static Real numberDensityFromEntropy(const Real dimensionlessX)
{
    int numOfIterations = 10000;

    return rungeKutta4(1.0, initialNumberDensityAstronomical, dimensionlessX, numOfIterations, &dN_dX);
}

/*
Gas pressure profile
*/

static Real pressureFromEntropyAndNumberDensity(const Real entropy, const Real numberDensity)
{
    return entropy * pow(numberDensity, gammaAdiabatic);
}

/*
Produce a table of normalized radial distances (dimensionlessX), densities, and pressure
*/

static void tabulatedNumberDensitiesAndPressures()
{
    Real a = log10(simulationBoxWidth * codeLength / scaledRadius);                                                           // set largest radius = box width to account for ghost zones
    Real b = log10(simulationBoxWidth / pow(2., numberOfRefinementLevels) / (numberOfZones * 2) * codeLength / scaledRadius); // set smallest radius = 1/2 width of the smallest cell
    Real spacing = (a - b) / (numberOfTabulatedEntries - 1);

    for (int i = 0; i < numberOfTabulatedEntries; ++i)
    {
        dimensionlessXArray[i] = pow(10., b + spacing * i); // basically the python np.logspace() function
        numberDensityFuncOfXArray[i] = numberDensityFromEntropy(dimensionlessXArray[i]);
        entropyFuncOfXArray[i] = coreEntropy(dimensionlessXArray[i]);
        pressureFuncOfXArray[i] = pressureFromEntropyAndNumberDensity(entropyFuncOfXArray[i], numberDensityFuncOfXArray[i]);
    }
}

// END INITIAL CONDITIONS-------------------------------------------------------------------------

// COOLING----------------------------------------------------------------------------------------

/*
Cooling function, part of the source function
*/

// Source function copied from Guo+23, added 03/06
void coolingSourceFunction(const Real dt, const AthenaArray<Real> &prim,
                           AthenaArray<Real> &cons, int k, int j, int i)
{
    Real primitiveNumberDensityCode = prim(IDN, k, j, i) / (mu * hydrogenMassAstronomical / codeMass);
    Real primitiveNumberDensityAstronomical = primitiveNumberDensityCode * codeNumberDensity;
    Real primitiveTemperatureCode = prim(IPR, k, j, i) / (primitiveNumberDensityCode * codeBoltzmannConst); // trying to calculate temperature via pressure instead of energy to bypass thermal energy requirement

    Real conservedDensityCode = fmax(cons(IDN, k, j, i), densityFloorAstronomical / codeDensity); // take the bigger between cell density and floor density
    Real conservedNumberDensityCode = conservedDensityCode / (mu * hydrogenMassAstronomical / codeMass);
    Real kineticEnergyDensityCode = 0.5 * (SQR(cons(IM1, k, j, i)) + SQR(cons(IM2, k, j, i)) + SQR(cons(IM3, k, j, i))) / conservedDensityCode; // unit: mass / length / time^2 (checked)
    Real thermalEnergyDensityCode = cons(IEN, k, j, i) - kineticEnergyDensityCode;
    Real conservedTemperatureCode = gammaMinus1 * thermalEnergyDensityCode / (conservedNumberDensityCode * codeBoltzmannConst);

    Real emissivityAstronomical = emissivityFromTemperature(primitiveTemperatureCode * codeTemperature);
    // softening function for emissivity, from Guo+23 [added 03/24/2023]
    emissivityAstronomical *= exp(-1.0e1 * pow(temperatureFloor / codeTemperature / primitiveTemperatureCode, 4.0));
    Real coolingRateCode = coolingFunction(primitiveNumberDensityAstronomical, emissivityAstronomical) / (codeEnergyDensity / codeTime);
    Real coolingLossCode = coolingRateCode * dt;

    // Make sure cooling is not too fast
    if (conservedTemperatureCode <= temperatureFloor / codeTemperature)
    {
        coolingLossCode = 0.; // no cooling if temperature at the start is lower than floor temperature
    }
    else if (gammaMinus1 * (thermalEnergyDensityCode - coolingLossCode) / (conservedNumberDensityCode * codeBoltzmannConst) < temperatureFloor / codeTemperature)
    {
        coolingLossCode = thermalEnergyDensityCode - (conservedNumberDensityCode * codeBoltzmannConst) * (temperatureFloor / codeTemperature) / gammaMinus1; // cooling just enough to reach floor temperature if naive cooling goes below floor
    }

    cons(IEN, k, j, i) -= coolingLossCode;
}

/*
Emissivity Lambda as a function of temperature, a piece-wise function
*/

static Real emissivityFromTemperature(Real temperature)
{
    Real emissivityCGS, emissivityAstronomical;
    Real logTemperature = log10(temperature);

    // Modified 08/02/2023: interpolate all datapoints instead of if-else clause to speed up computation
    emissivityCGS = pow(10., interpolate1D(logTemperature, logTemperatureArray, logEmissivityHydroArray));

    emissivityAstronomical = emissivityCGS / (solarMassCGS * pow(MpcCGS, 5) * pow(MyrCGS, -3));

    return emissivityAstronomical;
}

/*
Cooling function as a function of number density and emissivity
*/

static Real coolingFunction(Real numberDensity, Real emissivity)
{
    return SQR(numberDensity) * emissivity;
}

// END COOLING-------------------------------------------------------------------------------------

// JET FEEDBACK------------------------------------------------------------------------------------

/*
Rewrite jet feedback source function, 04/19/2023
*/

void jetFeedbackSourceFunction(MeshBlock *pmb, AthenaArray<Real> &cons, const AthenaArray<Real> &prim,
                               int k, int j, int i, const Real dt,
                               Real z, Real y, Real x, Real polarDistance)
{
    Real cellHeightCode = pmb->pcoord->GetEdge3Length(k, j, i);
    pmb->pmy_mesh->ruser_mesh_data[7](0) = cellHeightCode;
    Real &jetPrecessionAnglePhi = pmb->pmy_mesh->ruser_mesh_data[7](1); // current azimuthal angle of jet
    Real &jetMassLaunchRateCode = pmb->pmy_mesh->ruser_mesh_data[5](0); // mass launch equals mass accreted

    if ((z < jetLaunchingHeight / codeLength + 0.5 * cellHeightCode) && (z > jetLaunchingHeight / codeLength - 0.5 * cellHeightCode))
    {                                                                                           // top jet
        Real jetPlatformVolumeCode = PI * SQR(jetLaunchingWidth / codeLength) * cellHeightCode; // density = mass / volume, where volume = pi * R^2 * h, and mass = mass rate * dt

        // Added 06/14/2023: adding the Gaussian profile of the jet, following Li & Bryan (2014)
        Real profileWeight = 1.;
        Real jetSmoothingRadiusCode = 0.35 * jetLaunchingWidth / codeLength;             // this gives a 99.5% CDF profile
        profileWeight = 4.152 * exp(-0.5 * SQR(polarDistance / jetSmoothingRadiusCode)); // the factor 4.152 is numerically calculated to ensure total mass of the jet base equal to accreted mass, will differs for different smoothing factor

        // Added 06/17/2023: rewrote again to be consistent with Li & Bryan (2014)
        Real cellKineticEnergyDensityCode = 0.5 * (SQR(cons(IM1, k, j, i)) + SQR(cons(IM2, k, j, i)) + SQR(cons(IM3, k, j, i))) / cons(IDN, k, j, i); // KE = p^2 / 2 rho
        Real cellThermalEnergyDensityCode = cons(IEN, k, j, i) - cellKineticEnergyDensityCode;

        // Added 07/04/2023: jet precession
        Real jetVelocityXAstronomical, jetVelocityYAstronomical, jetVelocityZAstronomical;

        jetVelocityXAstronomical = jetVelocityAstronomical * cos(jetPrecessionAnglePhi) * sin(jetPrecessionAngleTheta);
        jetVelocityYAstronomical = jetVelocityAstronomical * sin(jetPrecessionAnglePhi) * sin(jetPrecessionAngleTheta);
        jetVelocityZAstronomical = jetVelocityAstronomical * cos(jetPrecessionAngleTheta);

        Real jetDensityCode = 0.5 * dt * jetMassLaunchRateCode * profileWeight / jetPlatformVolumeCode;
        cons(IDN, k, j, i) += jetDensityCode;                                           // inject density
        cons(IM1, k, j, i) += jetDensityCode * jetVelocityXAstronomical / codeVelocity; // inject momentum
        cons(IM2, k, j, i) += jetDensityCode * jetVelocityYAstronomical / codeVelocity;
        cons(IM3, k, j, i) += jetDensityCode * jetVelocityZAstronomical / codeVelocity;

        Real jetKineticEnergyDensityCode = 0.5 * jetDensityCode * SQR(jetVelocityAstronomical / codeVelocity); // KE = 1/2 rho v^2
        Real jetThermalEnergyDensityCode = (1. - jetKineticFraction) / jetKineticFraction * jetKineticEnergyDensityCode;

        Real updatedKineticEnergyDensityCode = 0.5 * (SQR(cons(IM1, k, j, i)) + SQR(cons(IM2, k, j, i)) + SQR(cons(IM3, k, j, i))) / cons(IDN, k, j, i); // KE = p^2 / 2 rho
        Real updatedThermalEnergyDensityCode = cellThermalEnergyDensityCode + jetThermalEnergyDensityCode;

        cons(IEN, k, j, i) = updatedKineticEnergyDensityCode + updatedThermalEnergyDensityCode; // update energy
    }
    else if ((z > -1. * jetLaunchingHeight / codeLength - 0.5 * cellHeightCode) && (z < -1. * jetLaunchingHeight / codeLength + 0.5 * cellHeightCode))
    {                                                                                           // bottom jet
        Real jetPlatformVolumeCode = PI * SQR(jetLaunchingWidth / codeLength) * cellHeightCode; // density = mass / volume, where volume = pi * R^2 * h, and mass = mass rate * dt

        // Added 06/14/2023: adding the Gaussian profile of the jet, following Li & Bryan (2014)
        Real profileWeight = 1.;
        Real jetSmoothingRadiusCode = 0.35 * jetLaunchingWidth / codeLength;             // this gives a 99.5% CDF profile
        profileWeight = 4.152 * exp(-0.5 * SQR(polarDistance / jetSmoothingRadiusCode)); // the factor 4.152 is numerically calculated to ensure total mass of the jet base equal to accreted mass, will differs for different smoothing factor

        // Added 06/17/2023: rewrote again to be consistent with Li & Bryan (2014)
        Real cellKineticEnergyDensityCode = 0.5 * (SQR(cons(IM1, k, j, i)) + SQR(cons(IM2, k, j, i)) + SQR(cons(IM3, k, j, i))) / cons(IDN, k, j, i); // KE = p^2 / 2 rho
        Real cellThermalEnergyDensityCode = cons(IEN, k, j, i) - cellKineticEnergyDensityCode;

        // Added 07/04/2023: jet precession
        Real jetVelocityXAstronomical, jetVelocityYAstronomical, jetVelocityZAstronomical;

        jetVelocityXAstronomical = jetVelocityAstronomical * cos(jetPrecessionAnglePhi) * sin(jetPrecessionAngleTheta);
        jetVelocityYAstronomical = jetVelocityAstronomical * sin(jetPrecessionAnglePhi) * sin(jetPrecessionAngleTheta);
        jetVelocityZAstronomical = jetVelocityAstronomical * cos(jetPrecessionAngleTheta);

        Real jetDensityCode = 0.5 * dt * jetMassLaunchRateCode * profileWeight / jetPlatformVolumeCode;
        cons(IDN, k, j, i) += jetDensityCode;                                            // inject density
        cons(IM1, k, j, i) += -jetDensityCode * jetVelocityXAstronomical / codeVelocity; // inject momentum
        cons(IM2, k, j, i) += -jetDensityCode * jetVelocityYAstronomical / codeVelocity;
        cons(IM3, k, j, i) += -jetDensityCode * jetVelocityZAstronomical / codeVelocity;

        Real jetKineticEnergyDensityCode = 0.5 * jetDensityCode * SQR(jetVelocityAstronomical / codeVelocity); // KE = 1/2 rho v^2
        Real jetThermalEnergyDensityCode = (1. - jetKineticFraction) / jetKineticFraction * jetKineticEnergyDensityCode;

        Real updatedKineticEnergyDensityCode = 0.5 * (SQR(cons(IM1, k, j, i)) + SQR(cons(IM2, k, j, i)) + SQR(cons(IM3, k, j, i))) / cons(IDN, k, j, i); // KE = p^2 / 2 rho
        Real updatedThermalEnergyDensityCode = cellThermalEnergyDensityCode + jetThermalEnergyDensityCode;

        cons(IEN, k, j, i) = updatedKineticEnergyDensityCode + updatedThermalEnergyDensityCode; // update energy
    }
}

// END JET FEEDBACK--------------------------------------------------------------------------------

// MISC. FUNCTIONS---------------------------------------------------------------------------------

/*
CROSS PRODUCT
*/

static void crossProduct(std::vector<Real> r, std::vector<Real> p, Real *x, Real *y, Real *z)
{
    *x = r[1] * p[2] - r[2] * p[1];
    *y = r[2] * p[0] - r[0] * p[2];
    *z = r[0] * p[1] - r[1] * p[0];

    return;
}

/*
1D LINEAR INTERPOLATOR
*/

static Real interpolate1D(Real inputX, std::vector<Real> tabulatedX, std::vector<Real> tabulatedY)
{
    std::vector<Real>::iterator upperIndexPointer;

    upperIndexPointer = std::lower_bound(tabulatedX.begin(), tabulatedX.end(), inputX);
    int upperIndexX = (upperIndexPointer - tabulatedX.begin());
    int lowerIndexX = upperIndexX - 1;

    Real lowerBoundX = tabulatedX[lowerIndexX];
    Real upperBoundX = tabulatedX[upperIndexX];

    Real lowerBoundY = tabulatedY[lowerIndexX];
    Real upperBoundY = tabulatedY[upperIndexX];

    Real interpolatedY = (lowerBoundY * (upperBoundX - inputX) + upperBoundY * (inputX - lowerBoundX)) / (upperBoundX - lowerBoundX);

    return interpolatedY;
}

/*
Enforce a temperature and density floor if a cell has lower temperature than floor temperature
*/

void enforceFloors(MeshBlock *pmb, AthenaArray<Real> &cons, int k, int j, int i)
{
    Real numberDensityCode = cons(IDN, k, j, i) / (mu * hydrogenMassAstronomical / codeMass);
    Real kineticEnergyDensityCode = 0.5 * (SQR(cons(IM1, k, j, i)) + SQR(cons(IM2, k, j, i)) + SQR(cons(IM3, k, j, i))) / cons(IDN, k, j, i);
    Real temperatureCode = gammaMinus1 * (cons(IEN, k, j, i) - kineticEnergyDensityCode) / (numberDensityCode * codeBoltzmannConst);

    if (temperatureCode < temperatureFloor / codeTemperature)
    {
        Real floorThermalEnergyDensityCode = fmax(numberDensityCode, numberDensityFloorAstronomical / codeNumberDensity) * codeBoltzmannConst * (temperatureFloor / codeTemperature) / (gammaMinus1);
        cons(IEN, k, j, i) = floorThermalEnergyDensityCode + kineticEnergyDensityCode; // hard reset of the internal energy
    }

    if (numberDensityCode < numberDensityFloorAstronomical / codeNumberDensity)
    {
        cons(IDN, k, j, i) = numberDensityFloorAstronomical * (mu * hydrogenMassAstronomical) / codeDensity; // hard reset of the density -> might change temperature subsequently, to be observed
    }
}

/*
Cooling time, for refinement condition
*/

static Real coolingTime(AthenaArray<Real> &w, int k, int j, int i)
{ // function is 5/2 nkT / n^2 Lambda, return code units
    Real primitiveNumberDensityCode = w(IDN, k, j, i) / (mu * hydrogenMassAstronomical / codeMass);
    Real primitiveTemperatureAstronomical = w(IPR, k, j, i) / (primitiveNumberDensityCode * codeBoltzmannConst) * codeTemperature;

    return (5. / 2.) * boltzmannConstAstronomical * primitiveTemperatureAstronomical / (primitiveNumberDensityCode * codeNumberDensity * emissivityFromTemperature(primitiveTemperatureAstronomical)) / codeTime;
}

/*
Sound crossing time, for refinement condition
*/

static Real soundCrossingTime(AthenaArray<Real> &w, int k, int j, int i, Real cellWidth)
{ // function is delta_x / c_s, where c_s = sqrt(gamma * P / rho), return code units
    Real soundSpeedCode = sqrt(gammaAdiabatic * w(IPR, k, j, i) / w(IDN, k, j, i));
    return cellWidth / soundSpeedCode;
}

// END MISC. FUNCTIONS-----------------------------------------------------------------------------

// USER HISTORY OUTPUT-----------------------------------------------------------------------------

/*
Output accretion rate and jet power data
*/

Real coldAccretionRate(MeshBlock *pmb, int iout)
{ // unit is Msun yr^-1
    return (pmb->pmy_mesh->ruser_mesh_data[5](0)) * codeMass / codeTime / 1.e6;
}

Real hotAccretionRate(MeshBlock *pmb, int iout)
{ // unit is Msun yr^-1
    return (pmb->pmy_mesh->ruser_mesh_data[5](1)) * codeMass / codeTime / 1.e6;
}

Real jetPower(MeshBlock *pmb, int iout)
{ // unit is erg s^-1
    return pmb->pmy_mesh->ruser_mesh_data[6](0) * codePower * (solarMassCGS * pow(MpcCGS, 2) * pow(MyrCGS, -3));
}

/*
CALCULATE TOTAL ANGULAR MOMENTUM OF THE SYSTEM IN EACH DIRECTION, NOT WRITTEN YET
*/

Real angularMomentumX(MeshBlock *pmb, int iout)
{ // total angular momentum of the system inside the jet region (exclude jet region)

    return 0.;
}

Real angularMomentumY(MeshBlock *pmb, int iout)
{
    return 0.;
}

Real angularMomentumZ(MeshBlock *pmb, int iout)
{
    return 0.;
}
// END USER HISTORY OUTPUT-------------------------------------------------------------------------