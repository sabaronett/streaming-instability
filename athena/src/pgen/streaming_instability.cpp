//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file streaming_instability.cpp
//  \brief sets up a linear mode for the streaming instability.

// NOTE: In this setup, Y <-> Z.

// C++ standard libraries
#include <cmath>   // sqrt(), pow(), round(), sin(), cos()
#include <limits>  // numeric_limits
#include <random>
#include <vector>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../hydro/hydro.hpp"
#include "../parameter_input.hpp"
#include "../particles/particles.hpp"
#include "../particles/dust_particles.hpp"

// Input parameters
static Real omega(1.0);
static Real epsilon(0.0);

// Global variables
static Real vol(0.0);  // total volume of the domain
static Real cs0(0.0);  // speed of sound
static Real ux0(0.0), uy0(0.0);  // equilibrium gas velocity
static Real vx0(0.0), vy0(0.0);  // equilibrium particle velocity
static Real two_omega(2.0 * omega);
static Real omega_half(0.5 * omega);
static Real gas_accel_x(0.0);

//======================================================================================
//! \fn Real Drhog(MeshBlock *pmb, int iout)
//  \brief Finds the first moment of the gas density deviation.
//======================================================================================
Real Drhog(MeshBlock *pmb, int iout)
{
  // Get constants.
  const int is(pmb->is), ie(pmb->ie);
  const int js(pmb->js), je(pmb->je);
  const int ks(pmb->ks), ke(pmb->ke);

  // Get the gas density field.
  AthenaArray<Real> rho;
  rho.InitWithShallowSlice(pmb->phydro->u, 4, IDN, 1);

  // Integrate the deviation.
  Coordinates *pcoord(pmb->pcoord);
  Real integral(0.0);
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i)
        integral += (rho(k,j,i) - 1.0) * pcoord->GetCellVolume(k,j,i);

  // Divide the total volume and return the moment.
  return integral / vol;
}

//======================================================================================
//! \fn Real Drhog2(MeshBlock *pmb, int iout)
//  \brief Finds the second moment of the gas density deviation.
//======================================================================================
Real Drhog2(MeshBlock *pmb, int iout)
{
  // Get constants.
  const int is(pmb->is), ie(pmb->ie);
  const int js(pmb->js), je(pmb->je);
  const int ks(pmb->ks), ke(pmb->ke);

  // Get the gas density field.
  AthenaArray<Real> rho;
  rho.InitWithShallowSlice(pmb->phydro->u, 4, IDN, 1);

  // Integrate the deviation.
  Coordinates *pcoord(pmb->pcoord);
  Real integral(0.0);
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real drho(rho(k,j,i) - 1.0);
        integral += drho * drho * pcoord->GetCellVolume(k,j,i);
      }

  // Divide the total volume and return the moment.
  return integral / vol;
}

//======================================================================================
//! \fn Real Dux(MeshBlock *pmb, int iout)
//  \brief Finds the first moment of the radial velocity deviation of the gas.
//======================================================================================
Real Dux(MeshBlock *pmb, int iout)
{
  // Get constants.
  const int is(pmb->is), ie(pmb->ie);
  const int js(pmb->js), je(pmb->je);
  const int ks(pmb->ks), ke(pmb->ke);

  // Get the gas density and radial velocity.
  AthenaArray<Real> rho, ux;
  rho.InitWithShallowSlice(pmb->phydro->u, 4, IDN, 1);
  ux.InitWithShallowSlice(pmb->phydro->w, 4, IVX, 1);

  // Integrate the deviation.
  Coordinates *pcoord(pmb->pcoord);
  Real integral(0.0);
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i)
        integral += rho(k,j,i) * (ux(k,j,i) - ux0) * pcoord->GetCellVolume(k,j,i);

  // Return the moment.
  return integral;
}

//======================================================================================
//! \fn Real Dux2(MeshBlock *pmb, int iout)
//  \brief Finds the second moment of the radial velocity deviation of the gas.
//======================================================================================
Real Dux2(MeshBlock *pmb, int iout)
{
  // Get constants.
  const int is(pmb->is), ie(pmb->ie);
  const int js(pmb->js), je(pmb->je);
  const int ks(pmb->ks), ke(pmb->ke);

  // Get the gas density and radial velocity.
  AthenaArray<Real> rho, ux;
  rho.InitWithShallowSlice(pmb->phydro->u, 4, IDN, 1);
  ux.InitWithShallowSlice(pmb->phydro->w, 4, IVX, 1);

  // Integrate the deviation.
  Coordinates *pcoord(pmb->pcoord);
  Real integral(0.0);
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real dux(ux(k,j,i) - ux0);
        integral += rho(k,j,i) * (dux * dux) * pcoord->GetCellVolume(k,j,i);
      }

  // Return the moment.
  return integral;
}

//======================================================================================
//! \fn Real Duy(MeshBlock *pmb, int iout)
//  \brief Finds the first moment of the azimuthal velocity deviation of the gas.
//======================================================================================
Real Duy(MeshBlock *pmb, int iout)
{
  // Get constants.
  const int is(pmb->is), ie(pmb->ie);
  const int js(pmb->js), je(pmb->je);
  const int ks(pmb->ks), ke(pmb->ke);

  // Get the gas density and azimuthal velocity.
  AthenaArray<Real> rho, uy;
  rho.InitWithShallowSlice(pmb->phydro->u, 4, IDN, 1);
  uy.InitWithShallowSlice(pmb->phydro->w, 4, IVZ, 1);

  // Integrate the deviation.
  Coordinates *pcoord(pmb->pcoord);
  Real integral(0.0);
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i)
        integral += rho(k,j,i) * (uy(k,j,i) - uy0) * pcoord->GetCellVolume(k,j,i);

  // Return the moment.
  return integral;
}

//======================================================================================
//! \fn Real Duy2(MeshBlock *pmb, int iout)
//  \brief Finds the second moment of the azimuthal velocity deviation of the gas.
//======================================================================================
Real Duy2(MeshBlock *pmb, int iout)
{
  // Get constants.
  const int is(pmb->is), ie(pmb->ie);
  const int js(pmb->js), je(pmb->je);
  const int ks(pmb->ks), ke(pmb->ke);

  // Get the gas density and azimuthal velocity.
  AthenaArray<Real> rho, uy;
  rho.InitWithShallowSlice(pmb->phydro->u, 4, IDN, 1);
  uy.InitWithShallowSlice(pmb->phydro->w, 4, IVZ, 1);

  // Integrate the deviation.
  Coordinates *pcoord(pmb->pcoord);
  Real integral(0.0);
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real duy(uy(k,j,i) - uy0);
        integral += rho(k,j,i) * (duy * duy) * pcoord->GetCellVolume(k,j,i);
      }

  // Return the moment.
  return integral;
}

//======================================================================================
//! \fn Real DuxDuy(MeshBlock *pmb, int iout)
//  \brief Finds the product of the first moments of the radial and azimuthal velocity
//         deviations of the gas.
//======================================================================================
Real DuxDuy(MeshBlock *pmb, int iout)
{
  // Get constants.
  const int is(pmb->is), ie(pmb->ie);
  const int js(pmb->js), je(pmb->je);
  const int ks(pmb->ks), ke(pmb->ke);

  // Get the gas density and azimuthal velocity.
  AthenaArray<Real> rho, ux, uy;
  rho.InitWithShallowSlice(pmb->phydro->u, 4, IDN, 1);
  ux.InitWithShallowSlice(pmb->phydro->w, 4, IVX, 1);
  uy.InitWithShallowSlice(pmb->phydro->w, 4, IVZ, 1);

  // Integrate the deviation.
  Coordinates *pcoord(pmb->pcoord);
  Real integral(0.0);
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real dux(ux(k,j,i) - ux0);
        Real duy(uy(k,j,i) - uy0);
        integral += rho(k,j,i) * (dux * duy) * pcoord->GetCellVolume(k,j,i);
      }

  // Return the product of the moments.
  return integral;
}

//======================================================================================
//! \fn Real DuxDuz(MeshBlock *pmb, int iout)
//  \brief Finds the product of the first moments of the radial and vertical velocity
//         deviations of the gas.
//======================================================================================
Real DuxDuz(MeshBlock *pmb, int iout)
{
  // Get constants.
  const int is(pmb->is), ie(pmb->ie);
  const int js(pmb->js), je(pmb->je);
  const int ks(pmb->ks), ke(pmb->ke);

  // Get the gas density and azimuthal velocity.
  AthenaArray<Real> rho, ux, uz;
  rho.InitWithShallowSlice(pmb->phydro->u, 4, IDN, 1);
  ux.InitWithShallowSlice(pmb->phydro->w, 4, IVX, 1);
  uz.InitWithShallowSlice(pmb->phydro->w, 4, IVY, 1);

  // Integrate the deviation.
  Coordinates *pcoord(pmb->pcoord);
  Real integral(0.0);
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real dux(ux(k,j,i) - ux0);
        integral += rho(k,j,i) * (dux * uz(k,j,i)) * pcoord->GetCellVolume(k,j,i);
      }

  // Return the product of the moments.
  return integral;
}

//======================================================================================
//! \fn Real DuxDuy(MeshBlock *pmb, int iout)
//  \brief Finds the product of the first moments of the azimuthal and vertical velocity
//         deviations of the gas.
//======================================================================================
Real DuyDuz(MeshBlock *pmb, int iout)
{
  // Get constants.
  const int is(pmb->is), ie(pmb->ie);
  const int js(pmb->js), je(pmb->je);
  const int ks(pmb->ks), ke(pmb->ke);

  // Get the gas density and azimuthal velocity.
  AthenaArray<Real> rho, uy, uz;
  rho.InitWithShallowSlice(pmb->phydro->u, 4, IDN, 1);
  uy.InitWithShallowSlice(pmb->phydro->w, 4, IVZ, 1);
  uz.InitWithShallowSlice(pmb->phydro->w, 4, IVY, 1);

  // Integrate the deviation.
  Coordinates *pcoord(pmb->pcoord);
  Real integral(0.0);
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real duy(uy(k,j,i) - uy0);
        integral += rho(k,j,i) * (duy * uz(k,j,i)) * pcoord->GetCellVolume(k,j,i);
      }

  // Return the product of the moments.
  return integral;
}

//======================================================================================
//! \fn Real Drhop(MeshBlock *pmb, int iout)
//  \brief Finds the first moment of the particle density deviation.
//======================================================================================
Real Drhop(MeshBlock *pmb, int iout)
{
  // Get constants.
  const int is(pmb->is), ie(pmb->ie);
  const int js(pmb->js), je(pmb->je);
  const int ks(pmb->ks), ke(pmb->ke);

  // Get the gas density field.
  AthenaArray<Real> rhop;
  rhop = pmb->ppar->GetMassDensity();

  // Integrate the deviation.
  Coordinates *pcoord(pmb->pcoord);
  Real integral(0.0);
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i)
        integral += (rhop(k,j,i) - epsilon) * pcoord->GetCellVolume(k,j,i);

  // Divide the total volume and return the moment.
  return integral / vol;
}

//======================================================================================
//! \fn Real Drhop2(MeshBlock *pmb, int iout)
//  \brief Finds the second moment of the particle density deviation.
//======================================================================================
Real Drhop2(MeshBlock *pmb, int iout)
{
  // Get constants.
  const int is(pmb->is), ie(pmb->ie);
  const int js(pmb->js), je(pmb->je);
  const int ks(pmb->ks), ke(pmb->ke);

  // Get the particle density field.
  AthenaArray<Real> rhop;
  rhop = pmb->ppar->GetMassDensity();

  // Integrate the deviation.
  Coordinates *pcoord(pmb->pcoord);
  Real integral(0.0);
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real drhop(rhop(k,j,i) - epsilon);
        integral += drhop * drhop * pcoord->GetCellVolume(k,j,i);
      }

  // Divide the total volume and return the moment.
  return integral / vol;
}

//======================================================================================
//! \fn Real Dvpx(MeshBlock *pmb, int iout)
//  \brief Finds the first moment of the radial velocity deviation of the particles.
//======================================================================================
Real Dvpx(MeshBlock *pmb, int iout)
{
  // Get constants.
  const int is(pmb->is), ie(pmb->ie);
  const int js(pmb->js), je(pmb->je);
  const int ks(pmb->ks), ke(pmb->ke);

  // Get the particle density and radial velocity.
  AthenaArray<Real> rhop, vp, vx;
  rhop = pmb->ppar->GetMassDensity();
  vp = pmb->ppar->GetVelocityField();
  vx.InitWithShallowSlice(vp, 4, 0, 1);

  // Integrate the deviation.
  Coordinates *pcoord(pmb->pcoord);
  Real integral(0.0);
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i)
        integral += rhop(k,j,i) * (vx(k,j,i) - vx0) * pcoord->GetCellVolume(k,j,i);

  // Return the moment.
  return integral;
}

//======================================================================================
//! \fn Real Dvpx2(MeshBlock *pmb, int iout)
//  \brief Finds the second moment of the radial velocity deviation of the particles.
//======================================================================================
Real Dvpx2(MeshBlock *pmb, int iout)
{
  // Get constants.
  const int is(pmb->is), ie(pmb->ie);
  const int js(pmb->js), je(pmb->je);
  const int ks(pmb->ks), ke(pmb->ke);

  // Get the particle density and radial velocity.
  AthenaArray<Real> rhop, vp, vx;
  rhop = pmb->ppar->GetMassDensity();
  vp = pmb->ppar->GetVelocityField();
  vx.InitWithShallowSlice(vp, 4, 0, 1);

  // Integrate the deviation.
  Coordinates *pcoord(pmb->pcoord);
  Real integral(0.0);
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real dvx(vx(k,j,i) - vx0);
        integral += rhop(k,j,i) * (dvx * dvx) * pcoord->GetCellVolume(k,j,i);
      }

  // Return the moment.
  return integral;
}

//======================================================================================
//! \fn Real Dvpy(MeshBlock *pmb, int iout)
//  \brief Finds the first moment of the azimuthal velocity deviation of the particles.
//======================================================================================
Real Dvpy(MeshBlock *pmb, int iout)
{
  // Get constants.
  const int is(pmb->is), ie(pmb->ie);
  const int js(pmb->js), je(pmb->je);
  const int ks(pmb->ks), ke(pmb->ke);

  // Get the particle density and azimuthal velocity.
  AthenaArray<Real> rhop, vp, vy;
  rhop = pmb->ppar->GetMassDensity();
  vp = pmb->ppar->GetVelocityField();
  vy.InitWithShallowSlice(vp, 4, 2, 1);

  // Integrate the deviation.
  Coordinates *pcoord(pmb->pcoord);
  Real integral(0.0);
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i)
        integral += rhop(k,j,i) * (vy(k,j,i) - vy0) * pcoord->GetCellVolume(k,j,i);

  // Return the moment.
  return integral;
}

//======================================================================================
//! \fn Real Dvpy2(MeshBlock *pmb, int iout)
//  \brief Finds the second moment of the azimuthal velocity deviation of the particles.
//======================================================================================
Real Dvpy2(MeshBlock *pmb, int iout)
{
  // Get constants.
  const int is(pmb->is), ie(pmb->ie);
  const int js(pmb->js), je(pmb->je);
  const int ks(pmb->ks), ke(pmb->ke);

  // Get the particle density and azimuthal velocity.
  AthenaArray<Real> rhop, vp, vy;
  rhop = pmb->ppar->GetMassDensity();
  vp = pmb->ppar->GetVelocityField();
  vy.InitWithShallowSlice(vp, 4, 2, 1);

  // Integrate the deviation.
  Coordinates *pcoord(pmb->pcoord);
  Real integral(0.0);
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real dvy(vy(k,j,i) - vy0);
        integral += rhop(k,j,i) * (dvy * dvy) * pcoord->GetCellVolume(k,j,i);
      }

  // Return the moment.
  return integral;
}

//======================================================================================
//! \fn Real Dvpz(MeshBlock *pmb, int iout)
//  \brief Finds the first moment of the vertical velocity deviation of the particles.
//======================================================================================
Real Dvpz(MeshBlock *pmb, int iout)
{
  // Get constants.
  const int is(pmb->is), ie(pmb->ie);
  const int js(pmb->js), je(pmb->je);
  const int ks(pmb->ks), ke(pmb->ke);

  // Get the particle density and vertical velocity.
  AthenaArray<Real> rhop, vp, vz;
  rhop = pmb->ppar->GetMassDensity();
  vp = pmb->ppar->GetVelocityField();
  vz.InitWithShallowSlice(vp, 4, 1, 1);

  // Integrate the deviation.
  Coordinates *pcoord(pmb->pcoord);
  Real integral(0.0);
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i)
        integral += rhop(k,j,i) * vz(k,j,i) * pcoord->GetCellVolume(k,j,i);

  // Return the moment.
  return integral;
}

//======================================================================================
//! \fn Real Dvpz2(MeshBlock *pmb, int iout)
//  \brief Finds the second moment of the vertical velocity deviation of the particles.
//======================================================================================
Real Dvpz2(MeshBlock *pmb, int iout)
{
  // Get constants.
  const int is(pmb->is), ie(pmb->ie);
  const int js(pmb->js), je(pmb->je);
  const int ks(pmb->ks), ke(pmb->ke);

  // Get the particle density and vertical velocity.
  AthenaArray<Real> rhop, vp, vz;
  rhop = pmb->ppar->GetMassDensity();
  vp = pmb->ppar->GetVelocityField();
  vz.InitWithShallowSlice(vp, 4, 1, 1);

  // Integrate the deviation.
  Coordinates *pcoord(pmb->pcoord);
  Real integral(0.0);
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real dvz(vz(k,j,i));
        integral += rhop(k,j,i) * (dvz * dvz) * pcoord->GetCellVolume(k,j,i);
      }

  // Return the moment.
  return integral;
}

//======================================================================================
//! \fn void SourceTermsForGas(MeshBlock *pmb, const Real time, const Real dt,
//       const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
//       const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
//       AthenaArray<Real> &cons_scalar)
//  \brief Adds source terms to the gas.
//======================================================================================
void SourceTermsForGas(MeshBlock *pmb, const Real time, const Real dt,
         const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
         const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
         AthenaArray<Real> &cons_scalar) {
  const int is(pmb->is), js(pmb->js), ks(pmb->ks);
  const int ie(pmb->ie), je(pmb->je), ke(pmb->ke);

  // Apply the Coriolis and centrifugal forces, and linear gravity from the star, and
  // the background radial pressure gradient.
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real rho_dt = prim(IDN,k,j,i) * dt;
        cons(IM1,k,j,i) += rho_dt * (two_omega * prim(IVZ,k,j,i) + gas_accel_x);
        cons(IM3,k,j,i) -= rho_dt * omega_half * prim(IVX,k,j,i);
      }
    }
  }
}

//======================================================================================
//! \fn void DustParticles::UserSourceTerms(Real t, Real dt,
//                              const AthenaArray<Real>& meshsrc)
//  \brief Adds source terms to the particles.
//======================================================================================
void DustParticles::UserSourceTerms(Real t, Real dt, const AthenaArray<Real>& meshsrc) {
  // Apply the Coriolis and centrifugal forces, and linear gravity from the star.
  for (int k = 0; k < npar; ++k) {
    dvpx[k] += two_omega * vpz[k];
    dvpz[k] -= omega_half * vpx[k];
  }
}

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Initializes problem-specific data in Mesh class.
//========================================================================================
void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Preprocess constants.
  omega = pin->GetOrAddReal("problem", "omega", omega);
  epsilon = pin->GetReal("problem", "epsilon");
  cs0 = pin->GetReal("hydro", "iso_sound_speed");
  Real duy0(cs0 * pin->GetReal("problem", "duy0"));
  two_omega = 2.0 * omega;
  omega_half = 0.5 * omega;
  gas_accel_x = two_omega * duy0;

  // Find the total volume of the domain.
  vol = mesh_size.x1len * mesh_size.x2len * mesh_size.x3len;

  // Find the Nakagawa-Sekiya-Hayashi (1986) equilibrium solution.
  Real taus(omega * DustParticles::GetStoppingTime());
  Real v(duy0 / (std::pow(1.0 + epsilon, 2) + std::pow(taus, 2)));
  ux0 = 2.0 * epsilon * taus * v;
  uy0 = -((1.0 + epsilon) + std::pow(taus, 2)) * v;
  vx0 = -2.0 * taus * v;
  vy0 = -(1.0 + epsilon) * v;

  // Enroll source terms and time step function.
  EnrollUserExplicitSourceFunction(SourceTermsForGas);

  // Enroll additional history outputs.
  AllocateUserHistoryOutput(17);
  EnrollUserHistoryOutput(0, Drhog, "drhog");
  EnrollUserHistoryOutput(1, Drhog2, "drhog2");
  EnrollUserHistoryOutput(2, Dux, "dux");
  EnrollUserHistoryOutput(3, Dux2, "dux2");
  EnrollUserHistoryOutput(4, Duy, "duy");
  EnrollUserHistoryOutput(5, Duy2, "duy2");
  EnrollUserHistoryOutput(6, DuxDuy, "duxduy");
  EnrollUserHistoryOutput(7, DuxDuz, "duxduz");
  EnrollUserHistoryOutput(8, DuyDuz, "duyduz");
  EnrollUserHistoryOutput(9, Drhop, "drhop");
  EnrollUserHistoryOutput(10, Drhop2, "drhop2");
  EnrollUserHistoryOutput(11, Dvpx, "dvpx");
  EnrollUserHistoryOutput(12, Dvpx2, "dvpx2");
  EnrollUserHistoryOutput(13, Dvpy, "dvpy");
  EnrollUserHistoryOutput(14, Dvpy2, "dvpy2");
  EnrollUserHistoryOutput(15, Dvpz, "dvpz");
  EnrollUserHistoryOutput(16, Dvpz2, "dvpz2");
}

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Sets the initial conditions.
//======================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // Find the total number of particles in each direction.
  RegionSize& mesh_size = pmy_mesh->mesh_size;
  int npx1(block_size.nx1 > 1 ?
           pin->GetOrAddInteger("problem", "npx1", mesh_size.nx1) : 1),
      npx2(block_size.nx2 > 1 ?
           pin->GetOrAddInteger("problem", "npx2", mesh_size.nx2) : 1),
      npx3(block_size.nx3 > 1 ?
           pin->GetOrAddInteger("problem", "npx3", mesh_size.nx3) : 1);

  // Find the mass of each particle and the distance between adjacent particles.
  Real dx1(mesh_size.x1len / npx1),
       dx2(mesh_size.x2len / npx2),
       dx3(mesh_size.x3len / npx3);
  DustParticles::SetOneParticleMass(epsilon * vol / (npx1 * npx2 * npx3));

  // Determine number of particles in the block.
  int npx1_loc(static_cast<int>(std::round(block_size.x1len / dx1))),
      npx2_loc(static_cast<int>(std::round(block_size.x2len / dx2))),
      npx3_loc(static_cast<int>(std::round(block_size.x3len / dx3)));
  int npar(npx1_loc * npx2_loc * npx3_loc);
  ppar->Resize(npar);

  // Get the dimensionless stopping time.
  Real taus(DustParticles::GetStoppingTime());
  DustParticles *pdp(static_cast<DustParticles*>(ppar));
  if (DustParticles::GetVariableTaus()) {
    for (int k = 0; k < npar; ++k)
      pdp->taus[k] = taus;
  }
  taus *= omega;

  if (pin->GetOrAddBoolean("problem", "randparpos", false)) { // Nonlinear model
    // Make the gas in equilibrium.
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is; i <= ie; ++i) {
          phydro->u(IDN,k,j,i) = 1.0;
          phydro->u(IM1,k,j,i) = ux0;
          phydro->u(IM3,k,j,i) = uy0;
          phydro->u(IM2,k,j,i) = 0.0;
        }
      }
    }

    // Randomly distribute particles.
    Real zp1 = 0.5 * (block_size.x3max + block_size.x3min);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<Real> x1rand(block_size.x1min, block_size.x1max);
    std::uniform_real_distribution<Real> x2rand(block_size.x2min, block_size.x2max);
    std::uniform_real_distribution<Real> x3rand(block_size.x3min, block_size.x3max);
    for (int k = 0; k < npar; ++k) {
      ppar->xp[k] = x1rand(gen);
      ppar->yp[k] = x2rand(gen);
      ppar->zp[k] = block_size.nx3 > 1 ? x3rand(gen) : zp1;
      ppar->vpx[k] = vx0;
      ppar->vpy[k] = 0.0;
      ppar->vpz[k] = vy0;
    }

  } else { // Linear model
    // Get the wavenumber.
    Real duy0(cs0 * pin->GetReal("problem", "duy0"));
    Real length(duy0 / omega);
    Real kx(pin->GetReal("problem", "kx") / length),
         kz(pin->GetReal("problem", "kz") / length);

    // Find the scale factor for the eigenvector.
    Real amp(pin->GetOrAddReal("problem", "amplitude", 1e-6));
    Real drhop_re(pin->GetOrAddReal("problem", "drhop_re", 1.0)),
         drhop_im(pin->GetOrAddReal("problem", "drhop_im", 0.0));
    amp /= std::sqrt(drhop_re * drhop_re + drhop_im * drhop_im);

    // Perturb the gas.
    Real drhog_re(epsilon * amp * pin->GetReal("problem", "drhog_re")),
         drhog_im(epsilon * amp * pin->GetReal("problem", "drhog_im"));
    Real dv(epsilon * amp * duy0);
    Real dux_re(dv * pin->GetReal("problem", "dux_re")),
         dux_im(dv * pin->GetReal("problem", "dux_im")),
         duy_re(dv * pin->GetReal("problem", "duy_re")),
         duy_im(dv * pin->GetReal("problem", "duy_im")),
         duz_re(dv * pin->GetReal("problem", "duz_re")),
         duz_im(dv * pin->GetReal("problem", "duz_im"));

    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        Real coskz(std::cos(kz * pcoord->x2v(j))), sinkz(std::sin(kz * pcoord->x2v(j)));
        for (int i = is; i <= ie; ++i) {
          Real coskx(std::cos(kx * pcoord->x1v(i))), sinkx(std::sin(kx * pcoord->x1v(i)));
          Real rhog(1.0 + (drhog_re * coskx - drhog_im * sinkx) * coskz);
          phydro->u(IDN,k,j,i) = rhog;
          phydro->u(IM1,k,j,i) = rhog * (ux0 + (dux_re * coskx - dux_im * sinkx) * coskz);
          phydro->u(IM3,k,j,i) = rhog * (uy0 + (duy_re * coskx - duy_im * sinkx) * coskz);
          phydro->u(IM2,k,j,i) = -rhog * (duz_im * coskx + duz_re * sinkx) * sinkz;
        }
      }
    }

    // Uniformly distribute particles.
    std::vector<Real> xp0(npx1_loc), zp0(npx2_loc), yp0(npx3_loc);
    std::vector<Real> argx(npx1_loc), argz(npx2_loc);
    std::vector<Real> cos2kx(npx1_loc), sin2kx(npx1_loc), sin2kz(npx2_loc);
    for (int i = 0; i < npx1_loc; ++i) {
      Real x(block_size.x1min + (i + 0.5) * dx1), arg(kx * x);
      xp0[i] = x;
      argx[i] = arg;
      arg *= 2.0;
      cos2kx[i] = std::cos(arg);
      sin2kx[i] = std::sin(arg);
    }
    for (int j = 0; j < npx2_loc; ++j) {
      Real z(block_size.x2min + (j + 0.5) * dx2), arg(kz * z);
      zp0[j] = z;
      argz[j] = arg;
      sin2kz[j] = std::sin(2.0 * arg);
    }
    for (int k = 0; k < npx3_loc; ++k)
      yp0[k] = block_size.x3min + (k + 0.5) * dx3;

    // Compute repeated constants.
    Real c1x(kx * kx + kz * kz), c2x(c1x * c1x);
    if (c1x > 0.0) {
      c1x = 0.5 / c1x;
      c2x = 1.0 / c2x;
    }
    Real c1z(c1x * kz), c2z(c2x * std::pow(kz, 3));
    c1x *= kx;
    c2x *= std::pow(kx, 3);

    drhop_re *= amp;
    drhop_im *= amp;

    Real a1(0.25 * (drhop_re * drhop_re - drhop_im * drhop_im)),
         a2(0.5 * drhop_re * drhop_im),
         a3(0.25 * (drhop_re * drhop_re + drhop_im * drhop_im));

    // Assign particle positions.
    int ipar(0);
    for (int k = 0; k < npx3_loc; ++k) {
      for (int j = 0; j < npx2_loc; ++j) {
        for (int i = 0; i < npx1_loc; ++i) {
          Real kp(argx[i] + argz[j]), km(argx[i] - argz[j]);
          Real sinp(std::sin(kp)), sinm(std::sin(km));
          Real cosp(std::cos(kp)), cosm(std::cos(km));
          Real sinp2(std::sin(2.0 * kp)), sinm2(std::sin(2.0 * km));
          Real cosp2(std::cos(2.0 * kp)), cosm2(std::cos(2.0 * km));

          Real dxp(-c1x * (drhop_re * (sinp + sinm) + drhop_im * (cosp + cosm) -
                           a1 * (sinp2 + sinm2) - a2 * (cosp2 + cosm2)) +
                    c2x * (a2 * cos2kx[i] + a1 * sin2kx[i]));
          Real dzp(-c1z * (drhop_re * (sinp - sinm) + drhop_im * (cosp - cosm) -
                           a1 * (sinp2 - sinm2) - a2 * (cosp2 - cosm2)) +
                    c2z * a3 * sin2kz[j]);

          ppar->xp[ipar] = xp0[i] + dxp;
          ppar->zp[ipar] = yp0[k];
          ppar->yp[ipar] = zp0[j] + dzp;
          ++ipar;
        }
      }
    }

    // Assign the particle velocities.
    Real dvpx_re(dv * pin->GetReal("problem", "dvpx_re")),
         dvpx_im(dv * pin->GetReal("problem", "dvpx_im")),
         dvpy_re(dv * pin->GetReal("problem", "dvpy_re")),
         dvpy_im(dv * pin->GetReal("problem", "dvpy_im")),
         dvpz_re(dv * pin->GetReal("problem", "dvpz_re")),
         dvpz_im(dv * pin->GetReal("problem", "dvpz_im"));
    for (int k = 0; k < npar; ++k) {
      Real coskx(std::cos(kx * ppar->xp[k])), sinkx(std::sin(kx * ppar->xp[k]));
      Real coskz(std::cos(kz * ppar->yp[k])), sinkz(std::sin(kz * ppar->yp[k]));
      ppar->vpx[k] = vx0 + (dvpx_re * coskx - dvpx_im * sinkx) * coskz;
      ppar->vpz[k] = vy0 + (dvpy_re * coskx - dvpy_im * sinkx) * coskz;
      ppar->vpy[k] = -(dvpz_im * coskx + dvpz_re * sinkx) * sinkz;
    }
  }
}
