#ifndef PELELM_PROB_H
#define PELELM_PROB_H

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_SPACE.H>

#include <PeleLMeX_Index.H>
#include <pelelmex_prob_parm.H>
#include <PMF.H>
#include <PMFData.H>
#include <PelePhysics.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pelelmex_initdata(
  int i,
  int j,
  int k,
  int /*is_incompressible*/,
  amrex::Array4<amrex::Real> const& state,
  amrex::Array4<amrex::Real> const& /*aux*/,
  amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm,
  pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
{

    const amrex::Real *prob_lo = geomdata.ProbLo();
    const amrex::Real *prob_hi = geomdata.ProbHi();
    const amrex::Real *dx = geomdata.CellSize();

    AMREX_D_TERM(const amrex::Real Lxm = .5*(prob_hi[0] + prob_lo[0]);
                 , const amrex::Real Lym = .5*(prob_hi[1] + prob_lo[1]);
                 , const amrex::Real Lzm = .5*(prob_hi[2] + prob_lo[2]));

    amrex::Real x[3] = {
      prob_lo[0] + static_cast<amrex::Real>(i + 0.5) * dx[0],
      prob_lo[1] + static_cast<amrex::Real>(j + 0.5) * dx[1],
      prob_lo[2] + static_cast<amrex::Real>(k + 0.5) * dx[2]};

    amrex::Real D = 2.e-3;
    amrex::Real radius = std::sqrt(std::pow(x[0]-Lxm,2) + std::pow(x[1]-Lym,2) + std::pow(x[2]-Lzm,2));
    amrex::Real eta = 0.5*(1.0 - std::tanh(2.0*(radius-D*0.5)/(0.3*D*0.5)));

  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};


  massfrac[N2_ID]  = prob_parm.Y_N2;
  massfrac[O2_ID]  = prob_parm.Y_O2;
  massfrac[CO2_ID] = prob_parm.Y_CO2;
  massfrac[H2O_ID] = prob_parm.Y_H2O;

 // massfrac[N2_ID] = massfrac[N2_ID] - eta;
 // massfrac[POSF10325_ID] = eta;

  state(i, j, k, TEMP) = prob_parm.T0;

  state(i, j, k, VELX) = prob_parm.vel;
#if (AMREX_SPACEDIM == 2)
  state(i, j, k, VELY) = 0.0;
#elif (AMREX_SPACEDIM == 3)
  state(i, j, k, VELY) = 0.0;
  state(i, j, k, VELZ) = 0.0;
#endif

  amrex::Real P_cgs = prob_parm.P_mean * 10.0;

  // Density
  amrex::Real rho_cgs = 0.0;
  eos.PYT2R(P_cgs, massfrac, state(i, j, k, TEMP), rho_cgs);
  state(i, j, k, DENSITY) = rho_cgs * 1.0e3;

  // Enthalpy
  amrex::Real h_cgs = 0.0;
  eos.TY2H(state(i, j, k, TEMP), massfrac, h_cgs);
  state(i, j, k, RHOH) = h_cgs * 1.0e-4 * state(i, j, k, DENSITY);

  // Species mass
  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, FIRSTSPEC + n) = massfrac[n] * state(i, j, k, DENSITY);
  }

#ifdef PELELM_USE_MF
    state(i,j,k,FIRSTMFVAR) = 0.0; //eta*state(i, j, k, DENSITY);
    state(i,j,k,FIRSTMFVAR+1) = 0.;
#endif
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  const int /*m_nAux*/,
  amrex::Real s_ext[NVAR],
  const int /*idir*/,
  const int /*sgn*/,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParm const& prob_parm,
  pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
{
  AMREX_D_TERM(s_ext[VELX] = prob_parm.vel;, s_ext[VELY] = 0.;
               , s_ext[VELZ] = 0.;);
  // amrex::GpuArray<amrex::Real, NUM_SPECIES> molefrac = {{0.0}};
  // amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {{0.0}};
  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  
  auto eos = pele::physics::PhysicsType::eos();
  
  massfrac[N2_ID]  = prob_parm.Y_N2;
  massfrac[O2_ID]  = prob_parm.Y_O2;
  massfrac[CO2_ID] = prob_parm.Y_CO2;
  massfrac[H2O_ID] = prob_parm.Y_H2O;
  
  //eos.X2Y(molefrac,massfrac);

  s_ext[TEMP] = prob_parm.T0;
  amrex::Real rho_cgs, P_cgs;
  P_cgs = prob_parm.P_mean * 10.;

  //auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2R(P_cgs, massfrac, prob_parm.T0, rho_cgs);
  amrex::Real rho_mks = rho_cgs * 1.0e3; // CGS -> MKS conversion
  s_ext[DENSITY] = rho_mks;
  amrex::Real rhoH;
  eos.TY2H(prob_parm.T0, massfrac, rhoH);
  s_ext[RHOH] = rhoH * 1.0e-4 * rho_mks;
  for (int n = 0; n < NUM_SPECIES; n++) {
    s_ext[FIRSTSPEC + n] = massfrac[n] * rho_mks;
  }
#ifdef PELELM_USE_MF
    s_ext[FIRSTMFVAR] = 0.;
    s_ext[FIRSTMFVAR+1] = 0.;
#endif
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
zero_visc(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& beta,
  amrex::GeometryData const& geomdata,
  amrex::Box const& domainBox,
  const int dir,
  const int beta_comp,
  const int nComp)
{
  amrex::ignore_unused(
    i, j, k, beta, geomdata, domainBox, dir, beta_comp, nComp);
  // We treat species when beta_comp == 0 and nComp == NUM_SPECIES
  // otherwise this routine could be called for other face diffusivity (Temp,
  // velocity, ...)
}
#endif
