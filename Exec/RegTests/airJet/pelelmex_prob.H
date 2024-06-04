#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_SPACE.H>

#include <PeleLMeX_Index.H>
#include <pelelmex_prob_parm.H>
#include <PMFData.H>
#include <PelePhysics.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pelelmex_initdata (int i, int j, int k,
                 int is_incompressible,
                 amrex::Array4<amrex::Real> const& state,
                 amrex::Array4<amrex::Real> const& aux,
                 amrex::GeometryData const& geomdata,
                 ProbParm const& prob_parm,
                 pele::physics::PMF::PmfData::DataContainer const * /*pmf_data*/)
{
    auto eos = pele::physics::PhysicsType::eos();

    amrex::Real massfrac[NUM_SPECIES] = {0.0};

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

    amrex::Real radius = std::sqrt(std::pow(x[0]-0.0,2) + std::pow(x[1]-0.0,2) + std::pow(x[2]-Lzm,2));
    amrex::Real eta = 0.5*(1.0 - std::tanh(2.0*(radius-prob_parm.D*0.5)/(0.3*prob_parm.D*0.5)));

    amrex::Real enthalpy_cgs = (1.0-eta) * prob_parm.H_fuel + eta*prob_parm.H_ox;
    if(radius < prob_parm.D*0.5){
       enthalpy_cgs = prob_parm.H_ox;
    }
    else{
       enthalpy_cgs = prob_parm.H_fuel;
    }
    
    for (int n = 0; n < NUM_SPECIES; n++)
    {
        massfrac[n] = prob_parm.Y_ox[n];
    }

    amrex::Real T;
    eos.HY2T(enthalpy_cgs, massfrac, T);
    state(i,j,k,TEMP) = T;
    
    state(i,j,k,VELX) = 0.1;
    state(i,j,k,VELY) = 0.0;
    state(i,j,k,VELZ) = 0.0;

    amrex::Real rho_cgs, P_cgs;
    P_cgs = prob_parm.P_mean * 10.0;

    eos.PYT2R(P_cgs, massfrac, state(i,j,k,TEMP), rho_cgs);
    state(i,j,k,DENSITY) = rho_cgs * 1.0e3;            // CGS -> MKS conversion

    state(i,j,k,RHOH) = 1.0e-4 * enthalpy_cgs * state(i,j,k,DENSITY);   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      state(i,j,k,FIRSTSPEC+n) = massfrac[n] * state(i,j,k,DENSITY);
    }

#ifdef PELELM_USE_MF
    if(radius < prob_parm.D*0.5){
       state(i,j,k,FIRSTMFVAR) = state(i,j,k,DENSITY); 
    }
    else{
       state(i,j,k,FIRSTMFVAR) = 0.0; 
    }
    state(i,j,k,FIRSTMFVAR) = eta * state(i,j,k,DENSITY); 
    state(i,j,k,FIRSTMFVAR+1) = 0.; 
    state(i,j,k,FIRSTMFVAR+2) = 0.;
#endif
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
time_func(amrex::Real time, amrex::Real start, amrex::Real duration, amrex::Real tau)
{
  amrex::Real eta_time;
  eta_time = 0.5 * (std::tanh((time - start)/(duration*tau))*
    std::tanh((-time + start + duration)/(duration*tau))+1.0);
  eta_time = std::min(std::max(eta_time,0.0),1.0);
  return eta_time;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real *x, /*x[AMREX_SPACEDIM]*/
  const int m_nAux,
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm,
  pele::physics::PMF::PmfData::DataContainer const *pmf_data)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real radius = {0.0};
  amrex::Real enthalpy;
  int ierr;

  if (sgn == 1 && idir == 2) {

    // Jet function
    radius = std::sqrt(std::pow(x[0]-prob_parm.center_xy[0],2) + std::pow(x[1]-prob_parm.center_xy[1],2));
    amrex::Real eta = 0.5*(1.0 - std::tanh(2.0*(radius-prob_parm.D*0.5)/(0.1*prob_parm.D*0.5)));
    eta *= time_func(time,prob_parm.inj_start,prob_parm.inj_dur,prob_parm.tau);

    s_ext[VELX] = 0.0;
    s_ext[VELY] = 0.0;
    s_ext[VELZ] = eta*prob_parm.v_in;
    
    for (int n = 0; n < NUM_SPECIES; n++)
    {
      massfrac[n] = prob_parm.Y_ox[n];
    }
    enthalpy = (1.-eta) * prob_parm.H_fuel + eta*prob_parm.H_ox;

    GET_T_GIVEN_HY(enthalpy,massfrac,s_ext[TEMP],ierr);

    amrex::Real rho_cgs, P_cgs, RhoH_temp;
    P_cgs = prob_parm.P_mean * 10.0;

    eos.PYT2R(P_cgs, massfrac, s_ext[TEMP], rho_cgs);
    s_ext[DENSITY] = rho_cgs * 1.0e3;

    eos.TY2H(s_ext[TEMP], massfrac, RhoH_temp);
    s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[FIRSTSPEC+n] = massfrac[n] * s_ext[DENSITY];
    }

#ifdef PELELM_USE_MF
    s_ext[FIRSTMFVAR] = 0.;
    s_ext[FIRSTMFVAR+1] = eta * s_ext[DENSITY];
    s_ext[FIRSTMFVAR+2] = 0.;
#endif
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
zero_visc (int i, int j, int k,
           amrex::Array4<amrex::Real> const& beta,
           amrex::GeometryData const& geomdata,
           amrex::Box const& domainBox,
           const int  dir,
           const int  beta_comp,
           const int  nComp)
{
     amrex::ignore_unused(i,j,k,beta,geomdata,domainBox,dir,beta_comp,nComp);
    // We treat species when beta_comp == 0 and nComp == NUM_SPECIES
    // otherwise this routine could be called for other face diffusivity (Temp, velocity, ...)
}
#endif
