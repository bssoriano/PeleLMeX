#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_SPACE.H>

#include <PeleLMeX_Index.H>
#include <pelelmex_prob_parm.H>
#include <PMF.H>
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
                 pele::physics::PMF::PmfData::DataContainer const *pmf_data)
{
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx      = geomdata.CellSize();

    AMREX_D_TERM(const amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];,
                 const amrex::Real y = prob_lo[1] + (j+0.5)*dx[1];,
                 const amrex::Real z = prob_lo[2] + (k+0.5)*dx[2];);

    AMREX_D_TERM(const amrex::Real Lx = prob_hi[0] - prob_lo[0];,
                 const amrex::Real Ly = prob_hi[1] - prob_lo[1];,
                 const amrex::Real Lz = prob_hi[2] - prob_lo[2];);

    auto eos = pele::physics::PhysicsType::eos();

    amrex::Real molefrac[NUM_SPECIES] = {0.0};
    amrex::Real massfrac[NUM_SPECIES] = {0.0};

    molefrac[O2_ID] = prob_parm.X_O2;
    molefrac[N2_ID] = prob_parm.X_N2;
    molefrac[H2O_ID] = prob_parm.X_H2O;
 
    //amrex::Real c = std::min(z/500.e-3,1.0);
       
    //state(i,j,k,TEMP) = c*prob_parm.T0+(1.0-c)*prob_parm.T_bottom;
    state(i,j,k,TEMP) = prob_parm.T0;

    eos.X2Y(molefrac, massfrac);

    state(i,j,k,VELX) = 0;
    state(i,j,k,VELY) = 0;
    state(i,j,k,VELZ) = prob_parm.vel;

    amrex::Real rho_cgs, P_cgs;
    P_cgs = prob_parm.P_mean * 10.0;

    eos.PYT2R(P_cgs, massfrac, state(i,j,k,TEMP), rho_cgs);
    state(i,j,k,DENSITY) = rho_cgs * 1.0e3;            // CGS -> MKS conversion

    eos.TY2H(state(i,j,k,TEMP), massfrac, state(i,j,k,RHOH));
    state(i,j,k,RHOH) *= 1.0e-4 * state(i,j,k,DENSITY);   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
        state(i,j,k,FIRSTSPEC+n) = massfrac[n] * state(i,j,k,DENSITY);
    }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
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
  const amrex::Real* prob_hi = geomdata.ProbHi();


  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  amrex::Real molefrac[NUM_SPECIES] = {0.0};

  molefrac[O2_ID]  = prob_parm.X_O2;
  molefrac[N2_ID]  = prob_parm.X_N2;
  molefrac[H2O_ID] = prob_parm.X_H2O;

  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real temp = prob_parm.T0;
  
  s_ext[VELZ] = prob_parm.Vin;

  if(idir == 2){
    //if(sgn == 1){
      amrex::Real y_center = 0.0e-3;
      amrex::Real x_center = 0.0e-3;
      amrex::Real r = sqrt(pow((y_center-x[1]),2)+pow((x_center-x[0]),2));
      if(r < prob_parm.inj_diam){
        s_ext[VELZ] = prob_parm.Vin;
        temp = prob_parm.T_in;
      }
    //} 
    
    eos.X2Y(molefrac, massfrac);

    s_ext[VELY] = 0.0;
    s_ext[VELX] = 0.0;

    s_ext[TEMP] = temp;

    amrex::Real rho_cgs, P_cgs, RhoH_temp;
    P_cgs = prob_parm.P_mean * 10.0;

    eos.PYT2R(P_cgs, massfrac, s_ext[TEMP], rho_cgs);
    s_ext[DENSITY] = rho_cgs * 1.0e3;

    eos.TY2H(s_ext[TEMP], massfrac, RhoH_temp);
    s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[FIRSTSPEC+n] = massfrac[n] * s_ext[DENSITY];
    }
  }

}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
zero_visc (int /*i*/, int /*j*/, int /*k*/,
           amrex::Array4<amrex::Real> const& /*beta*/,
           amrex::GeometryData const& /*geomdata*/,
           amrex::Box const& /*dBox*/,
           const int  /*dir*/,
           const int  /*state_comp*/,
           const int  /*nComp*/)
{
}
#endif
