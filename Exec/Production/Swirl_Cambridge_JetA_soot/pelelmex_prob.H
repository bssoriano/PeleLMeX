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
#include "SootModel.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
set_Y_from_Phi(ProbParm const& prob_parm,
               const amrex::Real phi_in,
               amrex::Real Y[])
{
   auto eos = pele::physics::PhysicsType::eos();
   amrex::Real Xt[NUM_SPECIES] = {0.0};
   amrex::Real a = 0.0;
   // if ( prob_parm.fuelID == NXC7H16_ID ) {
   //   a = 11.0;
   // }
   Xt[prob_parm.oxidID] = 1.0 / ( 1.0 + phi_in / a + 0.79 / 0.21 );
   Xt[prob_parm.fuelID] = phi_in * Xt[prob_parm.oxidID] / a;
   Xt[prob_parm.bathID] = 1.0 - Xt[prob_parm.oxidID] - Xt[prob_parm.fuelID];

   eos.X2Y(Xt,Y);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
set_Y_from_Ksi(ProbParm const& prob_parm,
               const amrex::Real ksi_in,
               amrex::Real Y[])
{
   amrex::Real phi = ksi_in / ( 1.0 - amrex::min(ksi_in,0.999999) ) * (1.0 - prob_parm.Zst ) / prob_parm.Zst;

   set_Y_from_Phi(prob_parm,phi,Y);
}

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

    AMREX_D_TERM(,
                 const amrex::Real Lx = prob_hi[0] - prob_lo[0];,
                 const amrex::Real Ly = prob_hi[1] - prob_lo[1]);

    auto eos = pele::physics::PhysicsType::eos();
    amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};
    amrex::Real pert = 0.0;
    amrex::Real y1, y2, z1, z2;
    amrex::Real massfrac_eq[NUM_SPECIES] = {0.0};

    amrex::Real massfrac_fuel[NUM_SPECIES] = {0.0};
    amrex::Real massfrac_air[NUM_SPECIES] = {0.0};
    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    amrex::Real molefrac[NUM_SPECIES] = {0.0};
    amrex::Real molefrac_air[NUM_SPECIES] = {0.0};

    amrex::Real Do_swirler = prob_parm.Do_swirler;
    amrex::Real Di_swirler = prob_parm.Di_swirler;

    amrex::Real Ro_swirler = Do_swirler/2;
    amrex::Real Ri_swirler = Di_swirler/2;
    amrex::Real R = Ri_swirler/Ro_swirler;
    amrex::Real x0 = 0.0;
    amrex::Real y0 = 0.0;
    amrex::Real H = 0.01;
    amrex::Real D_coflow = Do_swirler;
    amrex::Real Ro_coflow = prob_hi[0] ;
    amrex::Real Ri_coflow = Ro_swirler;
    amrex::Real c1 = -0.0067;
    amrex::Real c2 = 3.39e-7;
    amrex::Real c3 = -21.081;
    amrex::Real A = 0.0096;
    amrex::Real B = -3.869e-11;
    amrex::Real F_theta = -7.1884e5;

    massfrac_air[prob_parm.oxidID] = 0.233;
    massfrac_air[prob_parm.bathID] = 0.767;

    molefrac_air[prob_parm.oxidID] = 0.21;
    molefrac_air[prob_parm.bathID] = 0.79;

    amrex::Real y0_1;
    amrex::Real z0_1 = 0.0;
    amrex::Real c1_1 = 50;
    amrex::Real c2_1 = 50;

    amrex::Real y0_2;
    amrex::Real z0_2 = 0.01;
    amrex::Real c1_2 = 150;
    amrex::Real c2_2 = 150;

    amrex::Real delta = 0.02;

    amrex::Real cut_off_height = 0.04;

    amrex::Real jet_1 = 0.5*(1-((std::tanh((std::sqrt(std::pow(x-x0,2)+std::pow(y-y0,2))-(Do_swirler/2))/(H))))*(tanh((std::sqrt(std::pow(x-x0,2)+std::pow(y-y0,2))+(Do_swirler/2))/(H))));

    amrex::Real jet_2 = 0.5*(1-std::tanh((z-cut_off_height)/H));

    amrex::Real jet_total = jet_1 * jet_2;

    z1 = (z - prob_parm.standoff - 0.5*dx[2] + pert)*prob_parm.amplification;
    z2 = (z - prob_parm.standoff + 0.5*dx[2] + pert)*prob_parm.amplification;

    // pele::physics::PMF::pmf(pmf_data, z1, z2, pmf_vals);
      
    state(i,j,k,TEMP) = (prob_parm.T0*jet_total) + ((1-jet_total)*(prob_parm.T_in));

    // state(i,j,k,TEMP) = prob_parm.T_in;

    for (int n = 0; n < NUM_SPECIES; ++n) {
        // molefrac[n] = (pmf_vals[3 + n]*jet_total)+((1-jet_total)*(molefrac_air[n]));
        molefrac[n] = molefrac_air[n];
    }

    eos.X2Y(molefrac, massfrac);

    state(i,j,k,VELX) = 0;
    state(i,j,k,VELY) = 0;
    state(i,j,k,VELZ) = prob_parm.U_b;

    amrex::Real rho_cgs, P_cgs;
    P_cgs = prob_parm.P_mean * 10.0;

    eos.PYT2R(P_cgs, massfrac, state(i,j,k,TEMP), rho_cgs);
    state(i,j,k,DENSITY) = rho_cgs * 1.0e3;            // CGS -> MKS conversion

    eos.TY2H(state(i,j,k,TEMP), massfrac, state(i,j,k,RHOH));
    state(i,j,k,RHOH) *= 1.0e-4 * state(i,j,k,DENSITY);   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
        state(i,j,k,FIRSTSPEC+n) = massfrac[n] * state(i,j,k,DENSITY);
    }

    for (int is = 0; is < NUM_SOOT_MOMENTS + 1; ++is) {
      state(i, j, k, FIRSTSOOT + is) = prob_parm.soot_vals[is];
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


  amrex::Real massfrac_air[NUM_SPECIES] = {0.0};
  massfrac_air[prob_parm.oxidID] = 0.233;
  massfrac_air[prob_parm.bathID] = 0.767;


  auto eos = pele::physics::PhysicsType::eos();

  if(idir == 2 and sgn == 1){
    amrex::GpuArray<amrex::Real, NUM_SPECIES+4> vel_from_file = {0.0};

    amrex::Real r = sqrt(pow((0.0-x[0]),2)+pow((0.0-x[1]),2));

    pele::physics::PMF::pmf(pmf_data, r, r, vel_from_file);

    amrex::Real vel_swirler_x = 0.0; 
    amrex::Real vel_swirler_y = 0.0; 
    amrex::Real vel_swirler_z = vel_from_file[5];

    double beta  = atan(x[1]/x[0]);
    double alpha = atan(x[0]/x[1]);

    vel_swirler_x = vel_from_file[3]*cos(beta) - vel_from_file[4]*cos(alpha);
    vel_swirler_y = vel_from_file[4]*cos(beta) + vel_from_file[3]*cos(alpha);


    if(x[1] < 0){
      vel_swirler_x = -vel_swirler_x;
    }

    if(x[0] < 0){
      vel_swirler_y = -vel_swirler_y;
    }

    // multiplied by 1.9 to increase bulk velocity to approx. 16 m/s
    // adding velocity components to the velocity fluctuation
    s_ext[VELX] += vel_swirler_x*1.9;
    s_ext[VELY] += vel_swirler_y*1.9;
    s_ext[VELZ] += vel_swirler_z*1.9;

    if(time < 13.5e-3){
      s_ext[TEMP] = 1800.;
    }
    else{
      s_ext[TEMP] = prob_parm.T_in;
    }


    amrex::Real rho_cgs, P_cgs, RhoH_temp;
    P_cgs = prob_parm.P_mean * 10.0;


    eos.PYT2R(P_cgs, massfrac_air, s_ext[TEMP], rho_cgs);
    s_ext[DENSITY] = rho_cgs * 1.0e3;

    eos.TY2H(s_ext[TEMP], massfrac_air, RhoH_temp);
    s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[FIRSTSPEC+n] = massfrac_air[n] * s_ext[DENSITY];
    }

    for (int is = 0; is < NUM_SOOT_MOMENTS + 1; ++is) {
      s_ext[FIRSTSOOT + is] = prob_parm.soot_vals[is];
    }
  }
  else{

    amrex::Real rho_cgs, P_cgs, RhoH_temp;
    P_cgs = prob_parm.P_mean * 10.0;

    s_ext[TEMP] = prob_parm.T_in;

    eos.PYT2R(P_cgs, massfrac_air, s_ext[TEMP], rho_cgs);
    s_ext[DENSITY] = rho_cgs * 1.0e3;

    eos.TY2H(s_ext[TEMP], massfrac_air, RhoH_temp);
    s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[FIRSTSPEC+n] = massfrac_air[n] * s_ext[DENSITY];
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
