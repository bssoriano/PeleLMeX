#ifndef PELELM_PROB_H
#define PELELM_PROB_H

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_SPACE.H>

#include <PeleLMeX_ProblemSpecificFunctions.H>
#include <PeleLMeX_Index.H>
#include <PMF.H>
#include <PMFData.H>
#include <PelePhysics.H>
#include <Utilities.H>

namespace m2c = pele::physics::utilities::mks2cgs;
namespace c2m = pele::physics::utilities::cgs2mks;

struct MyProbParm : public ProbParmDefault
{

  amrex::Real T_mean = 298.0;
  amrex::Real P_mean = 101325.0;
  amrex::Real meanFlowMag = 0.0;
  amrex::Real lowFlowMag = 3.0;
  amrex::Real jet_radius = 0.0;
  int meanFlowDir = 1;
  int Flame_IC = 0;
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> jet_loc = {{0.0}};
  amrex::Real jet_vel = 0.0;
  amrex::Real jet_T = 0.0;
  amrex::Real jet_Yfuel = 0.0;
  amrex::Real Ktr  = 0.0;
  amrex::Real t_tr = 0.0;
  amrex::Real start_time_injProd = 0.0;
  amrex::Real duration_injProd = 0.0;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y_prod = {{0.0}};

  pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type> const*
    eosparm;
  
  int Zvar_ID = 1;
  int Z_ID = 2;
  int C_ID = 0;
  
  //int Z_ID = 1;
  //int C_ID = 0;

};

struct MyProblemSpecificFunctions : public DefaultProblemSpecificFunctions
{
  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static void
  initdata(
    int i,
    int j,
    int k,
    int is_incompressible,
    amrex::Array4<amrex::Real> const& state,
    amrex::Array4<amrex::Real> const& /*aux*/,
    amrex::GeometryData const& geomdata,
    MyProbParm const& prob_parm,
    pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
  {
    auto eos = pele::physics::PhysicsType::eos(prob_parm.eosparm);
    
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* dx = geomdata.CellSize();
  
    const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
    const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
  
  
    switch (prob_parm.meanFlowDir) {
    case 1:
      AMREX_D_TERM(state(i, j, k, VELX) = prob_parm.meanFlowMag;, state(i, j, k, VELY) = 0.0;
                   , state(i, j, k, VELZ) = 0.0);
      break;
    case 2:
      AMREX_D_TERM(state(i, j, k, VELY) = prob_parm.meanFlowMag;, state(i, j, k, VELX) = 0.0 ;
                   , state(i, j, k, VELZ) = 0.0);
      break;
    case 3:
      AMREX_D_TERM(state(i, j, k, VELZ) = prob_parm.meanFlowMag;, state(i, j, k, VELX) = 0.0;
                   , state(i, j, k, VELY) = 0.0);
      break;
    default:
      break;
    }
  
    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    massfrac[prob_parm.Zvar_ID] = 0.0;
    massfrac[prob_parm.C_ID] = 0.0;
    massfrac[prob_parm.Z_ID] = 0.0;
  
    state(i, j, k, TEMP) = prob_parm.T_mean;
  
    amrex::Real rho_temp, tmp1, tmp2;
    eos.PYT2R(tmp1, massfrac, tmp2, rho_temp);
    state(i, j, k, DENSITY) = c2m::Rho(rho_temp);
    state(i, j, k, RHOH) = 0.0; // No RhoH transport in manifold model
    for (int n = 0; n < NUM_SPECIES - 1; n++) {
      state(i, j, k, FIRSTSPEC + n) = massfrac[n] * state(i, j, k, DENSITY);
    }
    state(i, j, k, FIRSTSPEC + NUM_SPECIES - 1) = state(i, j, k, DENSITY);  
  }
  
  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static void bcnormal(
    const amrex::Real x[AMREX_SPACEDIM],
    const amrex::CellData<amrex::Real>& /*s_in*/,
    amrex::Real s_ext[NVAR],
    const int /*idir*/,
    const int /*sgn*/,
    const amrex::Real time,
    amrex::GeometryData const& /*geomdata*/,
    MyProbParm const& prob_parm,
    pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
  {
  
    amrex::Real r = sqrt(pow((prob_parm.jet_loc[1]-x[1]),2)+pow((prob_parm.jet_loc[2]-x[2]),2));
    
    amrex::Real massfrac[NUM_SPECIES] = {0.0};
  
    auto eos = pele::physics::PhysicsType::eos(prob_parm.eosparm);
  
    switch (prob_parm.meanFlowDir) {
    case 1:
      AMREX_D_TERM(s_ext[VELX] = prob_parm.meanFlowMag;, s_ext[VELY] = 0.0;
                   , s_ext[VELZ] = 0.0);
      break;
    case 2:
      AMREX_D_TERM(s_ext[VELY] = prob_parm.meanFlowMag;,if (r <= prob_parm.jet_radius) {s_ext[VELX] = prob_parm.jet_vel;} ;
                   , s_ext[VELZ] = 0.0);
      break;
    case 3:
      AMREX_D_TERM(s_ext[VELZ] = prob_parm.meanFlowMag;, s_ext[VELX] = 0.0;
                   , s_ext[VELY] = 0.0);
      break;
    default:
      break;
    }
    
    //if(time < prob_parm.t_tr){
    //  s_ext[VELY] = 0.0;
    //}
    massfrac[prob_parm.Z_ID] = 0.0;
    massfrac[prob_parm.Zvar_ID] = 0.0;
    massfrac[prob_parm.C_ID] = 0.0;
  
    s_ext[TEMP] = prob_parm.T_mean;
  
    amrex::Real P_dummy, rho_cgs;
    eos.PYT2R(P_dummy, massfrac, s_ext[TEMP], rho_cgs);
    s_ext[DENSITY] = rho_cgs * 1.0e3;
    s_ext[RHOH] = 0.0;
    for (int n = 0; n < NUM_SPECIES - 1; n++) {
          s_ext[FIRSTSPEC + n] = massfrac[n] * s_ext[DENSITY];
    }
    s_ext[FIRSTSPEC + NUM_SPECIES - 1] = s_ext[DENSITY];
  
  }
//};
  
      AMREX_GPU_DEVICE
      AMREX_FORCE_INLINE
      static void bctype_eb(
        const amrex::Real xEBface[AMREX_SPACEDIM],
        int& EBflagType,
        amrex::Real& EBfaceFrac,
        amrex::GeometryData const& /*geomdata*/,
        MyProbParm const& prob_parm)

      {
  #if (AMREX_SPACEDIM == 2)
          amrex::Real r = sqrt(pow((prob_parm.jet_loc[1]-xEBface[1]),2));
  #elif (AMREX_SPACEDIM == 3)
          amrex::Real r = sqrt(pow((prob_parm.jet_loc[1]-xEBface[1]),2)+pow((prob_parm.jet_loc[2]-xEBface[2]),2));
  #endif
  	if (r <= prob_parm.jet_radius) {
              EBflagType = pelelmex::BCTypeEB::inflow;
	      EBfaceFrac = 1.0;
          }
	else{
	      EBflagType = pelelmex::BCTypeEB::wall_adiab;
	      EBfaceFrac = 0.0;
	}
          
      }
    
      AMREX_GPU_DEVICE
      AMREX_FORCE_INLINE
      static void bcnormal_eb(
      const amrex::Real xEBface[AMREX_SPACEDIM],
      const amrex::Real* /*xEBfaceNormal[AMREX_SPACEDIM]*/,
      const amrex::CellData<const amrex::Real>& /*s_in*/,
      amrex::Real s_ext[NVAR],
      const amrex::Real time,
      amrex::GeometryData const& /*geomdata*/,
      MyProbParm const& prob_parm)

      {
  #if (AMREX_SPACEDIM == 2)
          amrex::Real r = sqrt(pow((prob_parm.jet_loc[1]-xEBface[1]),2));
  #elif (AMREX_SPACEDIM == 3)
          amrex::Real r = sqrt(pow((prob_parm.jet_loc[1]-xEBface[1]),2)+pow((prob_parm.jet_loc[2]-xEBface[2]),2));
  #endif
          // Isothermal wall
          //s_ext[TEMP] = prob_parm.Twall;
  
          // Inflow
          if (r <= prob_parm.jet_radius) {
  
              auto eos = pele::physics::PhysicsType::eos(prob_parm.eosparm);
              s_ext[VELX] = prob_parm.jet_vel;
              s_ext[VELY] = 0.0;
              s_ext[VELZ] = 0.0;
              
  	      amrex::Real massfrac[NUM_SPECIES] = {0.0};
              

              massfrac[prob_parm.Z_ID] = prob_parm.jet_Yfuel;
	      massfrac[prob_parm.Zvar_ID] = 0.0;
              
	      if(time > prob_parm.start_time_injProd and time < prob_parm.start_time_injProd+prob_parm.duration_injProd){
                  massfrac[prob_parm.C_ID] = 1.0;
	      }
	      else{
                  massfrac[prob_parm.C_ID] = 0.0;
	      }

              amrex::Real P_cgs = prob_parm.P_mean * 10.0;
              
              //s_ext[TEMP] = prob_parm.jet_T;
  	    //amrex::Real fb = 0.5*(tanh((time-prob_parm.t_tr)/prob_parm.Ktr)+1);
  	    //amrex::Real T = (1.-fb)*2266.64 + fb*prob_parm.jet_T;
  	    s_ext[TEMP] = prob_parm.jet_T;
  
  	    //for (int n = 0; n < NUM_SPECIES; n++) {
              //    massfrac[n] = fb * massfrac[n] + (1.0 - fb) * prob_parm.Y_prod[n];
              //}
              amrex::Real P_dummy, rho_cgs;
              eos.PYT2R(P_dummy, massfrac, s_ext[TEMP], rho_cgs);
              s_ext[DENSITY] = rho_cgs * 1.0e3;
              s_ext[RHOH] = 0.0;
              for (int n = 0; n < NUM_SPECIES - 1; n++) {
                    s_ext[FIRSTSPEC + n] = massfrac[n] * s_ext[DENSITY];
              }
              s_ext[FIRSTSPEC + NUM_SPECIES - 1] = s_ext[DENSITY];
  
          }
      }
};


using ProblemSpecificFunctions = MyProblemSpecificFunctions;
using ProbParm = MyProbParm;

#endif
