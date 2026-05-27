#ifndef PELELM_PROB_H
#define PELELM_PROB_H

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_SPACE.H>

#ifdef AMREX_USE_EB
#include <AMReX_ParmParse.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#endif

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

  // Inflow duct x/z extents (read from eb2.inflow_* in input file)
  amrex::Real inflow_x_lo = 0.0;
  amrex::Real inflow_x_hi = 0.0;
  amrex::Real inflow_z_lo = 0.0;
  amrex::Real inflow_z_hi = 0.0;
  
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
  
    if (!is_incompressible) {
      amrex::Real massfrac[NUM_SPECIES] = {0.0};
      massfrac[prob_parm.Zvar_ID] = 0.0;
      massfrac[prob_parm.C_ID] = 0.0;
      massfrac[prob_parm.Z_ID] = 0.0;

      state(i, j, k, TEMP) = prob_parm.T_mean;

      amrex::Real rho_temp, tmp1, tmp2;
      eos.PYT2R(tmp1, massfrac, tmp2, rho_temp);
      state(i, j, k, DENSITY) = c2m::Rho(rho_temp);
      state(i, j, k, RHOH) = 0.0;
      for (int n = 0; n < NUM_SPECIES - 1; n++) {
        state(i, j, k, FIRSTSPEC + n) = massfrac[n] * state(i, j, k, DENSITY);
      }
      state(i, j, k, FIRSTSPEC + NUM_SPECIES - 1) = state(i, j, k, DENSITY);
    }
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
  
    const bool in_inflow_region = (x[0] >= prob_parm.inflow_x_lo &&
                                   x[0] <= prob_parm.inflow_x_hi &&
                                   x[2] >= prob_parm.inflow_z_lo &&
                                   x[2] <= prob_parm.inflow_z_hi);

    switch (prob_parm.meanFlowDir) {
    case 1:
      AMREX_D_TERM(s_ext[VELX] = prob_parm.meanFlowMag;, s_ext[VELY] = 0.0;
                   , s_ext[VELZ] = 0.0);
      break;
    case 2:
      if (in_inflow_region) {
        AMREX_D_TERM(s_ext[VELY] = prob_parm.meanFlowMag;,
                     if (r <= prob_parm.jet_radius) {s_ext[VELX] = prob_parm.jet_vel;};,
                     s_ext[VELZ] = 0.0);
      } else {
        AMREX_D_TERM(s_ext[VELX] = 0.0;, s_ext[VELY] = 0.0;, s_ext[VELZ] = 0.0);
      }
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

#ifdef AMREX_USE_EB
  static void EBUserDefined(
    const amrex::Geometry& geom,
    const int required_coarsening_level,
    const int max_coarsening_level)
  {
    amrex::ParmParse pp("eb2");

    amrex::Real radius;
    pp.get("cylinder_radius", radius);

    amrex::Real height;
    pp.get("cylinder_height", height);

    int direction;
    pp.get("cylinder_direction", direction);

    amrex::Vector<amrex::Real> center_v;
    pp.getarr("cylinder_center", center_v);
    amrex::RealArray center = {AMREX_D_DECL(center_v[0], center_v[1], center_v[2])};

    // fluid_inside = false: cylinder is a solid obstacle (fluid outside)
    bool fluid_inside = false;
    amrex::EB2::CylinderIF cyl(radius, height, direction, center, fluid_inside);

    // Thin walls separating the inflow region from the zero-velocity buffer at ylo.
    // Boundaries of the inflow region (must match bcnormal constants).
    amrex::Real inflow_x_lo, inflow_x_hi, inflow_z_lo, inflow_z_hi;
    pp.get("inflow_x_lo", inflow_x_lo);
    pp.get("inflow_x_hi", inflow_x_hi);
    pp.get("inflow_z_lo", inflow_z_lo);
    pp.get("inflow_z_hi", inflow_z_hi);

    // Wall thickness (into the buffer region) and depth (into the domain, +y).
    amrex::Real wall_thickness, wall_depth;
    pp.get("thin_wall_thickness", wall_thickness);
    pp.get("thin_wall_depth",     wall_depth);

    // Extend walls slightly below ylo to ensure full coverage at the boundary face.
    const amrex::Real y_lo_wall = geom.ProbLo()[1] - 0.01;
    const amrex::Real y_hi_wall = geom.ProbLo()[1] + wall_depth;

    // Wall at x = inflow_x_lo: thickness extends into x < inflow_x_lo buffer.
    amrex::EB2::BoxIF wall_xlo(
      {AMREX_D_DECL(inflow_x_lo - wall_thickness, y_lo_wall, inflow_z_lo)},
      {AMREX_D_DECL(inflow_x_lo,                  y_hi_wall, inflow_z_hi)},
      false);

    // Wall at x = inflow_x_hi: thickness extends into x > inflow_x_hi buffer.
    amrex::EB2::BoxIF wall_xhi(
      {AMREX_D_DECL(inflow_x_hi,                  y_lo_wall, inflow_z_lo)},
      {AMREX_D_DECL(inflow_x_hi + wall_thickness,  y_hi_wall, inflow_z_hi)},
      false);

    // Wall at z = inflow_z_lo: thickness extends into z < inflow_z_lo buffer.
    // x span extended to cover corners where x-walls sit.
    amrex::EB2::BoxIF wall_zlo(
      {AMREX_D_DECL(inflow_x_lo - wall_thickness, y_lo_wall, inflow_z_lo - wall_thickness)},
      {AMREX_D_DECL(inflow_x_hi + wall_thickness,  y_hi_wall, inflow_z_lo)},
      false);

    // Wall at z = inflow_z_hi: thickness extends into z > inflow_z_hi buffer.
    amrex::EB2::BoxIF wall_zhi(
      {AMREX_D_DECL(inflow_x_lo - wall_thickness, y_lo_wall, inflow_z_hi)},
      {AMREX_D_DECL(inflow_x_hi + wall_thickness,  y_hi_wall, inflow_z_hi + wall_thickness)},
      false);

    // Union all solid objects.
    auto union_walls_x = amrex::EB2::makeUnion(wall_xlo,      wall_xhi);
    auto union_walls_z = amrex::EB2::makeUnion(wall_zlo,      wall_zhi);
    auto union_walls   = amrex::EB2::makeUnion(union_walls_x, union_walls_z);
    auto all_solid     = amrex::EB2::makeUnion(cyl,           union_walls);

    auto gshop = amrex::EB2::makeShop(all_solid);
    amrex::EB2::Build(gshop, geom, required_coarsening_level, max_coarsening_level);
  }
#endif

};


using ProblemSpecificFunctions = MyProblemSpecificFunctions;
using ProbParm = MyProbParm;

#endif
