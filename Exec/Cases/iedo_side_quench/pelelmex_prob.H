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

// -----------------------------------------------------------
// Search for the closest index in an array to a given value
// using the bisection technique.
// INPUTS/OUTPUTS:
// xtable(0:n-1) => array to search in (ascending order)
// n             => array size
// x             => x location
// idxlo        <=> output st. xtable(idxlo) <= x < xtable(idxlo+1)
// -----------------------------------------------------------
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
locate(const amrex::Real* xtable, const int n, const amrex::Real& x, int& idxlo)
{
  // If x is out of bounds, return boundary index
  if (x >= xtable[n - 1]) {
    idxlo = n - 1;
    return;
  }
  if (x <= xtable[0]) {
    idxlo = 0;
    return;
  }

  // Do the bisection
  idxlo = 0;
  int idxhi = n - 1;
  bool notdone = true;
  while (notdone) {
    if (idxhi - idxlo <= 1) {
      notdone = false;
    } else {
      const int idxmid = (idxhi + idxlo) / 2;
      if (x >= xtable[idxmid]) {
        idxlo = idxmid;
      } else {
        idxhi = idxmid;
      }
    }
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pelelmex_initdata(
  int i,
  int j,
  int k,
  int /*is_incompressible*/,
  amrex::Array4<amrex::Real> const& state,
  amrex::Array4<amrex::Real> const& aux,
  amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm,
  pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
{

  // ------ Initialization of flow field -----
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();

  amrex::Real x[3] = {
    prob_lo[0] + static_cast<amrex::Real>(i + 0.5) * dx[0],
    prob_lo[1] + static_cast<amrex::Real>(j + 0.5) * dx[1],
    prob_lo[2] + static_cast<amrex::Real>(k + 0.5) * dx[2]};

  // -----------------------------------------
  state(i, j, k, VELX) = 0.0;
  state(i, j, k, VELY) = 0.0;
  state(i, j, k, VELZ) = 0.0;

  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real massfrac[NUM_SPECIES] = {0.0};

  state(i, j, k, TEMP) = prob_parm.T_mean;
  
  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] = prob_parm.Y_chamber[n];
  }

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

  // Auxiliary scalars
  for (int n = 0; n < aux.nComp(); n++) {
    aux(i, j, k, n) = 0.0; 
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  const int /*m_nAux*/,
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& /*geomdata*/,
  ProbParm const& prob_parm,
  pele::physics::PMF::PmfData::DataContainer const* pmf_data)
{
  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  amrex::Real X_air[NUM_SPECIES] = {0.0};
  amrex::Real Y_air[NUM_SPECIES] = {0.0};
  amrex::Real X_jet[NUM_SPECIES] = {0.0};
  amrex::Real Y_jet[NUM_SPECIES] = {0.0};
  amrex::Real eta = 0.0;
  constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};

  // Simple, air into air injection
  X_air[O2_ID] = 0.21;
  X_air[N2_ID] = 1.0 - X_air[O2_ID];

  eos.X2Y(X_air, Y_air);
  //eos.X2Y(X_jet, Y_jet);


#if (AMREX_SPACEDIM == 2)
  if (idir == 1 and sgn == 1) {
      amrex::Real r = sqrt(pow((0.0-x[0]),2));
#elif (AMREX_SPACEDIM == 3)
  if (idir == 0 and sgn == 1) {
      amrex::Real r = sqrt(pow((prob_parm.jet_center[2]-x[2]),2)/pow(std::cos(prob_parm.jet_angle),2)+pow((prob_parm.jet_center[0]-x[1]),2));
#endif

      amrex::Real N = 20.;
      amrex::Real radial_factor = 0.0;
      //streamwise velocity as a function of radial location
      if(r <= prob_parm.jet_rad){
          radial_factor = (1-pow(r/prob_parm.jet_rad,N));
      }
      
      //amrex::Real vel_jet = prob_parm.V_jet*radial_factor;
      pele::physics::PMF::pmf(pmf_data, time, time, pmf_vals);

      amrex::Real vel_jet = 0.5*pmf_vals[0]*radial_factor;

      //temporal velocity function
      amrex::Real fu = 0.5*(tanh((time-prob_parm.t_si)/prob_parm.Ku)-tanh((time-prob_parm.t_ei)/prob_parm.Ku));
      vel_jet *= fu;


      //temporal burnt gas function
      amrex::Real fb = 0.5*(tanh((time-prob_parm.t_tr)/prob_parm.Ktr)+1);
      amrex::Real Tb = prob_parm.T_mean+(prob_parm.T_jet-prob_parm.T_mean)*radial_factor;
      // This is to enforce T_mean away from the nozzle region
      amrex::Real T = prob_parm.T_mean;
      if(r <= prob_parm.jet_rad){
          T = fb*Tb + (1-fb)*prob_parm.T_mean;
	  //T = fu*prob_parm.T_jet+(1-fu)*prob_parm.T_mean;

      }
      
      if(r <= prob_parm.jet_rad){
          for (int n = 0; n < NUM_SPECIES; n++) {
            massfrac[n] = fb * prob_parm.Y_jet[n] + (1.0 - fb) * prob_parm.Y_prechamber[n];
          }
      }
      else{
	  for (int n = 0; n < NUM_SPECIES; n++) {
            massfrac[n] = prob_parm.Y_chamber[n];
	  }
      }

      amrex::Real rho_cgs, P_cgs, RhoH_temp;
      P_cgs = prob_parm.P_mean * 10.0;

      s_ext[TEMP] = T;

      eos.PYT2R(P_cgs, massfrac, s_ext[TEMP], rho_cgs);
      s_ext[DENSITY] = rho_cgs * 1.0e3;

      eos.TY2H(s_ext[TEMP], massfrac, RhoH_temp);
      s_ext[RHOH] =
        RhoH_temp * 1.0e-4 * s_ext[DENSITY]; // CGS -> MKS conversion

      for (int n = 0; n < NUM_SPECIES; n++) {
        s_ext[FIRSTSPEC + n] = massfrac[n] * s_ext[DENSITY];
      }

      // Velocity
      if(r <= prob_parm.jet_rad){ //and time < (prob_parm.t_ei+prob_parm.Ku)){
#if (AMREX_SPACEDIM == 2)
          s_ext[VELX] += 0.0;
          s_ext[VELY] += std::sin(prob_parm.jet_angle)*vel_jet;
#elif (AMREX_SPACEDIM == 3)
	  //Scale vel. fluctuation based on bulk jet velocity
	  amrex::Real scaled_velx_fluct = s_ext[VELX]*vel_jet/100.; //100. is the maximum jet velocity
	  amrex::Real scaled_vely_fluct = s_ext[VELY]*vel_jet/100.;
	  amrex::Real scaled_velz_fluct = s_ext[VELZ]*vel_jet/100.;

          s_ext[VELX] = scaled_velx_fluct + std::cos(prob_parm.jet_angle)*vel_jet;
          s_ext[VELY] = scaled_vely_fluct;
          s_ext[VELZ] = scaled_velz_fluct + std::sin(prob_parm.jet_angle)*vel_jet;
      }
      else{
          s_ext[VELX] = 0.0;
          s_ext[VELY] = 0.0;
          s_ext[VELZ] = 0.0;
      }
#endif
  }
  
#if (AMREX_SPACEDIM == 2)
  // Set domain boundary Temp on the y-hi face when using NoSlipWallIsotherm
  if (idir == 1 && sgn == -1) {
    s_ext[TEMP] = prob_parm.Twall;
  }
#elif (AMREX_SPACEDIM == 3)
  // Set domain boundary Temp on the z-hi face when using NoSlipWallIsotherm
  if (idir == 2 && sgn == -1) {
    s_ext[TEMP] = prob_parm.Twall;
  }
  // Set domain boundary Temp on the x-hi face when using NoSlipWallIsotherm
  if (idir == 0 && sgn == -1) {
    s_ext[TEMP] = prob_parm.Twall;
  }
#endif
}

#define PELE_CUSTOM_BCNORMAL_AUX
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal_aux(
  const amrex::Real x[AMREX_SPACEDIM],
  amrex::CellData<amrex::Real>& aux_ext,
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm,
  pele::physics::PMF::PmfData::DataContainer const* pmf_data)
{
  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  amrex::Real X_air[NUM_SPECIES] = {0.0};
  amrex::Real Y_air[NUM_SPECIES] = {0.0};
  amrex::Real X_jet[NUM_SPECIES] = {0.0};
  amrex::Real Y_jet[NUM_SPECIES] = {0.0};
  amrex::Real eta = 0.0;
  constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;

  // Simple, air into air injection
  X_air[O2_ID] = 0.21;
  X_air[N2_ID] = 1.0 - X_air[O2_ID];

  eos.X2Y(X_air, Y_air);
  //eos.X2Y(X_jet, Y_jet);
  
  aux_ext[0] = 0.0;

#if (AMREX_SPACEDIM == 2)
  if (idir == 1 and sgn == 1) {
      amrex::Real r = sqrt(pow((0.0-x[0]),2));
#elif (AMREX_SPACEDIM == 3)
  if (idir == 0 and sgn == 1) {
      amrex::Real r = sqrt(pow((prob_parm.jet_center[2]-x[2]),2)/pow(std::cos(prob_parm.jet_angle),2)+pow((prob_parm.jet_center[0]-x[1]),2));
#endif

      amrex::Real N = 20.;
      amrex::Real radial_factor = 0.0;
      //streamwise velocity as a function of radial location
      if(r <= prob_parm.jet_rad){
          radial_factor = (1-pow(r/prob_parm.jet_rad,N));
          amrex::Real fu = 0.5*(tanh((time-prob_parm.t_si)/prob_parm.Ku)-tanh((time-prob_parm.t_ei)/prob_parm.Ku));
          aux_ext[0] = fu;
      }
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
zero_visc(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& beta,
  amrex::GeometryData const& /*geomdata*/,
  amrex::Box const& domainBox,
  const int dir,
  const int beta_comp,
  const int nComp)
{
  // amrex::ignore_unused(i,j,k,beta,geomdata,domainBox,dir,beta_comp,nComp);
  // We treat species when beta_comp == 0 and nComp == NUM_SPECIES
  // otherwise this routine could be called for other face diffusivity (Temp,
  // velocity, ...)
  if (beta_comp == 0 && nComp == NUM_SPECIES) {
    if (dir == 2 && k <= domainBox.smallEnd(dir)) {

      // Just over whole inlet
      for (int n = 0; n < nComp; n++) {
        beta(i, j, k, n) = 0.0;
      }
    }
  }
}


// struct EBhandle {
//     AMREX_GPU_DEVICE
//     AMREX_FORCE_INLINE
//     void bctype_eb(const amrex::Real xEBface[AMREX_SPACEDIM],
//                    int & EBflagType,
//                    amrex::GeometryData const& /*geomdata*/,
//                    ProbParm const& prob_parm) const
//     {
// 	//amrex::Real r = sqrt(pow((prob_parm.jet_center[0]-xEBface[0]),2)+pow((prob_parm.jet_center[1]-xEBface[1]),2)+pow((prob_parm.jet_center[2]-xEBface[2]),2));
// 	amrex::Real r = sqrt(pow((prob_parm.jet_center[0]-xEBface[0]),2)+pow((prob_parm.jet_center[1]-xEBface[1]),2));
//         if (r < prob_parm.jet_rad) {
//             EBflagType = pelelmex::BCTypeEB::inflow;
//         }
//     }
// 
//     AMREX_GPU_DEVICE
//     AMREX_FORCE_INLINE
//     void bcnormal_eb(const amrex::Real xEBface[AMREX_SPACEDIM],
//                      const amrex::Real EBfaceNormal[AMREX_SPACEDIM], // Outward normal
//                      const amrex::Real s_in[NVAR],
//                      amrex::Real s_ext[NVAR],
//                      const amrex::Real time,
//                      amrex::GeometryData const& /*geomdata*/,
//                      ProbParm const& prob_parm) const
//     {
// 
//         // Inflow
// 	amrex::Real r = sqrt(pow((prob_parm.jet_center[0]-xEBface[0]),2)+pow((prob_parm.jet_center[1]-xEBface[1]),2));
//         if (r < prob_parm.jet_rad) {
//             amrex::Real N = 20.;
//             amrex::Real radial_factor = 0.0;
//             
// 	    //streamwise velocity as a function of radial location
//             radial_factor = (1-pow(r/prob_parm.jet_rad,N));
//   
//             amrex::Real vel_jet = prob_parm.V_jet*radial_factor;
//       
// 	    //temporal velocity function
//             amrex::Real fu = 0.5*(tanh((time-prob_parm.t_si)/prob_parm.Ku)-tanh((time-prob_parm.t_ei)/prob_parm.Ku));
//             vel_jet *= fu;
// 
//             auto eos = pele::physics::PhysicsType::eos();
//             s_ext[VELX] = std::cos(prob_parm.jet_angle)*vel_jet;
//             s_ext[VELY] = 0.0;
//             s_ext[VELZ] = std::sin(prob_parm.jet_angle)*vel_jet;
// 
//             amrex::Real massfrac[NUM_SPECIES] = {0.0};
//             massfrac[O2_ID] = 0.233;
//             massfrac[N2_ID] = 1.0 - massfrac[O2_ID];
// 
//             amrex::Real P_cgs = prob_parm.P_mean * 10.0;
// 
//             s_ext[TEMP] = prob_parm.T_jet;
//             amrex::Real rho_cgs, RhoH_temp;
//             eos.PYT2R(P_cgs, massfrac, s_ext[TEMP], rho_cgs);
//             s_ext[DENSITY] = rho_cgs * 1.0e3;
// 
//             eos.TY2H(s_ext[TEMP], massfrac, RhoH_temp);
//             s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];
// 
//             for (int n = 0; n < NUM_SPECIES; n++) {
//                 s_ext[FIRSTSPEC+n] = massfrac[n] * s_ext[DENSITY];
//             }
//         }
//     }
// };


#endif
