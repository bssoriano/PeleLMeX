#ifndef PELELM_PROB_H
#define PELELM_PROB_H

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_SPACE.H>

#include <PeleLMeX_Index.H>
#include <pelelmex_prob_parm.H>
#include <PMFData.H>
#include <PelePhysics.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pelelmex_initdata(
  int i,
  int j,
  int k,
  int is_incompressible,
  amrex::Array4<amrex::Real> const& state,
  amrex::Array4<amrex::Real> const& /*aux*/,
  amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm,
  pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
{

  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();

  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];

  auto eos = pele::physics::PhysicsType::eos();

  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[O2_ID] = 0.233;
  massfrac[N2_ID] = 0.767;

  switch (prob_parm.meanFlowDir) {
  case 1:
    AMREX_D_TERM(state(i, j, k, VELX) = prob_parm.meanFlowMag;
                 , state(i, j, k, VELY) = 0.0;, state(i, j, k, VELZ) = 0.0);
    break;
  case -1:
    AMREX_D_TERM(state(i, j, k, VELX) = -prob_parm.meanFlowMag;
                 , state(i, j, k, VELY) = 0.0;, state(i, j, k, VELZ) = 0.0);
    break;
  case 2:
    AMREX_D_TERM(state(i, j, k, VELX) = 0.0;
                 , state(i, j, k, VELY) = prob_parm.meanFlowMag;
                 , state(i, j, k, VELZ) = 0.0);
    break;
  case -2:
    AMREX_D_TERM(state(i, j, k, VELX) = 0.0;
                 , state(i, j, k, VELY) = -prob_parm.meanFlowMag;
                 , state(i, j, k, VELZ) = 0.0);
    break;
  case 3:
    AMREX_D_TERM(state(i, j, k, VELX) = prob_parm.meanFlowMag;
                 , state(i, j, k, VELY) = prob_parm.meanFlowMag;
                 , state(i, j, k, VELZ) = 0.0);
    break;
  case -3:
    AMREX_D_TERM(state(i, j, k, VELX) = -prob_parm.meanFlowMag;
                 , state(i, j, k, VELY) = -prob_parm.meanFlowMag;
                 , state(i, j, k, VELZ) = 0.0);
    break;
  default:
    break;
  }

  if (!(is_incompressible == 0)) {
    return;
  }

  state(i, j, k, TEMP) = prob_parm.T_mean;
  
//------------- ADD A FLAME TO THE INITIAL CONDITION (TSUJI BURNER ASYMPTOTIC SOLUTION) ----------------
 switch (prob_parm.Flame_IC) {
    case 0:
        // Disabled
        //amrex::Print() << "Flame initialization disabled (0)\n";
        state(i, j, k, TEMP) = prob_parm.T_mean;
        massfrac[CH4_ID] = 0.0;
        massfrac[O2_ID] = (1.0 - massfrac[CH4_ID]) * 0.233;
        massfrac[N2_ID] =  1.0 - massfrac[CH4_ID] - massfrac[O2_ID];
        break;
    case 1:
        // Standard flame
        //amrex::Print() << "Standard flame initialization (1)\n";
        
        double y_sf = 0.4;
        double V_asymp = 10.0;
        double Pe_asymp = 55.0;
  
        double y_or = 0.272;
        double x_or = 0.136;

        //Temperature (Z>1): T = (q/(S+1.0)) * ( H - ( Z - 1.0)/S)
        //Temperature (Z<1): T = (q/(S+1.0)) * ( H - (-Z + 1.0)  )
  
        //S =  ( nu * YF_b) / YO_oo
          //nu = 3.99d0,    !for Methane , 3.51d0 for n-Heptane
        //YF_b = 1.0,
        //YO_oo = 0.21d0,
        double YF_b = 0.06;
        double S_IC = 3.99 * YF_b / 0.21;

        //Q = (q_dim * YF_b) / (cp_tot * Too) !Heat of combustion
        double Q_IC = (50.15e6 * YF_b) / (1937.35 * 298.0);
        //Entalphy:          H = (S+1.d0)*Tinf/q + 1.d0
        double H_IC = (S_IC+1.0) * 1.0 / Q_IC + 1.0;

        if (y < y_or) {
        massfrac[CH4_ID] = 0.0;
        state(i, j, k, TEMP) = prob_parm.T_mean;
          } else {

        // Calculate the unnormalized value
        double Z_unnormalized = (sqrt(y_sf/(y-y_or)) * exp(-(V_asymp*Pe_asymp/4.0) *((x-x_or)*(x-x_or))/(y-y_or)));
          
        // Find the maximum possible value (occurs at x=0.2)
        double Z_max_value = sqrt(y_sf/(y-y_or)); // When x=0.2, the exponential term becomes 1
          
        // Normalize between 0 and 1
        massfrac[CH4_ID] = YF_b * Z_unnormalized / Z_max_value;

        if (Z_unnormalized > 1) {
        state(i, j, k, TEMP) = prob_parm.T_mean * (Q_IC/(S_IC+1.0))*(H_IC - ( Z_unnormalized - 1.0)/S_IC);
        } else {
        state(i, j, k, TEMP) = prob_parm.T_mean * (Q_IC/(S_IC+1.0))*(H_IC - (-Z_unnormalized + 1.0));
        }

        }
        massfrac[O2_ID] = (1.0 - massfrac[CH4_ID]) * 0.233;
        massfrac[N2_ID] =  1.0 - massfrac[CH4_ID] - massfrac[O2_ID];


        if (state(i, j, k, TEMP) < 298.0) {
        state(i, j, k, TEMP) = prob_parm.T_mean;
        }

        break;
 }
 
//------------------------------- END IC FLAME ----------------------------------------------------------------      


//----------------- BACKWARD FACE STEP IC FLAME TYPE -------------------
  // massfrac[CH4_ID] = 0.0;
  // massfrac[O2_ID] = 0.233;
  // massfrac[N2_ID] = 1 - massfrac[CH4_ID] - massfrac[O2_ID];
  // 
  // if (x < 0.136 && y > 0.273) {
  //   state(i, j, k, TEMP) = 1600;
  //   massfrac[CH4_ID] = 0.0445;
  //   massfrac[O2_ID] = 0.233;
  //   massfrac[N2_ID] = 1 - massfrac[CH4_ID] - massfrac[O2_ID];
  //   state(i, j, k, VELY) = 0.0;
  //   state(i, j, k, VELX) = 0.0;
  // } else {
  //   state(i, j, k, TEMP) = prob_parm.T_mean;
  // }
//-------------------------------------------------------------------

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
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  //const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  const amrex::Real x[AMREX_SPACEDIM],
  const int /*m_nAux*/,
  amrex::Real s_ext[NVAR],
  const int /*idir*/,
  const int /*sgn*/,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParm const& prob_parm,
  pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
{

  amrex::Real r = sqrt(pow((prob_parm.jet_loc[1]-x[1]),2)+pow((prob_parm.jet_loc[2]-x[2]),2));
  
  amrex::Real massfrac[NUM_SPECIES] = {0.0};

  auto eos = pele::physics::PhysicsType::eos();

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

  massfrac[O2_ID] = 0.333;
  massfrac[N2_ID] = 0.667;

  s_ext[TEMP] = prob_parm.T_mean;

  amrex::Real rho_cgs, P_cgs, RhoH_temp;
  P_cgs = prob_parm.P_mean * 10.0;

  eos.PYT2R(P_cgs, massfrac, s_ext[TEMP], rho_cgs);
  s_ext[DENSITY] = rho_cgs * 1.0e3;

  eos.TY2H(s_ext[TEMP], massfrac, RhoH_temp);
  s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY]; // CGS -> MKS conversion

  for (int n = 0; n < NUM_SPECIES; n++) {
    s_ext[FIRSTSPEC + n] = massfrac[n] * s_ext[DENSITY];
  }
}

struct EBhandle {
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void bctype_eb(const amrex::Real xEBface[AMREX_SPACEDIM],
                   int & EBflagType,
                   amrex::GeometryData const& /*geomdata*/,
                   ProbParm const& prob_parm) const
    {
#if (AMREX_SPACEDIM == 2)
        amrex::Real r = sqrt(pow((prob_parm.jet_loc[1]-xEBface[1]),2));
#elif (AMREX_SPACEDIM == 3)
        amrex::Real r = sqrt(pow((prob_parm.jet_loc[1]-xEBface[1]),2)+pow((prob_parm.jet_loc[2]-xEBface[2]),2));
#endif
	if (r <= prob_parm.jet_radius) {
            EBflagType = pelelmex::BCTypeEB::inflow;
        }
        
    }
  
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void bcnormal_eb(const amrex::Real xEBface[AMREX_SPACEDIM],
                     const amrex::Real EBfaceNormal[AMREX_SPACEDIM], // Outward normal
                     const amrex::Real s_in[NVAR],
                     amrex::Real s_ext[NVAR],
                     const amrex::Real /*time*/,
                     amrex::GeometryData const& /*geomdata*/,
                     ProbParm const& prob_parm) const
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

            auto eos = pele::physics::PhysicsType::eos();
            s_ext[VELX] = prob_parm.jet_vel;
            s_ext[VELY] = 0.0;
            s_ext[VELZ] = 0.0;
            
	    amrex::Real massfrac[NUM_SPECIES] = {0.0};
            massfrac[CH4_ID] = prob_parm.jet_Yfuel;
            massfrac[O2_ID] = (1.0 - massfrac[CH4_ID]) * 0.233;
            massfrac[N2_ID] =  1.0 - massfrac[CH4_ID] - massfrac[O2_ID];

            amrex::Real P_cgs = prob_parm.P_mean * 10.0;

            s_ext[TEMP] = prob_parm.jet_T;
            amrex::Real rho_cgs, RhoH_temp;
            eos.PYT2R(P_cgs, massfrac, s_ext[TEMP], rho_cgs);
            s_ext[DENSITY] = rho_cgs * 1.0e3;

            eos.TY2H(s_ext[TEMP], massfrac, RhoH_temp);
            s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];

            for (int n = 0; n < NUM_SPECIES; n++) {
                s_ext[FIRSTSPEC+n] = massfrac[n] * s_ext[DENSITY];
            }
        }
    }
};

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
