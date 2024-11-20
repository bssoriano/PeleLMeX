#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{

  amrex::Real P_mean = 101325.0;
  amrex::Real T_in     = 300.0;
  amrex::Real T0       = 300.0;
  amrex::Real T_bottom = 300.0;
  amrex::Real T_wall   = 300.0;
  amrex::Real vel = 0.0;
  amrex::Real X_O2 = 0.0;
  amrex::Real X_N2 = 0.0;
  amrex::Real X_H2O = 0.0;
  amrex::Real X_H2O_inj = 0.0;
  amrex::Real Vin = 0.0;
  amrex::Real inj_diam = 0.0;


};
#endif
