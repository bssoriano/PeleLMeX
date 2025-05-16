#ifndef PELELM_PROB_PARM_H
#define PELELM_PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm : amrex::Gpu::Managed
{
  // Shared params
  amrex::Real T_mean = 298.0_rt;
  amrex::Real P_mean = 101325.0_rt;
  amrex::Real meanFlowMag = 0.0;
  amrex::Real lowFlowMag = 3.0;
  amrex::Real injection_diam = 0.0;
  int meanFlowDir = 1;
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> inject_loc = {0.0, 0.0, 0.0};

};
#endif
