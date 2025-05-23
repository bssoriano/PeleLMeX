#ifndef PELELM_PROB_PARM_H
#define PELELM_PROB_PARM_H

#include <AMReX_REAL.H>
#include <PelePhysics.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm : amrex::Gpu::Managed
{
  // Shared params
  amrex::Real T_mean = 298.0_rt;
  amrex::Real P_mean = 101325.0_rt;
  amrex::Real meanFlowMag = 0.0;
  amrex::Real lowFlowMag = 3.0;
  amrex::Real jet_radius = 0.0;
  int meanFlowDir = 1;
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> jet_loc = {{0.0}};
  amrex::Real jet_vel = 0.0;
  amrex::Real jet_T = 0.0;
  amrex::Real jet_Yfuel = 0.0;
  amrex::Real Ktr  = 0.0_rt;
  amrex::Real t_tr = 0.0_rt;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y_prod = {{0.0}};

};
#endif
