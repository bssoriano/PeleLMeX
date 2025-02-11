#ifndef PELELM_PROB_PARM_H
#define PELELM_PROB_PARM_H

#include <AMReX_REAL.H>
#include <PelePhysics.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
  amrex::Real P_mean = 101325.0_rt;
  amrex::Real T_mean = 300.0_rt;
  amrex::Real Twall  = 300.0_rt;
  amrex::Real T_jet = 1000.0_rt;
  amrex::Real V_jet = 1.0_rt;
  amrex::Real jet_rad = 0.25_rt;
  amrex::Real bl_thickness = 0.01_rt;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y_chamber = {{0.0}};
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y_jet     = {{0.0}};
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y_prechamber = {{0.0}};
  amrex::Real t_si = 0.0_rt;
  amrex::Real t_ei = 0.0_rt;
  amrex::Real Ku   = 0.0_rt;
  amrex::Real t_tr = 0.0_rt;
  amrex::Real Ktr  = 0.0_rt;
  amrex::Real jet_angle = 0.0_rt; //in rad. pi/2 correspond to wall normal
  amrex::GpuArray<amrex::Real, 3> jet_center = {{0.0}};


  int input_resolution = 0;
  amrex::Real urms0 = 1.0;
  amrex::Real uin_norm = 1.0;

  amrex::Real Linput = 0.0;

  amrex::Real* d_uinput = nullptr;
  amrex::Real* d_vinput = nullptr;
  amrex::Real* d_winput = nullptr;
  amrex::Real* d_xarray = nullptr;
  amrex::Real* d_xdiff = nullptr;
};
#endif
