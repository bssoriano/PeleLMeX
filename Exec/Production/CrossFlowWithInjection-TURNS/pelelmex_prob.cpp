#include <PeleLMeX.H>
#include <AMReX_ParmParse.H>

void
PeleLM::readProbParm() // NOLINT(readability-make-member-function-const)
{
  amrex::ParmParse pp("prob");

  pp.query("injection_diam", prob_parm->injection_diam);
  pp.query("T_mean", prob_parm->T_mean);
  pp.query("P_mean", prob_parm->P_mean);
  pp.query("meanFlowDir", prob_parm->meanFlowDir);
  pp.query("meanFlowMag", prob_parm->meanFlowMag);

  amrex::Vector<amrex::Real> local_inject_loc(
    AMREX_SPACEDIM, std::numeric_limits<amrex::Real>::lowest());
  pp.queryarr("inject_loc", local_inject_loc, 0, AMREX_SPACEDIM);

  for (int i = 0; i < AMREX_SPACEDIM; i++) {
    PeleLM::prob_parm->nject_loc[i] = local_inject_loc[i];
  }

  // if (!m_incompressible) {
  //    auto& trans_parm = PeleLM::trans_parms.host_parm();
  //    amrex::ParmParse pptr("transport");
  //    pp.query("const_viscosity", trans_parm.const_viscosity);
  //    pp.query("const_bulk_viscosity", trans_parm.const_bulk_viscosity);
  //    pp.query("const_conductivity", trans_parm.const_conductivity);
  //    pp.query("const_diffusivity", trans_parm.const_diffusivity);
  //    PeleLM::trans_parms.sync_to_device();
  // }
}
