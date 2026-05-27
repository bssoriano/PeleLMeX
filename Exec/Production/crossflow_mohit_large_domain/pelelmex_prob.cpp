#include <PeleLMeX.H>
#include <AMReX_ParmParse.H>

void
PeleLM::readProbParm() // NOLINT(readability-make-member-function-const)
{
  amrex::ParmParse ppeb("eb2");
  ppeb.get("inflow_x_lo", prob_parm->inflow_x_lo);
  ppeb.get("inflow_x_hi", prob_parm->inflow_x_hi);
  ppeb.get("inflow_z_lo", prob_parm->inflow_z_lo);
  ppeb.get("inflow_z_hi", prob_parm->inflow_z_hi);

  amrex::ParmParse pp("prob");

  pp.query("jet_radius", prob_parm->jet_radius);
  pp.query("T_mean", prob_parm->T_mean);
  pp.query("P_mean", prob_parm->P_mean);
  pp.query("meanFlowDir", prob_parm->meanFlowDir);
  pp.query("meanFlowMag", prob_parm->meanFlowMag);
  pp.query("jet_vel",prob_parm->jet_vel);
  pp.query("jet_vel",prob_parm->jet_vel);
  pp.query("jet_T",prob_parm->jet_T);
  pp.query("jet_Yfuel",prob_parm->jet_Yfuel);
  pp.query("t_tr",prob_parm->t_tr);
  pp.query("Ktr",prob_parm->Ktr);
  pp.query("Flame_IC",prob_parm->Flame_IC);
  pp.query("duration_injProd",prob_parm->duration_injProd);
  pp.query("start_time_injProd",prob_parm->start_time_injProd);
  
  //PeleLM::pmf_data.initialize();
  PeleLM::prob_parm->eosparm = PeleLM::eos_parms.device_parm();

  amrex::Vector<amrex::Real> local_inject_loc(
    AMREX_SPACEDIM, std::numeric_limits<amrex::Real>::lowest());
  pp.queryarr("jet_loc", local_inject_loc, 0, AMREX_SPACEDIM);

  for (int i = 0; i < AMREX_SPACEDIM; i++) {
    PeleLM::prob_parm->jet_loc[i] = local_inject_loc[i];
  }
//  PeleLM::prob_parm->Y_prod[O2_ID]  = 0.005135996986697483;
//  PeleLM::prob_parm->Y_prod[H2O_ID] = 0.12399164843694147;
//  PeleLM::prob_parm->Y_prod[CH4_ID] = 0.0;
//  PeleLM::prob_parm->Y_prod[CO_ID]  = 0.008763631883090546;
//  PeleLM::prob_parm->Y_prod[CO2_ID] = 0.1373434481760564;
//  PeleLM::prob_parm->Y_prod[N2_ID]  = 0.7247652745221234;


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

void
PeleLM::freeProbParm()
{
  PeleLM::pmf_data.deallocate();
}
