#include <PeleLMeX.H>
#include <AMReX_ParmParse.H>

void
PeleLM::readProbParm()
{
  amrex::ParmParse pp("prob");

  std::string type;
  pp.query("P_mean", PeleLM::prob_parm->P_mean);
  pp.query("init_T", PeleLM::prob_parm->T0);
  pp.query("init_vel", PeleLM::prob_parm->vel);
  pp.query("init_Y_N2", PeleLM::prob_parm->Y_N2);
  pp.query("init_Y_O2", PeleLM::prob_parm->Y_O2);
  pp.query("init_Y_CO2", PeleLM::prob_parm->Y_CO2);
  pp.query("init_Y_H2O", PeleLM::prob_parm->Y_H2O);
}
