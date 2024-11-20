#include <PeleLMeX.H>
#include <AMReX_ParmParse.H>

void PeleLM::readProbParm()
{
  amrex::ParmParse pp("prob");

   pp.query("P_mean"       , PeleLM::prob_parm->P_mean);
   pp.query("T_in"         ,   PeleLM::prob_parm->T_in);
   pp.query("T_init"       ,   PeleLM::prob_parm->T0);
   pp.query("T_bottom"     ,   PeleLM::prob_parm->T_bottom);
   pp.query("T_wall"       ,   PeleLM::prob_parm->T_wall);
   pp.query("velocity"     ,   PeleLM::prob_parm->vel);
   pp.query("X_O2"         ,   PeleLM::prob_parm->X_O2);
   pp.query("X_N2"         ,   PeleLM::prob_parm->X_N2);
   pp.query("X_H2O"        ,   PeleLM::prob_parm->X_H2O);
   pp.query("X_H2O_inj"    ,   PeleLM::prob_parm->X_H2O_inj);
   pp.query("Vin"          ,   PeleLM::prob_parm->Vin);
   pp.query("inj_diam"     ,   PeleLM::prob_parm->inj_diam);


}
