#include <PeleLMeX.H>
#include <AMReX_ParmParse.H>
#include <pelelmex_prob.H>

void PeleLM::readProbParm()
{
  amrex::ParmParse pp("prob");

   pp.query("P_mean"       , PeleLM::prob_parm->P_mean);
   pp.query("Zst"          ,    PeleLM::prob_parm->Zst);
   pp.query("T_in"         ,   PeleLM::prob_parm->T_in);
   pp.query("T_init"       ,   PeleLM::prob_parm->T0);
   pp.query("U_b"          ,    PeleLM::prob_parm->U_b);
   pp.query("standoff"     , PeleLM::prob_parm->standoff);


   PeleLM::prob_parm->bathID = N2_ID;  
   // PeleLM::prob_parm->fuelID = NXC7H16_ID;  
   PeleLM::prob_parm->fuelID = POSF10325_ID;  
   PeleLM::prob_parm->oxidID = O2_ID; 

   PeleLM::pmf_data.initialize();

  amrex::Real moments[NUM_SOOT_MOMENTS + 1] = {0.0};
  if (PeleLM::do_soot_solve) {
    SootData* const sd = PeleLM::soot_model->getSootData();
    sd->initialSmallMomVals(moments);
  }
  for (int n = 0; n < NUM_SOOT_MOMENTS + 1; ++n) {
    PeleLM::prob_parm->soot_vals[n] = moments[n];
  }
}
