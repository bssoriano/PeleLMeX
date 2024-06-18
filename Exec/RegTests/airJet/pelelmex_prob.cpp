#include <PeleLMeX.H>
#include <pelelmex_prob.H>

void PeleLM::readProbParm()
{
   amrex::ParmParse pp("prob");

   pp.query("P_mean",   prob_parm->P_mean);
   pp.query("inj_start",   prob_parm->inj_start);
   pp.query("inj_dur",   prob_parm->inj_dur);
   pp.query("v_in",  prob_parm->v_in);
   pp.query("D",   prob_parm->D);
   pp.query("Z",   prob_parm->Z);
   pp.query("T_fu",  prob_parm->T_fu);
   pp.query("T_ox",  prob_parm->T_ox);
   pp.query("tau",  prob_parm->tau);

   amrex::Vector<std::string> sname;
   pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(sname);
   amrex::Real Y_pure_fuel[NUM_SPECIES] = {0.0};
   int o2_indx = -1;
   int n2_indx = -1;
   for (int n=0; n<sname.size(); n++)
   {
     if ( sname[n] == "O2") o2_indx = n;
     if ( sname[n] == "N2") n2_indx = n;
   }

   prob_parm->Y_ox[o2_indx] = 0.;
   prob_parm->Y_ox[n2_indx] = 1.0;

   for (int n = 0; n < NUM_SPECIES; n++)
   {
       prob_parm->Y_fuel[n] = prob_parm->Y_ox[n];
   }

   auto eos = pele::physics::PhysicsType::eos();
   eos.TY2H(prob_parm->T_fu, prob_parm->Y_fuel, prob_parm->H_fuel);
   eos.TY2H(prob_parm->T_ox, prob_parm->Y_ox,   prob_parm->H_ox);

   auto problo = geom[0].ProbLo();
   auto probhi = geom[0].ProbHi();
   prob_parm->center_xy[0] = 0.5 * (probhi[0] + problo[0]);
   prob_parm->center_xy[1] = 0.5 * (probhi[1] + problo[1]);
}
