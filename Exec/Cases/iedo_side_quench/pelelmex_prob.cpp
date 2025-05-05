#include <PeleLMeX.H>
#include <AMReX_ParmParse.H>

AMREX_FORCE_INLINE
std::string
read_file(std::ifstream& in)
{
  return static_cast<std::stringstream const&>(
           std::stringstream() << in.rdbuf())
    .str();
}

// -----------------------------------------------------------
// Read a csv file
// INPUTS/OUTPUTS:
// iname => filename
// nx    => input resolution
// ny    => input resolution
// nz    => input resolution
// data  <= output data
// -----------------------------------------------------------
void
read_csv(
  const std::string& iname,
  const size_t nx,
  const size_t ny,
  const size_t nz,
  amrex::Vector<amrex::Real>& data)
{
  std::ifstream infile(iname, std::ios::in);
  const std::string memfile = read_file(infile);
  if (not infile.is_open()) {
    amrex::Abort("Unable to open input file " + iname);
  }
  infile.close();
  std::istringstream iss(memfile);

  // Read the file
  size_t nlines = 0;
  std::string firstline;
  std::string line;
  std::getline(iss, firstline); // skip header
  while (getline(iss, line)) {
    ++nlines;
  }

  // Quick sanity check
  if (nlines != nx * ny * nz) {
    amrex::Abort(
      "Number of lines in the input file (= " + std::to_string(nlines) +
      ") does not match the input resolution (=" + std::to_string(nx) + ")");
  }

  // Read the data from the file
  iss.clear();
  iss.seekg(0, std::ios::beg);
  std::getline(iss, firstline); // skip header
  int cnt = 0;
  while (std::getline(iss, line)) {
    std::istringstream linestream(line);
    std::string value;
    while (getline(linestream, value, ',')) {
      std::istringstream sinput(value);
      sinput >> data[cnt];
      cnt++;
    }
  }
}

void
PeleLM::readProbParm() // NOLINT(readability-make-member-function-const)
{
  amrex::ParmParse pp("prob");
  
  amrex::Real phi_main = 0.0;
  amrex::Real phi_prechamber = 0.0;
  constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;

  pp.query("P_mean", prob_parm->P_mean);
  pp.query("T_mean", prob_parm->T_mean);
  pp.query("T_wall", prob_parm->Twall);
  pp.query("T_jet", prob_parm->T_jet);
  pp.query("V_jet", prob_parm->V_jet);
  pp.query("phi_chamber", phi_main);
  pp.query("phi_prechamber", phi_prechamber);
  pp.query("jet_rad", prob_parm->jet_rad);
  pp.query("bl_thickness", prob_parm->bl_thickness);
  
  pp.query("t_si",prob_parm->t_si);
  pp.query("t_ei",prob_parm->t_ei);
  pp.query("Ku",prob_parm->Ku);
  pp.query("t_tr",prob_parm->t_tr);
  pp.query("Ktr",prob_parm->Ktr);
  amrex::Real jet_angle_deg = 0.0;
  pp.query("jet_angle",jet_angle_deg);
  
  PeleLM::pmf_data.initialize();

  prob_parm->jet_angle = Pi*jet_angle_deg/180; //conversion to radians
  prob_parm->jet_center[0] = 0.0;
  prob_parm->jet_center[1] = 0.0;
  amrex::Real z_over_D = 1.25; // jet is located in z at 1.25 jet diameters 
  prob_parm->jet_center[2] = prob_parm->jet_rad*2.*z_over_D;

  // ------ Initializing jet composition -------
  PeleLM::prob_parm->Y_jet[H_ID] = 6.428e-05;
  PeleLM::prob_parm->Y_jet[H2_ID] = 1.207e-03;
  PeleLM::prob_parm->Y_jet[O_ID] = 3.924e-04;
  PeleLM::prob_parm->Y_jet[OH_ID] = 6.746e-03;
  PeleLM::prob_parm->Y_jet[OHV_ID] = 2.147e-10;
  PeleLM::prob_parm->Y_jet[H2O_ID] = 2.393e-01;
  PeleLM::prob_parm->Y_jet[O2_ID] = 7.019e-03;
  PeleLM::prob_parm->Y_jet[N2_ID] = 7.452e-01;
  PeleLM::prob_parm->Y_jet[HO2_ID] = 7.768e-06;
  PeleLM::prob_parm->Y_jet[AR_ID] = 0.000e+00;
  PeleLM::prob_parm->Y_jet[H2O2_ID] = 1.580e-06;
  PeleLM::prob_parm->Y_jet[HE_ID] = 0.000e+00;

//  PeleLM::prob_parm->Y_jet[O2_ID] = 0.233;
//  PeleLM::prob_parm->Y_jet[N2_ID] = 0.767;

  // ------ Initializing H2 premixed composition for the main chamber -------
  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  amrex::Real a = 0.5;
  molefrac[O2_ID] = 1.0 / ( 1.0 + phi_main / a + 0.79 / 0.21 );
  molefrac[H2_ID] = phi_main * molefrac[O2_ID] / a;
  molefrac[N2_ID] = 1.0 - molefrac[O2_ID] - molefrac[H2_ID];
  
//  molefrac[N2_ID] = 0.79;
//  molefrac[O2_ID] = 0.21;

  auto eos = pele::physics::PhysicsType::eos();

  eos.X2Y(molefrac,massfrac);

  for (int n = 0; n < NUM_SPECIES; n++){
    (PeleLM::prob_parm->Y_chamber)[n] = massfrac[n];
  }
  
  // ------ Initializing H2 premixed composition for the pre-chamber -------
  for (int n = 0; n < NUM_SPECIES; n++){
     molefrac[n] = 0.0;
     massfrac[n] = 0.0;
  }

  molefrac[O2_ID] = 1.0 / ( 1.0 + phi_prechamber / a + 0.79 / 0.21 );
  molefrac[H2_ID] = phi_prechamber * molefrac[O2_ID] / a;
  molefrac[N2_ID] = 1.0 - molefrac[O2_ID] - molefrac[H2_ID];

  //molefrac[N2_ID] = 0.79;
  //molefrac[O2_ID] = 0.21;

  eos.X2Y(molefrac,massfrac);

  for (int n = 0; n < NUM_SPECIES; n++){
    (PeleLM::prob_parm->Y_prechamber)[n] = massfrac[n];
  }

}
