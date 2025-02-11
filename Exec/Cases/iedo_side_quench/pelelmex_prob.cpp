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

  prob_parm->jet_angle = Pi*jet_angle_deg/180; //conversion to rad

  // ------ Initializing jet composition -------
  PeleLM::prob_parm->Y_jet[H2_ID] = 1.310e-03;
  PeleLM::prob_parm->Y_jet[O2_ID] = 7.649e-03;
  PeleLM::prob_parm->Y_jet[H2O_ID] = 2.375e-01;
  PeleLM::prob_parm->Y_jet[H_ID] = 7.359e-05;
  PeleLM::prob_parm->Y_jet[O_ID] = 4.639e-04;
  PeleLM::prob_parm->Y_jet[OH_ID] = 7.724e-03;
  PeleLM::prob_parm->Y_jet[HO2_ID] = 1.057e-05;
  PeleLM::prob_parm->Y_jet[H2O2_ID] = 2.291e-06;
  PeleLM::prob_parm->Y_jet[N2_ID] = 7.452e-01;

  // ------ Initializing H2 premixed composition for the main chamber -------
  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  amrex::Real a = 0.5;
  molefrac[O2_ID] = 1.0 / ( 1.0 + phi_main / a + 0.79 / 0.21 );
  molefrac[H2_ID] = phi_main * molefrac[O2_ID] / a;
  molefrac[N2_ID] = 1.0 - molefrac[O2_ID] - molefrac[H2_ID];

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

  eos.X2Y(molefrac,massfrac);

  for (int n = 0; n < NUM_SPECIES; n++){
    (PeleLM::prob_parm->Y_prechamber)[n] = massfrac[n];
  }

  // ----- Read csv file with velocity fluctuations -----
  
//  ProbParm local_prob_parm;
//  std::string datafile;
//  pp.query("input_name", datafile);
//  pp.query("input_resolution", local_prob_parm.input_resolution);
//  int binfmt = 0; // Default is ASCII format
//  pp.query("urms0", local_prob_parm.urms0);
//
//  // Read initial velocity field
//  const size_t nx = local_prob_parm.input_resolution;
//  const size_t ny = local_prob_parm.input_resolution;
//  const size_t nz = local_prob_parm.input_resolution;
//  amrex::Vector<amrex::Real> data(
//    nx * ny * nz * 6); /* this needs to be double */
// 
//  read_csv(datafile, nx, ny, nz, data);
//
//  // Extract position and velocities
//  amrex::Vector<amrex::Real> xinput(nx * ny * nz);
//  amrex::Vector<amrex::Real> uinput(nx * ny * nz);
//  amrex::Vector<amrex::Real> vinput(nx * ny * nz);
//  amrex::Vector<amrex::Real> winput(nx * ny * nz);
//  amrex::Vector<amrex::Real> xdiff(nx);
//  amrex::Vector<amrex::Real> xarray(nx);
//
//  for (long i = 0; i < xinput.size(); i++) {
//    xinput[i] = data[0 + i * 6];
//    uinput[i] =
//      data[3 + i * 6] * local_prob_parm.urms0 / local_prob_parm.uin_norm;
//    vinput[i] =
//      data[4 + i * 6] * local_prob_parm.urms0 / local_prob_parm.uin_norm;
//    winput[i] =
//      data[5 + i * 6] * local_prob_parm.urms0 / local_prob_parm.uin_norm;
//  }
//
//  // Get the xarray table and the differences.
//  for (long i = 0; i < xarray.size(); i++) {
//    xarray[i] = xinput[i];
//  }
//  std::adjacent_difference(xarray.begin(), xarray.end(), xdiff.begin());
//  xdiff[0] = xdiff[1];
//
//  // Make sure the search array is increasing
//  if (not std::is_sorted(xarray.begin(), xarray.end())) {
//    amrex::Abort("Error: non ascending x-coordinate array.");
//  }
//
//  // Pass data to the local_prob_parm
//  local_prob_parm.Linput = xarray[nx - 1] + 0.5 * xdiff[nx - 1];
//
//  local_prob_parm.d_xarray =
//    (amrex::Real*)amrex::The_Arena()->alloc(nx * sizeof(amrex::Real));
//  local_prob_parm.d_xdiff =
//    (amrex::Real*)amrex::The_Arena()->alloc(nx * sizeof(amrex::Real));
//  local_prob_parm.d_uinput =
//    (amrex::Real*)amrex::The_Arena()->alloc(nx * ny * nz * sizeof(amrex::Real));
//  local_prob_parm.d_vinput =
//    (amrex::Real*)amrex::The_Arena()->alloc(nx * ny * nz * sizeof(amrex::Real));
//  local_prob_parm.d_winput =
//    (amrex::Real*)amrex::The_Arena()->alloc(nx * ny * nz * sizeof(amrex::Real));
//
//  for (unsigned long i = 0; i < nx; i++) {
//    local_prob_parm.d_xarray[i] = xarray[i];
//    local_prob_parm.d_xdiff[i] = xdiff[i];
//  }
//  for (unsigned long i = 0; i < nx * ny * nz; i++) {
//    local_prob_parm.d_uinput[i] = uinput[i];
//    local_prob_parm.d_vinput[i] = vinput[i];
//    local_prob_parm.d_winput[i] = winput[i];
//  }
//
//  // Initialize PeleLM::prob_parm container
//  PeleLM::prob_parm->d_xarray =
//    (amrex::Real*)amrex::The_Arena()->alloc(nx * sizeof(amrex::Real));
//  PeleLM::prob_parm->d_xdiff =
//    (amrex::Real*)amrex::The_Arena()->alloc(nx * sizeof(amrex::Real));
//  PeleLM::prob_parm->d_uinput =
//    (amrex::Real*)amrex::The_Arena()->alloc(nx * ny * nz * sizeof(amrex::Real));
//  PeleLM::prob_parm->d_vinput =
//    (amrex::Real*)amrex::The_Arena()->alloc(nx * ny * nz * sizeof(amrex::Real));
//  PeleLM::prob_parm->d_winput =
//    (amrex::Real*)amrex::The_Arena()->alloc(nx * ny * nz * sizeof(amrex::Real));
//
//  // Copy into PeleLM::prob_parm: CPU only for now
//  PeleLM::prob_parm->d_xarray =
//    (amrex::Real*)amrex::The_Arena()->alloc(nx * sizeof(amrex::Real));
//  std::memcpy(
//    &PeleLM::prob_parm->d_xarray, &local_prob_parm.d_xarray,
//    sizeof(local_prob_parm.d_xarray));
//  std::memcpy(
//    &PeleLM::prob_parm->d_xdiff, &local_prob_parm.d_xdiff,
//    sizeof(local_prob_parm.d_xdiff));
//  std::memcpy(
//    &PeleLM::prob_parm->d_uinput, &local_prob_parm.d_uinput,
//    sizeof(local_prob_parm.d_uinput));
//  std::memcpy(
//    &PeleLM::prob_parm->d_vinput, &local_prob_parm.d_vinput,
//    sizeof(local_prob_parm.d_vinput));
//  std::memcpy(
//    &PeleLM::prob_parm->d_winput, &local_prob_parm.d_winput,
//    sizeof(local_prob_parm.d_winput));
//  PeleLM::prob_parm->Linput = local_prob_parm.Linput;
//  PeleLM::prob_parm->input_resolution = local_prob_parm.input_resolution;
//  PeleLM::prob_parm->urms0 = local_prob_parm.urms0;
//  PeleLM::prob_parm->uin_norm = local_prob_parm.uin_norm;
}
