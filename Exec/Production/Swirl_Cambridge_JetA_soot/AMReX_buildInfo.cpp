
namespace amrex {

const char* buildInfoGetBuildDate() {

  static const char BUILD_DATE[] = "2024-03-25 15:42:13.731034";
  return BUILD_DATE;
}

const char* buildInfoGetBuildDir() {

  static const char BUILD_DIR[] = "/home/bsouzas/my_LMeX/PeleLMeX/Exec/Production/Swirl_Cambridge_JetA_soot";
  return BUILD_DIR;
}

const char* buildInfoGetBuildMachine() {

  static const char BUILD_MACHINE[] = "Linux blaze 4.18.0-193.6.3.el8_2.x86_64 #1 SMP Wed Jun 10 11:09:32 UTC 2020 x86_64 x86_64 x86_64 GNU/Linux";
  return BUILD_MACHINE;
}

const char* buildInfoGetAMReXDir() {

  static const char AMREX_DIR[] = "/home/bsouzas/my_LMeX/PeleLMeX/Submodules/PelePhysics/Submodules/amrex";
  return AMREX_DIR;
}

const char* buildInfoGetComp() {

  static const char COMP[] = "gnu";
  return COMP;
}

const char* buildInfoGetCompVersion() {

  static const char COMP_VERSION[] = "9.2.0";
  return COMP_VERSION;
}

// deprecated
const char* buildInfoGetFcomp() {

  static const char FCOMP[] = "";
  return FCOMP;
}

// deprecated
const char* buildInfoGetFcompVersion() {

  static const char FCOMP_VERSION[] = "";
  return FCOMP_VERSION;
}

const char* buildInfoGetCXXName() {

  static const char CXX_comp_name[] = "";
  return CXX_comp_name;
}

const char* buildInfoGetFName() {

  static const char F_comp_name[] = "";
  return F_comp_name;
}

const char* buildInfoGetCXXFlags() {

  static const char CXX_flags[] = "";
  return CXX_flags;
}

const char* buildInfoGetFFlags() {

  static const char F_flags[] = "";
  return F_flags;
}

const char* buildInfoGetLinkFlags() {

  static const char link_flags[] = "";
  return link_flags;
}

const char* buildInfoGetLibraries() {

  static const char libraries[] = "";
  return libraries;
}

const char* buildInfoGetAux(int i) {

  //static const char AUX1[] = "${AUX[1]}";

  static const char EMPT[] = "";

  switch(i)
  {

    default: return EMPT;
  }
}

int buildInfoGetNumModules() {
  // int const num_modules = X;
  int const num_modules = 0;

  return num_modules;
}

const char* buildInfoGetModuleName(int i) {

  //static const char MNAME1[] = "${MNAME[1]}";

  static const char EMPT[] = "";

  switch(i)
  {

    default: return EMPT;
  }
}

const char* buildInfoGetModuleVal(int i) {

  //static const char MVAL1[] = "${MVAL[1]}";

  static const char EMPT[] = "";

  switch(i)
  {

    default: return EMPT;
  }
}

const char* buildInfoGetGitHash(int i) {

  //static const char HASH1[] = "${GIT[1]}";
  static const char HASH1[] = "92d7705";
  static const char HASH2[] = "24.02-31-g97512176d";
  static const char HASH3[] = "v23.03-98-g50b990b8";
  static const char HASH4[] = "24.01-5-gd61ab60";
  static const char HASH5[] = "v6.6.2";

  static const char EMPT[] = "";

  switch(i)
  {
    case 1: return HASH1;
    case 2: return HASH2;
    case 3: return HASH3;
    case 4: return HASH4;
    case 5: return HASH5;

    default: return EMPT;
  }
}

const char* buildInfoGetBuildGitHash() {

  //static const char HASH[] = "${GIT}";
  static const char HASH[] = "";


  return HASH;
}

const char* buildInfoGetBuildGitName() {

  //static const char NAME[] = "";
  static const char NAME[] = "";


  return NAME;
}

#ifdef AMREX_USE_CUDA
const char* buildInfoGetCUDAVersion() {

  static const char CUDA_VERSION[] = "";
  return CUDA_VERSION;
}
#endif

}
