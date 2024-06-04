#ifndef MECHANISM_H
#define MECHANISM_H

#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

/* Elements
0  O
1  H
2  C
3  N
4  Ar
5  He
*/

// Species
#define POSF10325_ID 0
#define O2_ID 1
#define N2_ID 2
#define CO2_ID 3
#define H2O_ID 4

#define NUM_GAS_ELEMENTS 6  // Elements in the homogeneous phase
#define NUM_GAS_SPECIES 5   // Species in the homogeneous phase
#define NUM_GAS_REACTIONS 0 // Reactions in the homogeneous phase

#define SITE_DENSITY 0.000000E+00 // mol/cm^2

#define NUM_SURFACE_ELEMENTS 0  // Additional elements in heterogeneous phase
#define NUM_SURFACE_SPECIES 0   // Species in the heterogeneous phase
#define NUM_SURFACE_REACTIONS 0 // Reactions in the heterogeneous phase

#define NUM_ELEMENTS (NUM_GAS_ELEMENTS + NUM_SURFACE_ELEMENTS)
#define NUM_SPECIES (NUM_GAS_SPECIES + NUM_SURFACE_SPECIES)
#define NUM_REACTIONS (NUM_GAS_REACTIONS + NUM_SURFACE_REACTIONS)

#define NUM_IONS 0

#define NUM_FIT 4

//  ALWAYS on CPU stuff -- can have different def depending on if we are CPU or
//  GPU based. Defined in mechanism.cpp
void atomicWeight(amrex::Real* awt);
//  MISC
void CKAWT(amrex::Real* awt);
void CKNCF(int* ncf);
void CKSYME_STR(amrex::Vector<std::string>& ename);
void CKSYMS_STR(amrex::Vector<std::string>& kname);
void GET_RMAP(int* _rmap);
void CKINU(const int i, int& nspec, int* ki, int* nu);
void CKKFKR(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real* x,
  amrex::Real* q_f,
  amrex::Real* q_r);
void progressRateFR(
  amrex::Real* q_f, amrex::Real* q_r, amrex::Real* sc, amrex::Real T);
//  SPARSE INFORMATION
void SPARSITY_INFO(int* nJdata, const int* consP, int NCELLS);
void SPARSITY_INFO_SYST(int* nJdata, const int* consP, int NCELLS);
void SPARSITY_INFO_SYST_SIMPLIFIED(int* nJdata, const int* consP);
void
SPARSITY_PREPROC_CSC(int* rowVals, int* colPtrs, const int* consP, int NCELLS);
void SPARSITY_PREPROC_CSR(
  int* colVals, int* rowPtrs, const int* consP, int NCELLS, int base);
void SPARSITY_PREPROC_SYST_CSR(
  int* colVals, int* rowPtrs, const int* consP, int NCELLS, int base);
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(
  int* rowVals, int* colPtrs, int* indx, const int* consP);
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(
  int* colVals, int* rowPtr, const int* consP, int base);

// A few mechanism parameters
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKINDX(int& mm, int& kk, int& ii, int& nfit)
{
  mm = 6;
  kk = 5;
  ii = 0;
  nfit = -1; // Why do you need this anyway ?
}

//  inverse molecular weights
#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT const amrex::Real global_imw[5] = {
  0.0064810074077915, // POSF10325
  0.0312519532470779, // O2
  0.0356964374955379, // N2
  0.0227226249176305, // CO2
  0.0555092978073827, // H2O
};
#endif
const amrex::Real h_global_imw[5] = {
  0.0064810074077915, // POSF10325
  0.0312519532470779, // O2
  0.0356964374955379, // N2
  0.0227226249176305, // CO2
  0.0555092978073827, // H2O
};

//  molecular weights
#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT const amrex::Real global_mw[5] = {
  154.297000, // POSF10325
  31.998000,  // O2
  28.014000,  // N2
  44.009000,  // CO2
  18.015000,  // H2O
};
#endif
const amrex::Real h_global_mw[5] = {
  154.297000, // POSF10325
  31.998000,  // O2
  28.014000,  // N2
  44.009000,  // CO2
  18.015000,  // H2O
};

//  inverse molecular weights
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
get_imw(amrex::Real* imw_new)
{
  imw_new[0] = 0.0064810074077915; // POSF10325
  imw_new[1] = 0.0312519532470779; // O2
  imw_new[2] = 0.0356964374955379; // N2
  imw_new[3] = 0.0227226249176305; // CO2
  imw_new[4] = 0.0555092978073827; // H2O
}

//  inverse molecular weight
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
imw(const int n)
{
#if AMREX_DEVICE_COMPILE
  return global_imw[n];
#else
  return h_global_imw[n];
#endif
}
//  molecular weights
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
get_mw(amrex::Real* mw_new)
{
  mw_new[0] = 154.297000; // POSF10325
  mw_new[1] = 31.998000;  // O2
  mw_new[2] = 28.014000;  // N2
  mw_new[3] = 44.009000;  // CO2
  mw_new[4] = 18.015000;  // H2O
}

//  molecular weight
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
mw(const int n)
{
#if AMREX_DEVICE_COMPILE
  return global_mw[n];
#else
  return h_global_mw[n];
#endif
}

//  Returns R, Rc, Patm
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKRP(amrex::Real& ru, amrex::Real& ruc, amrex::Real& pa)
{
  ru = 8.31446261815324e+07;
  ruc = 1.98721558317399615845;
  pa = 1.01325e+06;
}

// compute Cv/R at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
cv_R(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: POSF10325
    species[0] = +1.57855180e+00 + 7.26248320e-02 * T + 7.32363700e-05 * T2 -
                 1.51154440e-07 * T3 + 6.89355790e-11 * T4;
    // species 1: O2
    species[1] = +2.78245636e+00 - 2.99673416e-03 * T + 9.84730201e-06 * T2 -
                 9.68129509e-09 * T3 + 3.24372837e-12 * T4;
    // species 2: N2
    species[2] = +2.29867700e+00 + 1.40824040e-03 * T - 3.96322200e-06 * T2 +
                 5.64151500e-09 * T3 - 2.44485400e-12 * T4;
    // species 3: CO2
    species[3] = +1.35677352e+00 + 8.98459677e-03 * T - 7.12356269e-06 * T2 +
                 2.45919022e-09 * T3 - 1.43699548e-13 * T4;
    // species 4: H2O
    species[4] = +3.19864056e+00 - 2.03643410e-03 * T + 6.52040211e-06 * T2 -
                 5.48797062e-09 * T3 + 1.77197817e-12 * T4;
  } else {
    // species 0: POSF10325
    species[0] = +2.03578910e+01 + 6.65890870e-02 * T - 2.59674690e-05 * T2 +
                 4.47750010e-09 * T3 - 2.36071520e-13 * T4;
    // species 1: O2
    species[1] = +2.28253784e+00 + 1.48308754e-03 * T - 7.57966669e-07 * T2 +
                 2.09470555e-10 * T3 - 2.16717794e-14 * T4;
    // species 2: N2
    species[2] = +1.92664000e+00 + 1.48797680e-03 * T - 5.68476000e-07 * T2 +
                 1.00970380e-10 * T3 - 6.75335100e-15 * T4;
    // species 3: CO2
    species[3] = +2.85746029e+00 + 4.41437026e-03 * T - 2.21481404e-06 * T2 +
                 5.23490188e-10 * T3 - 4.72084164e-14 * T4;
    // species 4: H2O
    species[4] = +2.03399249e+00 + 2.17691804e-03 * T - 1.64072518e-07 * T2 -
                 9.70419870e-11 * T3 + 1.68200992e-14 * T4;
  }
}

// compute Cp/R at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
cp_R(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: POSF10325
    species[0] = +2.57855180e+00 + 7.26248320e-02 * T + 7.32363700e-05 * T2 -
                 1.51154440e-07 * T3 + 6.89355790e-11 * T4;
    // species 1: O2
    species[1] = +3.78245636e+00 - 2.99673416e-03 * T + 9.84730201e-06 * T2 -
                 9.68129509e-09 * T3 + 3.24372837e-12 * T4;
    // species 2: N2
    species[2] = +3.29867700e+00 + 1.40824040e-03 * T - 3.96322200e-06 * T2 +
                 5.64151500e-09 * T3 - 2.44485400e-12 * T4;
    // species 3: CO2
    species[3] = +2.35677352e+00 + 8.98459677e-03 * T - 7.12356269e-06 * T2 +
                 2.45919022e-09 * T3 - 1.43699548e-13 * T4;
    // species 4: H2O
    species[4] = +4.19864056e+00 - 2.03643410e-03 * T + 6.52040211e-06 * T2 -
                 5.48797062e-09 * T3 + 1.77197817e-12 * T4;
  } else {
    // species 0: POSF10325
    species[0] = +2.13578910e+01 + 6.65890870e-02 * T - 2.59674690e-05 * T2 +
                 4.47750010e-09 * T3 - 2.36071520e-13 * T4;
    // species 1: O2
    species[1] = +3.28253784e+00 + 1.48308754e-03 * T - 7.57966669e-07 * T2 +
                 2.09470555e-10 * T3 - 2.16717794e-14 * T4;
    // species 2: N2
    species[2] = +2.92664000e+00 + 1.48797680e-03 * T - 5.68476000e-07 * T2 +
                 1.00970380e-10 * T3 - 6.75335100e-15 * T4;
    // species 3: CO2
    species[3] = +3.85746029e+00 + 4.41437026e-03 * T - 2.21481404e-06 * T2 +
                 5.23490188e-10 * T3 - 4.72084164e-14 * T4;
    // species 4: H2O
    species[4] = +3.03399249e+00 + 2.17691804e-03 * T - 1.64072518e-07 * T2 -
                 9.70419870e-11 * T3 + 1.68200992e-14 * T4;
  }
}

// compute the g/(RT) at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
gibbs(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;
  const amrex::Real invT = 1.0 / T;
  const amrex::Real logT = log(T);

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: POSF10325
    species[0] = -3.802985500000000e+04 * invT - 1.995518320000000e+01 -
                 2.578551800000000e+00 * logT - 3.631241600000000e-02 * T -
                 1.220606166666667e-05 * T2 + 1.259620333333333e-08 * T3 -
                 3.446778950000000e-12 * T4;
    // species 1: O2
    species[1] = -1.063943560000000e+03 * invT + 1.247806300000001e-01 -
                 3.782456360000000e+00 * logT + 1.498367080000000e-03 * T -
                 1.641217001666667e-06 * T2 + 8.067745908333334e-10 * T3 -
                 1.621864185000000e-13 * T4;
    // species 2: N2
    species[2] = -1.020899900000000e+03 * invT - 6.516950000000001e-01 -
                 3.298677000000000e+00 * logT - 7.041202000000000e-04 * T +
                 6.605369999999999e-07 * T2 - 4.701262500000001e-10 * T3 +
                 1.222427000000000e-13 * T4;
    // species 3: CO2
    species[3] = -4.837196970000000e+04 * invT - 7.544278700000000e+00 -
                 2.356773520000000e+00 * logT - 4.492298385000000e-03 * T +
                 1.187260448333333e-06 * T2 - 2.049325183333333e-10 * T3 +
                 7.184977399999999e-15 * T4;
    // species 4: H2O
    species[4] = -3.029372670000000e+04 * invT + 5.047672768000000e+00 -
                 4.198640560000000e+00 * logT + 1.018217050000000e-03 * T -
                 1.086733685000000e-06 * T2 + 4.573308850000000e-10 * T3 -
                 8.859890850000000e-14 * T4;
  } else {
    // species 0: POSF10325
    species[0] = -4.579699600000000e+04 * invT + 1.074939770000000e+02 -
                 2.135789100000000e+01 * logT - 3.329454350000000e-02 * T +
                 4.327911500000000e-06 * T2 - 3.731250083333333e-10 * T3 +
                 1.180357600000000e-14 * T4;
    // species 1: O2
    species[1] = -1.088457720000000e+03 * invT - 2.170693450000000e+00 -
                 3.282537840000000e+00 * logT - 7.415437700000000e-04 * T +
                 1.263277781666667e-07 * T2 - 1.745587958333333e-11 * T3 +
                 1.083588970000000e-15 * T4;
    // species 2: N2
    species[2] = -9.227977000000000e+02 * invT - 3.053888000000000e+00 -
                 2.926640000000000e+00 * logT - 7.439884000000000e-04 * T +
                 9.474600000000001e-08 * T2 - 8.414198333333333e-12 * T3 +
                 3.376675500000000e-16 * T4;
    // species 3: CO2
    species[3] = -4.875916600000000e+04 * invT + 1.585822230000000e+00 -
                 3.857460290000000e+00 * logT - 2.207185130000000e-03 * T +
                 3.691356733333334e-07 * T2 - 4.362418233333334e-11 * T3 +
                 2.360420820000000e-15 * T4;
    // species 4: H2O
    species[4] = -3.000429710000000e+04 * invT - 1.932777610000000e+00 -
                 3.033992490000000e+00 * logT - 1.088459020000000e-03 * T +
                 2.734541966666666e-08 * T2 + 8.086832250000000e-12 * T3 -
                 8.410049600000000e-16 * T4;
  }
}

// compute the a/(RT) at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
helmholtz(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;
  const amrex::Real invT = 1.0 / T;
  const amrex::Real logT = log(T);

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: POSF10325
    species[0] = -3.80298550e+04 * invT - 2.09551832e+01 -
                 2.57855180e+00 * logT - 3.63124160e-02 * T -
                 1.22060617e-05 * T2 + 1.25962033e-08 * T3 -
                 3.44677895e-12 * T4;
    // species 1: O2
    species[1] = -1.06394356e+03 * invT - 8.75219370e-01 -
                 3.78245636e+00 * logT + 1.49836708e-03 * T -
                 1.64121700e-06 * T2 + 8.06774591e-10 * T3 -
                 1.62186418e-13 * T4;
    // species 2: N2
    species[2] = -1.02089990e+03 * invT - 1.65169500e+00 -
                 3.29867700e+00 * logT - 7.04120200e-04 * T +
                 6.60537000e-07 * T2 - 4.70126250e-10 * T3 +
                 1.22242700e-13 * T4;
    // species 3: CO2
    species[3] = -4.83719697e+04 * invT - 8.54427870e+00 -
                 2.35677352e+00 * logT - 4.49229839e-03 * T +
                 1.18726045e-06 * T2 - 2.04932518e-10 * T3 +
                 7.18497740e-15 * T4;
    // species 4: H2O
    species[4] = -3.02937267e+04 * invT + 4.04767277e+00 -
                 4.19864056e+00 * logT + 1.01821705e-03 * T -
                 1.08673369e-06 * T2 + 4.57330885e-10 * T3 -
                 8.85989085e-14 * T4;
  } else {
    // species 0: POSF10325
    species[0] = -4.57969960e+04 * invT + 1.06493977e+02 -
                 2.13578910e+01 * logT - 3.32945435e-02 * T +
                 4.32791150e-06 * T2 - 3.73125008e-10 * T3 +
                 1.18035760e-14 * T4;
    // species 1: O2
    species[1] = -1.08845772e+03 * invT - 3.17069345e+00 -
                 3.28253784e+00 * logT - 7.41543770e-04 * T +
                 1.26327778e-07 * T2 - 1.74558796e-11 * T3 +
                 1.08358897e-15 * T4;
    // species 2: N2
    species[2] = -9.22797700e+02 * invT - 4.05388800e+00 -
                 2.92664000e+00 * logT - 7.43988400e-04 * T +
                 9.47460000e-08 * T2 - 8.41419833e-12 * T3 +
                 3.37667550e-16 * T4;
    // species 3: CO2
    species[3] = -4.87591660e+04 * invT + 5.85822230e-01 -
                 3.85746029e+00 * logT - 2.20718513e-03 * T +
                 3.69135673e-07 * T2 - 4.36241823e-11 * T3 +
                 2.36042082e-15 * T4;
    // species 4: H2O
    species[4] = -3.00042971e+04 * invT - 2.93277761e+00 -
                 3.03399249e+00 * logT - 1.08845902e-03 * T +
                 2.73454197e-08 * T2 + 8.08683225e-12 * T3 -
                 8.41004960e-16 * T4;
  }
}

// compute the e/(RT) at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
speciesInternalEnergy(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;
  const amrex::Real invT = 1.0 / T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: POSF10325
    species[0] = +1.57855180e+00 + 3.63124160e-02 * T + 2.44121233e-05 * T2 -
                 3.77886100e-08 * T3 + 1.37871158e-11 * T4 -
                 3.80298550e+04 * invT;
    // species 1: O2
    species[1] = +2.78245636e+00 - 1.49836708e-03 * T + 3.28243400e-06 * T2 -
                 2.42032377e-09 * T3 + 6.48745674e-13 * T4 -
                 1.06394356e+03 * invT;
    // species 2: N2
    species[2] = +2.29867700e+00 + 7.04120200e-04 * T - 1.32107400e-06 * T2 +
                 1.41037875e-09 * T3 - 4.88970800e-13 * T4 -
                 1.02089990e+03 * invT;
    // species 3: CO2
    species[3] = +1.35677352e+00 + 4.49229839e-03 * T - 2.37452090e-06 * T2 +
                 6.14797555e-10 * T3 - 2.87399096e-14 * T4 -
                 4.83719697e+04 * invT;
    // species 4: H2O
    species[4] = +3.19864056e+00 - 1.01821705e-03 * T + 2.17346737e-06 * T2 -
                 1.37199266e-09 * T3 + 3.54395634e-13 * T4 -
                 3.02937267e+04 * invT;
  } else {
    // species 0: POSF10325
    species[0] = +2.03578910e+01 + 3.32945435e-02 * T - 8.65582300e-06 * T2 +
                 1.11937502e-09 * T3 - 4.72143040e-14 * T4 -
                 4.57969960e+04 * invT;
    // species 1: O2
    species[1] = +2.28253784e+00 + 7.41543770e-04 * T - 2.52655556e-07 * T2 +
                 5.23676387e-11 * T3 - 4.33435588e-15 * T4 -
                 1.08845772e+03 * invT;
    // species 2: N2
    species[2] = +1.92664000e+00 + 7.43988400e-04 * T - 1.89492000e-07 * T2 +
                 2.52425950e-11 * T3 - 1.35067020e-15 * T4 -
                 9.22797700e+02 * invT;
    // species 3: CO2
    species[3] = +2.85746029e+00 + 2.20718513e-03 * T - 7.38271347e-07 * T2 +
                 1.30872547e-10 * T3 - 9.44168328e-15 * T4 -
                 4.87591660e+04 * invT;
    // species 4: H2O
    species[4] = +2.03399249e+00 + 1.08845902e-03 * T - 5.46908393e-08 * T2 -
                 2.42604967e-11 * T3 + 3.36401984e-15 * T4 -
                 3.00042971e+04 * invT;
  }
}

// compute the h/(RT) at the given temperature (Eq 20)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
speciesEnthalpy(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;
  const amrex::Real invT = 1.0 / T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: POSF10325
    species[0] = +2.57855180e+00 + 3.63124160e-02 * T + 2.44121233e-05 * T2 -
                 3.77886100e-08 * T3 + 1.37871158e-11 * T4 -
                 3.80298550e+04 * invT;
    // species 1: O2
    species[1] = +3.78245636e+00 - 1.49836708e-03 * T + 3.28243400e-06 * T2 -
                 2.42032377e-09 * T3 + 6.48745674e-13 * T4 -
                 1.06394356e+03 * invT;
    // species 2: N2
    species[2] = +3.29867700e+00 + 7.04120200e-04 * T - 1.32107400e-06 * T2 +
                 1.41037875e-09 * T3 - 4.88970800e-13 * T4 -
                 1.02089990e+03 * invT;
    // species 3: CO2
    species[3] = +2.35677352e+00 + 4.49229839e-03 * T - 2.37452090e-06 * T2 +
                 6.14797555e-10 * T3 - 2.87399096e-14 * T4 -
                 4.83719697e+04 * invT;
    // species 4: H2O
    species[4] = +4.19864056e+00 - 1.01821705e-03 * T + 2.17346737e-06 * T2 -
                 1.37199266e-09 * T3 + 3.54395634e-13 * T4 -
                 3.02937267e+04 * invT;
  } else {
    // species 0: POSF10325
    species[0] = +2.13578910e+01 + 3.32945435e-02 * T - 8.65582300e-06 * T2 +
                 1.11937502e-09 * T3 - 4.72143040e-14 * T4 -
                 4.57969960e+04 * invT;
    // species 1: O2
    species[1] = +3.28253784e+00 + 7.41543770e-04 * T - 2.52655556e-07 * T2 +
                 5.23676387e-11 * T3 - 4.33435588e-15 * T4 -
                 1.08845772e+03 * invT;
    // species 2: N2
    species[2] = +2.92664000e+00 + 7.43988400e-04 * T - 1.89492000e-07 * T2 +
                 2.52425950e-11 * T3 - 1.35067020e-15 * T4 -
                 9.22797700e+02 * invT;
    // species 3: CO2
    species[3] = +3.85746029e+00 + 2.20718513e-03 * T - 7.38271347e-07 * T2 +
                 1.30872547e-10 * T3 - 9.44168328e-15 * T4 -
                 4.87591660e+04 * invT;
    // species 4: H2O
    species[4] = +3.03399249e+00 + 1.08845902e-03 * T - 5.46908393e-08 * T2 -
                 2.42604967e-11 * T3 + 3.36401984e-15 * T4 -
                 3.00042971e+04 * invT;
  }
}

// compute the S/R at the given temperature (Eq 21)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
speciesEntropy(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;
  const amrex::Real logT = log(T);

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: POSF10325
    species[0] = +2.57855180e+00 * logT + 7.26248320e-02 * T +
                 3.66181850e-05 * T2 - 5.03848133e-08 * T3 +
                 1.72338947e-11 * T4 + 2.25337350e+01;
    // species 1: O2
    species[1] = +3.78245636e+00 * logT - 2.99673416e-03 * T +
                 4.92365101e-06 * T2 - 3.22709836e-09 * T3 +
                 8.10932092e-13 * T4 + 3.65767573e+00;
    // species 2: N2
    species[2] = +3.29867700e+00 * logT + 1.40824040e-03 * T -
                 1.98161100e-06 * T2 + 1.88050500e-09 * T3 -
                 6.11213500e-13 * T4 + 3.95037200e+00;
    // species 3: CO2
    species[3] = +2.35677352e+00 * logT + 8.98459677e-03 * T -
                 3.56178134e-06 * T2 + 8.19730073e-10 * T3 -
                 3.59248870e-14 * T4 + 9.90105222e+00;
    // species 4: H2O
    species[4] = +4.19864056e+00 * logT - 2.03643410e-03 * T +
                 3.26020105e-06 * T2 - 1.82932354e-09 * T3 +
                 4.42994543e-13 * T4 - 8.49032208e-01;
  } else {
    // species 0: POSF10325
    species[0] = +2.13578910e+01 * logT + 6.65890870e-02 * T -
                 1.29837345e-05 * T2 + 1.49250003e-09 * T3 -
                 5.90178800e-14 * T4 - 8.61360860e+01;
    // species 1: O2
    species[1] = +3.28253784e+00 * logT + 1.48308754e-03 * T -
                 3.78983334e-07 * T2 + 6.98235183e-11 * T3 -
                 5.41794485e-15 * T4 + 5.45323129e+00;
    // species 2: N2
    species[2] = +2.92664000e+00 * logT + 1.48797680e-03 * T -
                 2.84238000e-07 * T2 + 3.36567933e-11 * T3 -
                 1.68833775e-15 * T4 + 5.98052800e+00;
    // species 3: CO2
    species[3] = +3.85746029e+00 * logT + 4.41437026e-03 * T -
                 1.10740702e-06 * T2 + 1.74496729e-10 * T3 -
                 1.18021041e-14 * T4 + 2.27163806e+00;
    // species 4: H2O
    species[4] = +3.03399249e+00 * logT + 2.17691804e-03 * T -
                 8.20362590e-08 * T2 - 3.23473290e-11 * T3 +
                 4.20502480e-15 * T4 + 4.96677010e+00;
  }
}

// compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
dcvpRdT(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: POSF10325
    species[0] = +7.26248320e-02 + 1.46472740e-04 * T - 4.53463320e-07 * T2 +
                 2.75742316e-10 * T3;
    // species 1: O2
    species[1] = -2.99673416e-03 + 1.96946040e-05 * T - 2.90438853e-08 * T2 +
                 1.29749135e-11 * T3;
    // species 2: N2
    species[2] = +1.40824040e-03 - 7.92644400e-06 * T + 1.69245450e-08 * T2 -
                 9.77941600e-12 * T3;
    // species 3: CO2
    species[3] = +8.98459677e-03 - 1.42471254e-05 * T + 7.37757066e-09 * T2 -
                 5.74798192e-13 * T3;
    // species 4: H2O
    species[4] = -2.03643410e-03 + 1.30408042e-05 * T - 1.64639119e-08 * T2 +
                 7.08791268e-12 * T3;
  } else {
    // species 0: POSF10325
    species[0] = +6.65890870e-02 - 5.19349380e-05 * T + 1.34325003e-08 * T2 -
                 9.44286080e-13 * T3;
    // species 1: O2
    species[1] = +1.48308754e-03 - 1.51593334e-06 * T + 6.28411665e-10 * T2 -
                 8.66871176e-14 * T3;
    // species 2: N2
    species[2] = +1.48797680e-03 - 1.13695200e-06 * T + 3.02911140e-10 * T2 -
                 2.70134040e-14 * T3;
    // species 3: CO2
    species[3] = +4.41437026e-03 - 4.42962808e-06 * T + 1.57047056e-09 * T2 -
                 1.88833666e-13 * T3;
    // species 4: H2O
    species[4] = +2.17691804e-03 - 3.28145036e-07 * T - 2.91125961e-10 * T2 +
                 6.72803968e-14 * T3;
  }
}

// Returns the mean specific heat at CP (Eq. 33)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPBL(const amrex::Real T, const amrex::Real x[], amrex::Real& cpbl)
{
  amrex::Real result = 0;
  amrex::Real cpor[5]; //  temporary storage
  cp_R(cpor, T);

  // perform dot product
  for (int id = 0; id < 5; ++id) {
    result += x[id] * cpor[id];
  }

  cpbl = result * 8.31446261815324e+07;
}

// Returns the mean specific heat at CP (Eq. 34)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPBS(const amrex::Real T, const amrex::Real y[], amrex::Real& cpbs)
{
  amrex::Real result = 0.0;

  // compute Cp/R at the given temperature
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: POSF10325
    result += y[0] *
              (+2.57855180e+00 + 7.26248320e-02 * T + 7.32363700e-05 * T2 -
               1.51154440e-07 * T3 + 6.89355790e-11 * T4) *
              0.0064810074077915;
    // species 1: O2
    result += y[1] *
              (+3.78245636e+00 - 2.99673416e-03 * T + 9.84730201e-06 * T2 -
               9.68129509e-09 * T3 + 3.24372837e-12 * T4) *
              0.0312519532470779;
    // species 2: N2
    result += y[2] *
              (+3.29867700e+00 + 1.40824040e-03 * T - 3.96322200e-06 * T2 +
               5.64151500e-09 * T3 - 2.44485400e-12 * T4) *
              0.0356964374955379;
    // species 3: CO2
    result += y[3] *
              (+2.35677352e+00 + 8.98459677e-03 * T - 7.12356269e-06 * T2 +
               2.45919022e-09 * T3 - 1.43699548e-13 * T4) *
              0.0227226249176305;
    // species 4: H2O
    result += y[4] *
              (+4.19864056e+00 - 2.03643410e-03 * T + 6.52040211e-06 * T2 -
               5.48797062e-09 * T3 + 1.77197817e-12 * T4) *
              0.0555092978073827;
  } else {
    // species 0: POSF10325
    result += y[0] *
              (+2.13578910e+01 + 6.65890870e-02 * T - 2.59674690e-05 * T2 +
               4.47750010e-09 * T3 - 2.36071520e-13 * T4) *
              0.0064810074077915;
    // species 1: O2
    result += y[1] *
              (+3.28253784e+00 + 1.48308754e-03 * T - 7.57966669e-07 * T2 +
               2.09470555e-10 * T3 - 2.16717794e-14 * T4) *
              0.0312519532470779;
    // species 2: N2
    result += y[2] *
              (+2.92664000e+00 + 1.48797680e-03 * T - 5.68476000e-07 * T2 +
               1.00970380e-10 * T3 - 6.75335100e-15 * T4) *
              0.0356964374955379;
    // species 3: CO2
    result += y[3] *
              (+3.85746029e+00 + 4.41437026e-03 * T - 2.21481404e-06 * T2 +
               5.23490188e-10 * T3 - 4.72084164e-14 * T4) *
              0.0227226249176305;
    // species 4: H2O
    result += y[4] *
              (+3.03399249e+00 + 2.17691804e-03 * T - 1.64072518e-07 * T2 -
               9.70419870e-11 * T3 + 1.68200992e-14 * T4) *
              0.0555092978073827;
  }

  cpbs = result * 8.31446261815324e+07;
}

// Returns the mean specific heat at CV (Eq. 35)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCVBL(const amrex::Real T, const amrex::Real x[], amrex::Real& cvbl)
{
  amrex::Real result = 0;
  amrex::Real cvor[5]; //  temporary storage
  cv_R(cvor, T);

  // perform dot product
  for (int id = 0; id < 5; ++id) {
    result += x[id] * cvor[id];
  }

  cvbl = result * 8.31446261815324e+07;
}

// Returns the mean specific heat at CV (Eq. 36)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCVBS(const amrex::Real T, const amrex::Real y[], amrex::Real& cvbs)
{
  amrex::Real result = 0.0;
  // compute Cv/R at the given temperature
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: POSF10325
    result += y[0] *
              (+1.57855180e+00 + 7.26248320e-02 * T + 7.32363700e-05 * T2 -
               1.51154440e-07 * T3 + 6.89355790e-11 * T4) *
              0.0064810074077915;
    // species 1: O2
    result += y[1] *
              (+2.78245636e+00 - 2.99673416e-03 * T + 9.84730201e-06 * T2 -
               9.68129509e-09 * T3 + 3.24372837e-12 * T4) *
              0.0312519532470779;
    // species 2: N2
    result += y[2] *
              (+2.29867700e+00 + 1.40824040e-03 * T - 3.96322200e-06 * T2 +
               5.64151500e-09 * T3 - 2.44485400e-12 * T4) *
              0.0356964374955379;
    // species 3: CO2
    result += y[3] *
              (+1.35677352e+00 + 8.98459677e-03 * T - 7.12356269e-06 * T2 +
               2.45919022e-09 * T3 - 1.43699548e-13 * T4) *
              0.0227226249176305;
    // species 4: H2O
    result += y[4] *
              (+3.19864056e+00 - 2.03643410e-03 * T + 6.52040211e-06 * T2 -
               5.48797062e-09 * T3 + 1.77197817e-12 * T4) *
              0.0555092978073827;
  } else {
    // species 0: POSF10325
    result += y[0] *
              (+2.03578910e+01 + 6.65890870e-02 * T - 2.59674690e-05 * T2 +
               4.47750010e-09 * T3 - 2.36071520e-13 * T4) *
              0.0064810074077915;
    // species 1: O2
    result += y[1] *
              (+2.28253784e+00 + 1.48308754e-03 * T - 7.57966669e-07 * T2 +
               2.09470555e-10 * T3 - 2.16717794e-14 * T4) *
              0.0312519532470779;
    // species 2: N2
    result += y[2] *
              (+1.92664000e+00 + 1.48797680e-03 * T - 5.68476000e-07 * T2 +
               1.00970380e-10 * T3 - 6.75335100e-15 * T4) *
              0.0356964374955379;
    // species 3: CO2
    result += y[3] *
              (+2.85746029e+00 + 4.41437026e-03 * T - 2.21481404e-06 * T2 +
               5.23490188e-10 * T3 - 4.72084164e-14 * T4) *
              0.0227226249176305;
    // species 4: H2O
    result += y[4] *
              (+2.03399249e+00 + 2.17691804e-03 * T - 1.64072518e-07 * T2 -
               9.70419870e-11 * T3 + 1.68200992e-14 * T4) *
              0.0555092978073827;
  }

  cvbs = result * 8.31446261815324e+07;
}

// Returns the mean enthalpy of the mixture in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHBML(const amrex::Real T, const amrex::Real x[], amrex::Real& hbml)
{
  amrex::Real result = 0;
  amrex::Real hml[5];                        //  temporary storage
  amrex::Real RT = 8.31446261815324e+07 * T; // R*T
  speciesEnthalpy(hml, T);

  // perform dot product
  for (int id = 0; id < 5; ++id) {
    result += x[id] * hml[id];
  }

  hbml = result * RT;
}

// Returns mean enthalpy of mixture in mass units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHBMS(const amrex::Real T, const amrex::Real y[], amrex::Real& hbms)
{
  amrex::Real result = 0.0;
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;
  const amrex::Real invT = 1.0 / T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: POSF10325
    result +=
      y[0] *
      (+2.57855180e+00 + 3.63124160e-02 * T + 2.44121233e-05 * T2 -
       3.77886100e-08 * T3 + 1.37871158e-11 * T4 - 3.80298550e+04 * invT) *
      0.0064810074077915;
    // species 1: O2
    result +=
      y[1] *
      (+3.78245636e+00 - 1.49836708e-03 * T + 3.28243400e-06 * T2 -
       2.42032377e-09 * T3 + 6.48745674e-13 * T4 - 1.06394356e+03 * invT) *
      0.0312519532470779;
    // species 2: N2
    result +=
      y[2] *
      (+3.29867700e+00 + 7.04120200e-04 * T - 1.32107400e-06 * T2 +
       1.41037875e-09 * T3 - 4.88970800e-13 * T4 - 1.02089990e+03 * invT) *
      0.0356964374955379;
    // species 3: CO2
    result +=
      y[3] *
      (+2.35677352e+00 + 4.49229839e-03 * T - 2.37452090e-06 * T2 +
       6.14797555e-10 * T3 - 2.87399096e-14 * T4 - 4.83719697e+04 * invT) *
      0.0227226249176305;
    // species 4: H2O
    result +=
      y[4] *
      (+4.19864056e+00 - 1.01821705e-03 * T + 2.17346737e-06 * T2 -
       1.37199266e-09 * T3 + 3.54395634e-13 * T4 - 3.02937267e+04 * invT) *
      0.0555092978073827;
  } else {
    // species 0: POSF10325
    result +=
      y[0] *
      (+2.13578910e+01 + 3.32945435e-02 * T - 8.65582300e-06 * T2 +
       1.11937502e-09 * T3 - 4.72143040e-14 * T4 - 4.57969960e+04 * invT) *
      0.0064810074077915;
    // species 1: O2
    result +=
      y[1] *
      (+3.28253784e+00 + 7.41543770e-04 * T - 2.52655556e-07 * T2 +
       5.23676387e-11 * T3 - 4.33435588e-15 * T4 - 1.08845772e+03 * invT) *
      0.0312519532470779;
    // species 2: N2
    result +=
      y[2] *
      (+2.92664000e+00 + 7.43988400e-04 * T - 1.89492000e-07 * T2 +
       2.52425950e-11 * T3 - 1.35067020e-15 * T4 - 9.22797700e+02 * invT) *
      0.0356964374955379;
    // species 3: CO2
    result +=
      y[3] *
      (+3.85746029e+00 + 2.20718513e-03 * T - 7.38271347e-07 * T2 +
       1.30872547e-10 * T3 - 9.44168328e-15 * T4 - 4.87591660e+04 * invT) *
      0.0227226249176305;
    // species 4: H2O
    result +=
      y[4] *
      (+3.03399249e+00 + 1.08845902e-03 * T - 5.46908393e-08 * T2 -
       2.42604967e-11 * T3 + 3.36401984e-15 * T4 - 3.00042971e+04 * invT) *
      0.0555092978073827;
  }

  const amrex::Real RT = 8.31446261815324e+07 * T; // R*T

  hbms = result * RT;
}

// get mean internal energy in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUBML(const amrex::Real T, const amrex::Real x[], amrex::Real& ubml)
{
  amrex::Real result = 0;
  amrex::Real uml[5];                        //  temporary energy array
  amrex::Real RT = 8.31446261815324e+07 * T; // R*T
  speciesInternalEnergy(uml, T);

  // perform dot product
  for (int id = 0; id < 5; ++id) {
    result += x[id] * uml[id];
  }

  ubml = result * RT;
}

// get mean internal energy in mass units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUBMS(const amrex::Real T, const amrex::Real y[], amrex::Real& ubms)
{
  amrex::Real result = 0.0;
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;
  const amrex::Real invT = 1.0 / T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: POSF10325
    result +=
      y[0] *
      (+1.57855180e+00 + 3.63124160e-02 * T + 2.44121233e-05 * T2 -
       3.77886100e-08 * T3 + 1.37871158e-11 * T4 - 3.80298550e+04 * invT) *
      0.0064810074077915;
    // species 1: O2
    result +=
      y[1] *
      (+2.78245636e+00 - 1.49836708e-03 * T + 3.28243400e-06 * T2 -
       2.42032377e-09 * T3 + 6.48745674e-13 * T4 - 1.06394356e+03 * invT) *
      0.0312519532470779;
    // species 2: N2
    result +=
      y[2] *
      (+2.29867700e+00 + 7.04120200e-04 * T - 1.32107400e-06 * T2 +
       1.41037875e-09 * T3 - 4.88970800e-13 * T4 - 1.02089990e+03 * invT) *
      0.0356964374955379;
    // species 3: CO2
    result +=
      y[3] *
      (+1.35677352e+00 + 4.49229839e-03 * T - 2.37452090e-06 * T2 +
       6.14797555e-10 * T3 - 2.87399096e-14 * T4 - 4.83719697e+04 * invT) *
      0.0227226249176305;
    // species 4: H2O
    result +=
      y[4] *
      (+3.19864056e+00 - 1.01821705e-03 * T + 2.17346737e-06 * T2 -
       1.37199266e-09 * T3 + 3.54395634e-13 * T4 - 3.02937267e+04 * invT) *
      0.0555092978073827;
  } else {
    // species 0: POSF10325
    result +=
      y[0] *
      (+2.03578910e+01 + 3.32945435e-02 * T - 8.65582300e-06 * T2 +
       1.11937502e-09 * T3 - 4.72143040e-14 * T4 - 4.57969960e+04 * invT) *
      0.0064810074077915;
    // species 1: O2
    result +=
      y[1] *
      (+2.28253784e+00 + 7.41543770e-04 * T - 2.52655556e-07 * T2 +
       5.23676387e-11 * T3 - 4.33435588e-15 * T4 - 1.08845772e+03 * invT) *
      0.0312519532470779;
    // species 2: N2
    result +=
      y[2] *
      (+1.92664000e+00 + 7.43988400e-04 * T - 1.89492000e-07 * T2 +
       2.52425950e-11 * T3 - 1.35067020e-15 * T4 - 9.22797700e+02 * invT) *
      0.0356964374955379;
    // species 3: CO2
    result +=
      y[3] *
      (+2.85746029e+00 + 2.20718513e-03 * T - 7.38271347e-07 * T2 +
       1.30872547e-10 * T3 - 9.44168328e-15 * T4 - 4.87591660e+04 * invT) *
      0.0227226249176305;
    // species 4: H2O
    result +=
      y[4] *
      (+2.03399249e+00 + 1.08845902e-03 * T - 5.46908393e-08 * T2 -
       2.42604967e-11 * T3 + 3.36401984e-15 * T4 - 3.00042971e+04 * invT) *
      0.0555092978073827;
  }

  const amrex::Real RT = 8.31446261815324e+07 * T; // R*T

  ubms = result * RT;
}

// get mixture entropy in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKSBML(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real x[],
  amrex::Real& sbml)
{
  amrex::Real result = 0;
  // Log of normalized pressure in cgs units dynes/cm^2 by Patm
  amrex::Real logPratio = log(P / 1013250.0);
  amrex::Real sor[5]; //  temporary storage
  speciesEntropy(sor, T);

  // Compute Eq 42
  for (int id = 0; id < 5; ++id) {
    result += x[id] * (sor[id] - log((x[id] + 1e-100)) - logPratio);
  }

  sbml = result * 8.31446261815324e+07;
}

// get mixture entropy in mass units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKSBMS(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real y[],
  amrex::Real& sbms)
{
  amrex::Real result = 0;
  // Log of normalized pressure in cgs units dynes/cm^2 by Patm
  amrex::Real logPratio = log(P / 1013250.0);
  amrex::Real sor[5];  //  temporary storage
  amrex::Real x[5];    //  need a ytx conversion
  amrex::Real YOW = 0; // See Eq 4, 6 in CK Manual

  // Compute inverse of mean molecular wt first
  for (int i = 0; i < 5; i++) {
    YOW += y[i] * imw(i);
  }
  // Now compute y to x conversion
  x[0] = y[0] / (154.297000 * YOW);
  x[1] = y[1] / (31.998000 * YOW);
  x[2] = y[2] / (28.014000 * YOW);
  x[3] = y[3] / (44.009000 * YOW);
  x[4] = y[4] / (18.015000 * YOW);
  speciesEntropy(sor, T);
  // Perform computation in Eq 42 and 43
  for (int i = 0; i < 5; i++) {
    result += x[i] * (sor[i] - log((x[i] + 1e-100)) - logPratio);
  }
  // Scale by R/W
  sbms = result * 8.31446261815324e+07 * YOW;
}

//  get temperature given internal energy in mass units and mass fracs
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
GET_T_GIVEN_EY(
  const amrex::Real e, const amrex::Real y[], amrex::Real& t, int& ierr)
{
#ifdef CONVERGENCE
  const int maxiter = 5000;
  const amrex::Real tol = 1.e-12;
#else
  const int maxiter = 200;
  const amrex::Real tol = 1.e-6;
#endif
  amrex::Real tmin = 90;   // max lower bound for thermo def
  amrex::Real tmax = 4000; // min upper bound for thermo def
  amrex::Real e1, emin, emax, cv, t1, dt;
  CKUBMS(tmin, y, emin);
  CKUBMS(tmax, y, emax);
  if (e < emin) {
    // Linear Extrapolation below tmin
    CKCVBS(tmin, y, cv);
    t = tmin - (emin - e) / cv;
    ierr = 1;
    return;
  }
  if (e > emax) {
    // Linear Extrapolation above tmax
    CKCVBS(tmax, y, cv);
    t = tmax - (emax - e) / cv;
    ierr = 1;
    return;
  }
  t1 = t;
  if (t1 < tmin || t1 > tmax) {
    t1 = tmin + (tmax - tmin) / (emax - emin) * (e - emin);
  }
  for (int i = 0; i < maxiter; ++i) {
    CKUBMS(t1, y, e1);
    CKCVBS(t1, y, cv);
    dt = (e - e1) / cv;
    if (dt > 100.) {
      dt = 100.;
    } else if (dt < -100.) {
      dt = -100.;
    } else if (fabs(dt) < tol) {
      break;
    }
    t1 += dt;
  }
  t = t1;
  ierr = 0;
}

//  get temperature given enthalpy in mass units and mass fracs
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
GET_T_GIVEN_HY(
  const amrex::Real h, const amrex::Real y[], amrex::Real& t, int& ierr)
{
#ifdef CONVERGENCE
  const int maxiter = 5000;
  const amrex::Real tol = 1.e-12;
#else
  const int maxiter = 200;
  const amrex::Real tol = 1.e-6;
#endif
  amrex::Real tmin = 90;   // max lower bound for thermo def
  amrex::Real tmax = 4000; // min upper bound for thermo def
  amrex::Real h1, hmin, hmax, cp, t1, dt;
  CKHBMS(tmin, y, hmin);
  CKHBMS(tmax, y, hmax);
  if (h < hmin) {
    // Linear Extrapolation below tmin
    CKCPBS(tmin, y, cp);
    t = tmin - (hmin - h) / cp;
    ierr = 1;
    return;
  }
  if (h > hmax) {
    // Linear Extrapolation above tmax
    CKCPBS(tmax, y, cp);
    t = tmax - (hmax - h) / cp;
    ierr = 1;
    return;
  }
  t1 = t;
  if (t1 < tmin || t1 > tmax) {
    t1 = tmin + (tmax - tmin) / (hmax - hmin) * (h - hmin);
  }
  for (int i = 0; i < maxiter; ++i) {
    CKHBMS(t1, y, h1);
    CKCPBS(t1, y, cp);
    dt = (h - h1) / cp;
    if (dt > 100.) {
      dt = 100.;
    } else if (dt < -100.) {
      dt = -100.;
    } else if (fabs(dt) < tol) {
      break;
    }
    t1 += dt;
  }
  t = t1;
  ierr = 0;
}

// Compute P = rhoRT/W(x)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKPX(
  const amrex::Real rho,
  const amrex::Real T,
  const amrex::Real x[],
  amrex::Real& P)
{
  amrex::Real XW = 0;                      //  To hold mean molecular wt
  XW += x[0] * 154.297000;                 // POSF10325
  XW += x[1] * 31.998000;                  // O2
  XW += x[2] * 28.014000;                  // N2
  XW += x[3] * 44.009000;                  // CO2
  XW += x[4] * 18.015000;                  // H2O
  P = rho * 8.31446261815324e+07 * T / XW; // P = rho*R*T/W
}

// Compute P = rhoRT/W(y)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKPY(
  const amrex::Real rho,
  const amrex::Real T,
  const amrex::Real y[],
  amrex::Real& P)
{
  amrex::Real YOW = 0; //  for computing mean MW

  for (int i = 0; i < 5; i++) {
    YOW += y[i] * imw(i);
  }
  P = rho * 8.31446261815324e+07 * T * YOW; // P = rho*R*T/W
}

// Compute P = rhoRT/W(c)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKPC(
  const amrex::Real rho,
  const amrex::Real T,
  const amrex::Real c[],
  amrex::Real& P)
{
  // See Eq 5 in CK Manual
  amrex::Real W = 0;
  amrex::Real sumC = 0;
  W += c[0] * 154.297000; // POSF10325
  W += c[1] * 31.998000;  // O2
  W += c[2] * 28.014000;  // N2
  W += c[3] * 44.009000;  // CO2
  W += c[4] * 18.015000;  // H2O

  for (int id = 0; id < 5; ++id) {
    sumC += c[id];
  }
  P = rho * 8.31446261815324e+07 * T * sumC / W; // P = rho*R*T/W
}

// Compute rho = PW(x)/RT
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKRHOX(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real x[],
  amrex::Real& rho)
{
  amrex::Real XW = 0;                        //  To hold mean molecular wt
  XW += x[0] * 154.297000;                   // POSF10325
  XW += x[1] * 31.998000;                    // O2
  XW += x[2] * 28.014000;                    // N2
  XW += x[3] * 44.009000;                    // CO2
  XW += x[4] * 18.015000;                    // H2O
  rho = P * XW / (8.31446261815324e+07 * T); // rho = P*W/(R*T)
}

// Compute rho = P*W(y)/RT
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKRHOY(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real y[],
  amrex::Real& rho)
{
  amrex::Real YOW = 0;

  for (int i = 0; i < 5; i++) {
    YOW += y[i] * imw(i);
  }

  rho = P / (8.31446261815324e+07 * T * YOW); // rho = P*W/(R*T)
}

// Compute rho = P*W(c)/(R*T)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKRHOC(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real c[],
  amrex::Real& rho)
{
  // See Eq 5 in CK Manual
  amrex::Real W = 0;
  amrex::Real sumC = 0;
  W += c[0] * 154.297000; // POSF10325
  W += c[1] * 31.998000;  // O2
  W += c[2] * 28.014000;  // N2
  W += c[3] * 44.009000;  // CO2
  W += c[4] * 18.015000;  // H2O

  for (int id = 0; id < 5; ++id) {
    sumC += c[id];
  }
  rho = P * W / (sumC * T * 8.31446261815324e+07); // rho = PW/(R*T)
}

// get molecular weight for all species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWT(amrex::Real wt[])
{
  get_mw(wt);
}

// given y[species]: mass fractions
// s mean molecular weight (gm/mole)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKMMWY(const amrex::Real y[], amrex::Real& wtm)
{
  amrex::Real YOW = 0;

  for (int i = 0; i < 5; i++) {
    YOW += y[i] * imw(i);
  }

  wtm = 1.0 / YOW;
}

// given x[species]: mole fractions
// returns mean molecular weight (gm/mole)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKMMWX(const amrex::Real x[], amrex::Real& wtm)
{
  amrex::Real XW = 0;      //  see Eq 4 in CK Manual
  XW += x[0] * 154.297000; // POSF10325
  XW += x[1] * 31.998000;  // O2
  XW += x[2] * 28.014000;  // N2
  XW += x[3] * 44.009000;  // CO2
  XW += x[4] * 18.015000;  // H2O
  wtm = XW;
}

// given c[species]: molar concentration
// returns mean molecular weight (gm/mole)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKMMWC(const amrex::Real c[], amrex::Real& wtm)
{
  // See Eq 5 in CK Manual
  amrex::Real W = 0;
  amrex::Real sumC = 0;
  W += c[0] * 154.297000; // POSF10325
  W += c[1] * 31.998000;  // O2
  W += c[2] * 28.014000;  // N2
  W += c[3] * 44.009000;  // CO2
  W += c[4] * 18.015000;  // H2O

  for (int id = 0; id < 5; ++id) {
    sumC += c[id];
  }
  //  CK provides no guard against division by zero
  wtm = W / sumC;
}

// get Cp/R as a function of T
// for all species (Eq 19)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPOR(const amrex::Real T, amrex::Real cpor[])
{
  cp_R(cpor, T);
}

// get H/RT as a function of T
// for all species (Eq 20)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHORT(const amrex::Real T, amrex::Real hort[])
{
  speciesEnthalpy(hort, T);
}

// get S/R as a function of T
// for all species (Eq 21)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKSOR(const amrex::Real T, amrex::Real sor[])
{
  speciesEntropy(sor, T);
}

// convert y[species] (mass fracs) to x[species] (mole fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKYTX(const amrex::Real y[], amrex::Real x[])
{
  amrex::Real YOW = 0;

  for (int i = 0; i < 5; i++) {
    YOW += y[i] * imw(i);
  }

  amrex::Real YOWINV = 1.0 / YOW;

  for (int i = 0; i < 5; i++) {
    x[i] = y[i] * imw(i) * YOWINV;
  }
}

// convert y[species] (mass fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKYTCP(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real y[],
  amrex::Real c[])
{
  amrex::Real YOW = 0;
  amrex::Real PWORT;

  // Compute inverse of mean molecular wt first
  for (int i = 0; i < 5; i++) {
    c[i] = y[i] * imw(i);
  }
  for (int i = 0; i < 5; i++) {
    YOW += c[i];
  }

  // PW/RT (see Eq. 7)
  PWORT = P / (YOW * 8.31446261815324e+07 * T);
  // Now compute conversion

  for (int i = 0; i < 5; i++) {
    c[i] = PWORT * y[i] * imw(i);
  }
}

// convert y[species] (mass fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKYTCR(
  const amrex::Real rho,
  amrex::Real /*T*/,
  const amrex::Real y[],
  amrex::Real c[])
{

  for (int i = 0; i < 5; i++) {
    c[i] = rho * y[i] * imw(i);
  }
}

// convert x[species] (mole fracs) to y[species] (mass fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKXTY(const amrex::Real x[], amrex::Real y[])
{
  amrex::Real XW = 0; // See Eq 4, 9 in CK Manual
  // Compute mean molecular wt first
  XW += x[0] * 154.297000; // POSF10325
  XW += x[1] * 31.998000;  // O2
  XW += x[2] * 28.014000;  // N2
  XW += x[3] * 44.009000;  // CO2
  XW += x[4] * 18.015000;  // H2O
  // Now compute conversion
  amrex::Real XWinv = 1.0 / XW;
  y[0] = x[0] * 154.297000 * XWinv;
  y[1] = x[1] * 31.998000 * XWinv;
  y[2] = x[2] * 28.014000 * XWinv;
  y[3] = x[3] * 44.009000 * XWinv;
  y[4] = x[4] * 18.015000 * XWinv;
}

// convert x[species] (mole fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKXTCP(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real x[],
  amrex::Real c[])
{
  amrex::Real PORT = P / (8.31446261815324e+07 * T); // P/RT

  // Compute conversion, see Eq 10
  for (int id = 0; id < 5; ++id) {
    c[id] = x[id] * PORT;
  }
}

// convert x[species] (mole fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKXTCR(
  const amrex::Real rho,
  const amrex::Real /*T*/,
  const amrex::Real x[],
  amrex::Real c[])
{
  amrex::Real XW = 0; // See Eq 4, 11 in CK Manual
  amrex::Real ROW;
  // Compute mean molecular wt first
  XW += x[0] * 154.297000; // POSF10325
  XW += x[1] * 31.998000;  // O2
  XW += x[2] * 28.014000;  // N2
  XW += x[3] * 44.009000;  // CO2
  XW += x[4] * 18.015000;  // H2O
  ROW = rho / XW;

  // Compute conversion, see Eq 11
  for (int id = 0; id < 5; ++id) {
    c[id] = x[id] * ROW;
  }
}

// convert c[species] (molar conc) to x[species] (mole fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCTX(const amrex::Real c[], amrex::Real x[])
{
  amrex::Real sumC = 0;

  // compute sum of c
  for (int id = 0; id < 5; ++id) {
    sumC += c[id];
  }

  //  See Eq 13
  amrex::Real sumCinv = 1.0 / sumC;
  for (int id = 0; id < 5; ++id) {
    x[id] = c[id] * sumCinv;
  }
}

// convert c[species] (molar conc) to y[species] (mass fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCTY(const amrex::Real c[], amrex::Real y[])
{
  amrex::Real CW = 0; // See Eq 12 in CK Manual
  // compute denominator in eq 12 first
  CW += c[0] * 154.297000; // POSF10325
  CW += c[1] * 31.998000;  // O2
  CW += c[2] * 28.014000;  // N2
  CW += c[3] * 44.009000;  // CO2
  CW += c[4] * 18.015000;  // H2O
  // Now compute conversion
  amrex::Real CWinv = 1.0 / CW;
  y[0] = c[0] * 154.297000 * CWinv;
  y[1] = c[1] * 31.998000 * CWinv;
  y[2] = c[2] * 28.014000 * CWinv;
  y[3] = c[3] * 44.009000 * CWinv;
  y[4] = c[4] * 18.015000 * CWinv;
}

// get specific heat at constant volume as a function
// of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCVML(const amrex::Real T, amrex::Real cvml[])
{
  cv_R(cvml, T);

  // convert to chemkin units
  for (int id = 0; id < 5; ++id) {
    cvml[id] *= 8.31446261815324e+07;
  }
}

// get specific heat at constant pressure as a
// function of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPML(const amrex::Real T, amrex::Real cpml[])
{
  cp_R(cpml, T);

  // convert to chemkin units
  for (int id = 0; id < 5; ++id) {
    cpml[id] *= 8.31446261815324e+07;
  }
}

// get internal energy as a function
// of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUML(const amrex::Real T, amrex::Real uml[])
{
  amrex::Real RT = 8.31446261815324e+07 * T; // R*T
  speciesInternalEnergy(uml, T);

  // convert to chemkin units
  for (int id = 0; id < 5; ++id) {
    uml[id] *= RT;
  }
}

// get enthalpy as a function
// of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHML(const amrex::Real T, amrex::Real hml[])
{
  amrex::Real RT = 8.31446261815324e+07 * T; // R*T
  speciesEnthalpy(hml, T);

  // convert to chemkin units
  for (int id = 0; id < 5; ++id) {
    hml[id] *= RT;
  }
}

// Returns the standard-state entropies in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKSML(const amrex::Real T, amrex::Real sml[])
{
  speciesEntropy(sml, T);

  // convert to chemkin units
  for (int id = 0; id < 5; ++id) {
    sml[id] *= 8.31446261815324e+07;
  }
}

// Returns the specific heats at constant volume
// in mass units (Eq. 29)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCVMS(const amrex::Real T, amrex::Real cvms[])
{
  cv_R(cvms, T);
  // multiply by R/molecularweight
  cvms[0] *= 5.388609382005639e+05; // POSF10325
  cvms[1] *= 2.598431970171023e+06; // O2
  cvms[2] *= 2.967966951578939e+06; // N2
  cvms[3] *= 1.889264154639560e+06; // CO2
  cvms[4] *= 4.615299815794193e+06; // H2O
}

// Returns the specific heats at constant pressure
// in mass units (Eq. 26)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPMS(const amrex::Real T, amrex::Real cpms[])
{
  cp_R(cpms, T);
  // multiply by R/molecularweight
  cpms[0] *= 5.388609382005639e+05; // POSF10325
  cpms[1] *= 2.598431970171023e+06; // O2
  cpms[2] *= 2.967966951578939e+06; // N2
  cpms[3] *= 1.889264154639560e+06; // CO2
  cpms[4] *= 4.615299815794193e+06; // H2O
}

// Returns internal energy in mass units (Eq 30.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUMS(const amrex::Real T, amrex::Real ums[])
{
  amrex::Real RT = 8.31446261815324e+07 * T; // R*T

  speciesInternalEnergy(ums, T);

  for (int i = 0; i < 5; i++) {
    ums[i] *= RT * imw(i);
  }
}

// Returns enthalpy in mass units (Eq 27.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHMS(const amrex::Real T, amrex::Real hms[])
{
  amrex::Real RT = 8.31446261815324e+07 * T; // R*T

  speciesEnthalpy(hms, T);

  for (int i = 0; i < 5; i++) {
    hms[i] *= RT * imw(i);
  }
}

// Returns the entropies in mass units (Eq 28.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKSMS(const amrex::Real T, amrex::Real sms[])
{
  speciesEntropy(sms, T);
  // multiply by R/molecularweight
  sms[0] *= 5.388609382005639e+05; // POSF10325
  sms[1] *= 2.598431970171023e+06; // O2
  sms[2] *= 2.967966951578939e+06; // N2
  sms[3] *= 1.889264154639560e+06; // CO2
  sms[4] *= 4.615299815794193e+06; // H2O
}

// GPU version of productionRate: no more use of thermo namespace vectors
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
comp_qfqr(
  amrex::Real* /*qf*/,
  amrex::Real* /*qr*/,
  const amrex::Real* /*sc*/,
  const amrex::Real* /*sc_qss*/,
  const amrex::Real /*T*/,
  const amrex::Real /*invT*/,
  const amrex::Real /*logT*/)
{
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
productionRate(
  amrex::Real* wdot, const amrex::Real* /*sc*/, const amrex::Real /*T*/)
{

  for (int i = 0; i < 5; ++i) {
    wdot[i] = 0.0;
  }
}

// compute the production rate for each species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWC(const amrex::Real T, amrex::Real C[], amrex::Real wdot[])
{

  // convert to SI
  for (int id = 0; id < 5; ++id) {
    C[id] *= 1.0e6;
  }

  // convert to chemkin units
  productionRate(wdot, C, T);

  // convert to chemkin units
  for (int id = 0; id < 5; ++id) {
    C[id] *= 1.0e-6;
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given P, T, and mass fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWYP(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real y[],
  amrex::Real wdot[])
{
  amrex::Real c[5]; // temporary storage
  amrex::Real YOW = 0;
  amrex::Real PWORT;

  // Compute inverse of mean molecular wt first
  for (int i = 0; i < 5; i++) {
    YOW += y[i] * imw(i);
  }
  // PW/RT (see Eq. 7)
  PWORT = P / (YOW * 8.31446261815324e+07 * T);
  // multiply by 1e6 so c goes to SI
  PWORT *= 1e6;
  // Now compute conversion (and go to SI)
  for (int i = 0; i < 5; i++) {
    c[i] = PWORT * y[i] * imw(i);
  }

  // convert to chemkin units
  productionRate(wdot, c, T);

  // convert to chemkin units
  for (int id = 0; id < 5; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given P, T, and mole fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWXP(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real x[],
  amrex::Real wdot[])
{
  amrex::Real c[5]; // temporary storage
  amrex::Real PORT =
    1e6 * P / (8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (int id = 0; id < 5; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  productionRate(wdot, c, T);

  // convert to chemkin units
  for (int id = 0; id < 5; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given rho, T, and mass fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWYR(
  const amrex::Real rho,
  const amrex::Real T,
  const amrex::Real y[],
  amrex::Real wdot[])
{
  amrex::Real c[5]; // temporary storage

  // See Eq 8 with an extra 1e6 so c goes to SI
  for (int i = 0; i < 5; i++) {
    c[i] = 1e6 * rho * y[i] * imw(i);
  }

  // call productionRate
  productionRate(wdot, c, T);

  // convert to chemkin units
  for (int id = 0; id < 5; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given rho, T, and mole fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWXR(
  const amrex::Real rho,
  const amrex::Real T,
  const amrex::Real x[],
  amrex::Real wdot[])
{
  amrex::Real c[5];   // temporary storage
  amrex::Real XW = 0; // See Eq 4, 11 in CK Manual
  amrex::Real ROW;
  // Compute mean molecular wt first
  XW += x[0] * 154.297000; // POSF10325
  XW += x[1] * 31.998000;  // O2
  XW += x[2] * 28.014000;  // N2
  XW += x[3] * 44.009000;  // CO2
  XW += x[4] * 18.015000;  // H2O
  // Extra 1e6 factor to take c to SI
  ROW = 1e6 * rho / XW;

  // Compute conversion, see Eq 11
  for (int id = 0; id < 5; ++id) {
    c[id] = x[id] * ROW;
  }

  // convert to chemkin units
  productionRate(wdot, c, T);

  // convert to chemkin units
  for (int id = 0; id < 5; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

//  species unit charge number
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCHRG(int kcharge[])
{
  kcharge[0] = 0; // POSF10325
  kcharge[1] = 0; // O2
  kcharge[2] = 0; // N2
  kcharge[3] = 0; // CO2
  kcharge[4] = 0; // H2O
}

//  species charge per unit mass
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCHRGMASS(amrex::Real zk[])
{

  int kchrg[5];
  CKCHRG(kchrg);

  for (int id = 0; id < 5; ++id) {
    zk[id] = 6.02214076e+23 * 1.60217663e-19 * kchrg[id] * imw(id);
  }
}

// compute an approx to the reaction Jacobian
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
aJacobian_precond(
  amrex::Real* J,
  const amrex::Real* /*sc*/,
  const amrex::Real /*T*/,
  const int /*HP*/)
{

#if defined(PELE_COMPILE_AJACOBIAN) || !defined(AMREX_USE_HIP)
  for (int i = 0; i < 36; i++) {
    J[i] = 0.0;
  }
#else
  amrex::Abort();
#endif
}

// compute an approx to the reaction Jacobian (for preconditioning)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
DWDOT_SIMPLIFIED(
  amrex::Real* J, const amrex::Real* sc, const amrex::Real* Tp, const int* HP)
{
  amrex::Real c[5];

  for (int k = 0; k < 5; k++) {
    c[k] = 1.e6 * sc[k];
  }

  aJacobian_precond(J, c, *Tp, *HP);

  // dwdot[k]/dT
  // dTdot/d[X]
  for (int k = 0; k < 5; k++) {
    J[30 + k] *= 1.e-6;
    J[k * 6 + 5] *= 1.e6;
  }
}

// compute the reaction Jacobian
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
aJacobian(
  amrex::Real* J,
  const amrex::Real* /*sc*/,
  const amrex::Real /*T*/,
  const int /*consP*/)
{

#if defined(PELE_COMPILE_AJACOBIAN) || !defined(AMREX_USE_HIP)
  for (int i = 0; i < 36; i++) {
    J[i] = 0.0;
  }
#else
  amrex::Abort();
#endif
}

// compute the reaction Jacobian
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
DWDOT(
  amrex::Real* J,
  const amrex::Real* sc,
  const amrex::Real* Tp,
  const int* consP)
{
  amrex::Real c[5];

  for (int k = 0; k < 5; k++) {
    c[k] = 1.e6 * sc[k];
  }

  aJacobian(J, c, *Tp, *consP);

  // dwdot[k]/dT
  // dTdot/d[X]
  for (int k = 0; k < 5; k++) {
    J[30 + k] *= 1.e-6;
    J[k * 6 + 5] *= 1.e6;
  }
}

// Transport function declarations

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetLENIMC(int* LENIMC)
{
  *LENIMC = 20;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetLENRMC(int* LENRMC)
{
  *LENRMC = 610;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetNO(int* NO)
{
  *NO = 4;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetKK(int* KK)
{
  *KK = 5;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetNLITE(int* NLITE)
{
  *NLITE = 0;
}

// Patm in ergs/cm3
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetPATM(amrex::Real* PATM)
{
  *PATM = 0.1013250000000000E+07;
}

// the molecular weights in g/mol
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetWT(amrex::Real* WT)
{
  WT[0] = 1.54297000E+02;
  WT[1] = 3.19980000E+01;
  WT[2] = 2.80140000E+01;
  WT[3] = 4.40090000E+01;
  WT[4] = 1.80150000E+01;
}

// the lennard-jones potential well depth eps/kb in K
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetEPS(amrex::Real* EPS)
{
  EPS[0] = 1.77645000E+02;
  EPS[1] = 1.07400000E+02;
  EPS[2] = 9.75300000E+01;
  EPS[3] = 2.44000000E+02;
  EPS[4] = 5.72400000E+02;
}

// the lennard-jones collision diameter in Angstroms
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetSIG(amrex::Real* SIG)
{
  SIG[0] = 9.60500000E+00;
  SIG[1] = 3.45800000E+00;
  SIG[2] = 3.62100000E+00;
  SIG[3] = 3.76300000E+00;
  SIG[4] = 2.60500000E+00;
}

// the dipole moment in Debye
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetDIP(amrex::Real* DIP)
{
  DIP[0] = 0.00000000E+00;
  DIP[1] = 0.00000000E+00;
  DIP[2] = 0.00000000E+00;
  DIP[3] = 0.00000000E+00;
  DIP[4] = 1.84400000E+00;
}

// the polarizability in cubic Angstroms
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetPOL(amrex::Real* POL)
{
  POL[0] = 0.00000000E+00;
  POL[1] = 1.60000000E+00;
  POL[2] = 1.76000000E+00;
  POL[3] = 2.65000000E+00;
  POL[4] = 0.00000000E+00;
}

// the rotational relaxation collision number at 298 K
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetZROT(amrex::Real* ZROT)
{
  ZROT[0] = 1.00000000E+00;
  ZROT[1] = 3.80000000E+00;
  ZROT[2] = 4.00000000E+00;
  ZROT[3] = 2.10000000E+00;
  ZROT[4] = 4.00000000E+00;
}

// 0: monoatomic, 1: linear, 2: nonlinear
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetNLIN(int* NLIN)
{
  NLIN[0] = 2;
  NLIN[1] = 1;
  NLIN[2] = 1;
  NLIN[3] = 1;
  NLIN[4] = 2;
}

// Poly fits for the viscosities, dim NO*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFETA(amrex::Real* COFETA)
{
  COFETA[0] = -2.27773842E+01;
  COFETA[1] = 4.31076442E+00;
  COFETA[2] = -4.79725560E-01;
  COFETA[3] = 2.09551504E-02;
  COFETA[4] = -1.71618439E+01;
  COFETA[5] = 2.68036374E+00;
  COFETA[6] = -2.72570227E-01;
  COFETA[7] = 1.21650964E-02;
  COFETA[8] = -1.65695492E+01;
  COFETA[9] = 2.39056562E+00;
  COFETA[10] = -2.34558144E-01;
  COFETA[11] = 1.05024037E-02;
  COFETA[12] = -2.40015087E+01;
  COFETA[13] = 5.14359547E+00;
  COFETA[14] = -5.74269731E-01;
  COFETA[15] = 2.44937679E-02;
  COFETA[16] = -1.05420930E+01;
  COFETA[17] = -1.37777266E+00;
  COFETA[18] = 4.20502561E-01;
  COFETA[19] = -2.40627344E-02;
}

// Poly fits for the conductivities, dim NO*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFLAM(amrex::Real* COFLAM)
{
  COFLAM[0] = -2.28849559E+01;
  COFLAM[1] = 9.72792992E+00;
  COFLAM[2] = -1.00067773E+00;
  COFLAM[3] = 3.67484895E-02;
  COFLAM[4] = -1.93717634E+00;
  COFLAM[5] = 2.89110257E+00;
  COFLAM[6] = -2.71096988E-01;
  COFLAM[7] = 1.15344939E-02;
  COFLAM[8] = 1.29305954E+01;
  COFLAM[9] = -3.52816965E+00;
  COFLAM[10] = 6.45498434E-01;
  COFLAM[11] = -3.19375023E-02;
  COFLAM[12] = -1.13649201E+01;
  COFLAM[13] = 5.88177354E+00;
  COFLAM[14] = -5.68651792E-01;
  COFLAM[15] = 2.03561484E-02;
  COFLAM[16] = 2.33729836E+01;
  COFLAM[17] = -8.96536178E+00;
  COFLAM[18] = 1.52828630E+00;
  COFLAM[19] = -7.58551616E-02;
}

// Poly fits for the diffusion coefficients, dim NO*KK*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFD(amrex::Real* COFD)
{
  COFD[0] = -2.18108044E+01;
  COFD[1] = 4.73185706E+00;
  COFD[2] = -3.97136602E-01;
  COFD[3] = 1.71287819E-02;
  COFD[4] = -1.87890298E+01;
  COFD[5] = 4.11991398E+00;
  COFD[6] = -3.22514659E-01;
  COFD[7] = 1.40877517E-02;
  COFD[8] = -1.84342672E+01;
  COFD[9] = 4.00051993E+00;
  COFD[10] = -3.07414279E-01;
  COFD[11] = 1.34475779E-02;
  COFD[12] = -2.18171057E+01;
  COFD[13] = 5.12625274E+00;
  COFD[14] = -4.43796461E-01;
  COFD[15] = 1.89713627E-02;
  COFD[16] = -2.30218519E+01;
  COFD[17] = 5.57736917E+00;
  COFD[18] = -4.77603538E-01;
  COFD[19] = 1.94525527E-02;
  COFD[20] = -1.87890298E+01;
  COFD[21] = 4.11991398E+00;
  COFD[22] = -3.22514659E-01;
  COFD[23] = 1.40877517E-02;
  COFD[24] = -1.55511233E+01;
  COFD[25] = 3.48070094E+00;
  COFD[26] = -2.40859499E-01;
  COFD[27] = 1.05972514E-02;
  COFD[28] = -1.52414498E+01;
  COFD[29] = 3.35922578E+00;
  COFD[30] = -2.25181399E-01;
  COFD[31] = 9.92132878E-03;
  COFD[32] = -1.84688303E+01;
  COFD[33] = 4.49330851E+00;
  COFD[34] = -3.68208715E-01;
  COFD[35] = 1.59565402E-02;
  COFD[36] = -2.12652444E+01;
  COFD[37] = 5.59961822E+00;
  COFD[38] = -4.91624856E-01;
  COFD[39] = 2.05035546E-02;
  COFD[40] = -1.84342672E+01;
  COFD[41] = 4.00051993E+00;
  COFD[42] = -3.07414279E-01;
  COFD[43] = 1.34475779E-02;
  COFD[44] = -1.52414498E+01;
  COFD[45] = 3.35922578E+00;
  COFD[46] = -2.25181399E-01;
  COFD[47] = 9.92132878E-03;
  COFD[48] = -1.49828551E+01;
  COFD[49] = 3.25781069E+00;
  COFD[50] = -2.12199367E-01;
  COFD[51] = 9.36657283E-03;
  COFD[52] = -1.81432499E+01;
  COFD[53] = 4.37565431E+00;
  COFD[54] = -3.53906025E-01;
  COFD[55] = 1.53760786E-02;
  COFD[56] = -2.10643264E+01;
  COFD[57] = 5.53614864E+00;
  COFD[58] = -4.86046753E-01;
  COFD[59] = 2.03659193E-02;
  COFD[60] = -2.18171057E+01;
  COFD[61] = 5.12625274E+00;
  COFD[62] = -4.43796461E-01;
  COFD[63] = 1.89713627E-02;
  COFD[64] = -1.84688303E+01;
  COFD[65] = 4.49330851E+00;
  COFD[66] = -3.68208715E-01;
  COFD[67] = 1.59565402E-02;
  COFD[68] = -1.81432499E+01;
  COFD[69] = 4.37565431E+00;
  COFD[70] = -3.53906025E-01;
  COFD[71] = 1.53760786E-02;
  COFD[72] = -2.13425604E+01;
  COFD[73] = 5.40460130E+00;
  COFD[74] = -4.72718910E-01;
  COFD[75] = 1.99362717E-02;
  COFD[76] = -2.07653629E+01;
  COFD[77] = 5.01091988E+00;
  COFD[78] = -3.77985581E-01;
  COFD[79] = 1.40968617E-02;
  COFD[80] = -2.30218519E+01;
  COFD[81] = 5.57736917E+00;
  COFD[82] = -4.77603538E-01;
  COFD[83] = 1.94525527E-02;
  COFD[84] = -2.12652444E+01;
  COFD[85] = 5.59961822E+00;
  COFD[86] = -4.91624856E-01;
  COFD[87] = 2.05035546E-02;
  COFD[88] = -2.10643264E+01;
  COFD[89] = 5.53614864E+00;
  COFD[90] = -4.86046753E-01;
  COFD[91] = 2.03659193E-02;
  COFD[92] = -2.07653629E+01;
  COFD[93] = 5.01091988E+00;
  COFD[94] = -3.77985581E-01;
  COFD[95] = 1.40968617E-02;
  COFD[96] = -1.19157788E+01;
  COFD[97] = 9.28952533E-01;
  COFD[98] = 2.42107483E-01;
  COFD[99] = -1.59824150E-02;
}

// List of specs with small weight, dim NLITE
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetKTDIF(int* /*KTDIF*/)
{
}

// Poly fits for thermal diff ratios, dim NO*NLITE*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFTD(amrex::Real* /*COFTD*/)
{
}

// compute the critical parameters for each species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
GET_CRITPARAMS(
  amrex::Real* Tci, amrex::Real* ai, amrex::Real* bi, amrex::Real* acentric_i)
{

  amrex::Real EPS[5];
  amrex::Real SIG[5];
  amrex::Real wt[5];
  amrex::Real Rcst = 83.144598; // in bar [CGS] !
  amrex::Real avogadro = 6.02214199e23;
  amrex::Real boltzmann = 1.3806503e-16; // we work in CGS

  egtransetEPS(EPS);
  egtransetSIG(SIG);
  get_mw(wt);

  // species 0: POSF10325
  Tci[0] = 1.316 * EPS[0];
  ai[0] = (5.55 * avogadro * avogadro * EPS[0] * boltzmann * 1e-24 * SIG[0] *
           SIG[0] * SIG[0]) /
          (wt[0] * wt[0]);
  bi[0] = 0.855 * avogadro * 1e-24 * SIG[0] * SIG[0] * SIG[0] / (wt[0]);
  acentric_i[0] = 0.0;

  // species 1: O2
  // Imported from NIST
  Tci[1] = 154.581000;
  ai[1] = 1e6 * 0.42748 * Rcst * Rcst * Tci[1] * Tci[1] /
          (31.998800 * 31.998800 * 50.430466);
  bi[1] = 0.08664 * Rcst * Tci[1] / (31.998800 * 50.430466);
  acentric_i[1] = 0.022200;

  // species 2: N2
  // Imported from NIST
  Tci[2] = 126.192000;
  ai[2] = 1e6 * 0.42748 * Rcst * Rcst * Tci[2] * Tci[2] /
          (28.013400 * 28.013400 * 33.958000);
  bi[2] = 0.08664 * Rcst * Tci[2] / (28.013400 * 33.958000);
  acentric_i[2] = 0.037200;

  // species 3: CO2
  // Imported from NIST
  Tci[3] = 304.120000;
  ai[3] = 1e6 * 0.42748 * Rcst * Rcst * Tci[3] * Tci[3] /
          (44.009950 * 44.009950 * 73.740000);
  bi[3] = 0.08664 * Rcst * Tci[3] / (44.009950 * 73.740000);
  acentric_i[3] = 0.225000;

  // species 4: H2O
  // Imported from NIST
  Tci[4] = 647.096000;
  ai[4] = 1e6 * 0.42748 * Rcst * Rcst * Tci[4] * Tci[4] /
          (18.015340 * 18.015340 * 220.640000);
  bi[4] = 0.08664 * Rcst * Tci[4] / (18.015340 * 220.640000);
  acentric_i[4] = 0.344300;
}

// compute the critical parameter quantities for each species for SRK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
GET_CRITPARAMS_SRK(
  amrex::Real* sqrtOneOverTc,
  amrex::Real* sqrtAsti,
  amrex::Real* Bi,
  amrex::Real* Fomega)
{

  // species 0: POSF10325
  // Computed from Lennard-Jones
  sqrtOneOverTc[0] = 6.5402682478301e-02;
  sqrtAsti[0] = 4.2865154585193e+04;
  Bi[0] = 2.9570027693584e+00;
  Fomega[0] = 4.8508000000000e-01;

  // species 1: O2
  // Imported from NIST
  sqrtOneOverTc[1] = 8.0430717653241e-02;
  sqrtAsti[1] = 3.6980081924740e+04;
  Bi[1] = 6.9005216444999e-01;
  Fomega[1] = 5.1945301904908e-01;

  // species 2: N2
  // Imported from NIST
  sqrtOneOverTc[2] = 8.9019282240563e-02;
  sqrtAsti[2] = 4.2022983929562e+04;
  Bi[2] = 9.5560052221563e-01;
  Fomega[2] = 5.4259343186608e-01;

  // species 3: CO2
  // Imported from NIST
  sqrtOneOverTc[3] = 5.7342616962522e-02;
  sqrtAsti[3] = 4.3745610142919e+04;
  Bi[3] = 6.7506167841745e-01;
  Fomega[3] = 8.2653709187500e-01;

  // species 4: H2O
  // Imported from NIST
  sqrtOneOverTc[4] = 3.9311140369715e-02;
  sqrtAsti[4] = 1.3145468419592e+05;
  Bi[4] = 1.1727204136223e+00;
  Fomega[4] = 1.0013577274636e+00;
}

// gauss-jordan solver external routine
// Replace this routine with the one generated by the Gauss Jordan solver of DW
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
sgjsolve(amrex::Real* /*A*/, amrex::Real* /*x*/, amrex::Real* /*b*/)
{
  amrex::Abort("sgjsolve not implemented, choose a different solver ");
}

// Replace this routine with the one generated by the Gauss Jordan solver of DW
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
sgjsolve_simplified(amrex::Real* /*A*/, amrex::Real* /*x*/, amrex::Real* /*b*/)
{
  amrex::Abort(
    "sgjsolve_simplified not implemented, choose a different solver ");
}

#endif
