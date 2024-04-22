#include <PeleLMeX.H>
#include <AMReX_TagBox.H>
#ifdef AMREX_USE_EB
#include <AMReX_EBAmrUtil.H>
#endif

using namespace amrex;

void
PeleLM::ErrorEst(int lev, TagBoxArray& tags, Real time, int /*ng*/)
{
  BL_PROFILE("PeleLMeX::ErrorEst()");

#ifdef AMREX_USE_EB
  // Tag EB up to m_EB_refine_LevMax-1 if Static or
  //              m_EB_refine_LevAdapt-1 if Adaptive
  if (
    (m_EB_refine_type == "Static" && lev < m_EB_refine_LevMax) ||
    (m_EB_refine_type == "Adaptive" && lev < m_EB_refine_LevAdapt)) {
    const MultiFab& state = (getLevelDataPtr(lev, AmrNewTime))->state;
    TagCutCells(tags, state);
  }
#endif

  for (const auto& errTag : errTags) {
    std::unique_ptr<MultiFab> mf;
    if (!errTag.Field().empty()) {
      mf = deriveComp(errTag.Field(), time, lev, errTag.NGrow());
    }
    errTag(tags, mf.get(), TagBox::CLEAR, TagBox::SET, time, lev, geom[lev]);
  }

#ifdef AMREX_USE_EB
  // Untag covered cells
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(tags, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const auto& bx = mfi.tilebox();
    auto tag = tags.array(mfi);
    auto vfrac = EBFactory(lev).getVolFrac().const_array(mfi);
    amrex::ParallelFor(bx, [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) {
      if (vfrac(i, j, k) <= 0.0) {
        tag(i, j, k) = TagBox::CLEAR;
      }
    });
  }

  // Untag cell close to EB
  if (m_EB_refine_type == "Static" && lev >= m_EB_refine_LevMax) {
    // Get distance function at current level
    MultiFab signDist(grids[lev], dmap[lev], 1, 0, MFInfo(), EBFactory(lev));
    getEBDistance(lev, signDist);

    // Estimate how far I need to derefine
    Real diagFac = std::sqrt(2.0) * m_derefineEBBuffer;
    Real clearTagDist = Geom(m_EB_refine_LevMax).CellSize(0) *
                        static_cast<Real>(nErrorBuf(m_EB_refine_LevMax)) *
                        diagFac;
    for (int ilev = m_EB_refine_LevMax + 1; ilev <= finest_level; ++ilev) {
      clearTagDist += static_cast<Real>(nErrorBuf(ilev)) *
                      Geom(m_EB_refine_LevMax).CellSize(0) * diagFac;
    }

    // Untag cells too close to EB
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(tags, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const auto& bx = mfi.tilebox();
      const auto& dist = signDist.const_array(mfi);
      auto tag = tags.array(mfi);
      amrex::ParallelFor(bx, [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) {
        if (dist(i, j, k) < clearTagDist) {
          tag(i, j, k) = TagBox::CLEAR;
        }
      });
    }
  }
#endif

   if (lev >= m_Zdir_refine_LevMax ) {
      // Untag cells when z is larger than m_MaxZdir_refine

      const amrex::Real* prob_lo = geom[lev].ProbLo();
      const amrex::Real* dx = geom[lev].CellSize();
      const amrex::Real prob_lo0 = prob_lo[0];
      const amrex::Real dx0 = dx[0];
      const amrex::Real prob_lo1 = prob_lo[1];
      const amrex::Real dx1 = dx[1];
      const amrex::Real prob_lo2 = prob_lo[2];
      const amrex::Real dx2 = dx[2];

      const amrex::Real MaxZdir_refine = m_MaxZdir_refine;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(tags,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
          const auto& bx    = mfi.tilebox();
          auto tag          = tags.array(mfi);
          // const auto& dist  = signDist.const_array(mfi);

          amrex::ParallelFor(bx,
          [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
          {
              const amrex::Real x = prob_lo0 + (i+0.5)*dx0;
              const amrex::Real y = prob_lo1 + (j+0.5)*dx1;
              const amrex::Real z = prob_lo2 + (k+0.5)*dx2;
              const amrex::Real r = sqrt(pow((0.0-x),2)+pow((0.0-y),2));

              if (z > MaxZdir_refine and r < 20.e-3) {
                  tag(i,j,k) = TagBox::CLEAR;
              }
          });
      }
   }


}
