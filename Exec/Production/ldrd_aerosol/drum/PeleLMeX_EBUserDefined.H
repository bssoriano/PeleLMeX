#ifndef _EBUSERDEFINED_H_
#define _EBUSERDEFINED_H_

using namespace amrex;

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EB2_IF_Lathe.H>
void
EBUserDefined(const Geometry& geom,
              const int required_coarsening_level,
              const int max_coarsening_level)
{
    // ParmParse your geometry parameters
    // --> fixed geometry of the challenge problem

    const amrex::Real *problo = geom.ProbLo();
    const amrex::Real *probhi = geom.ProbHi();
    amrex::Real dy = geom.CellSize()[1]; //* pow(2.0,max_coarsening_level);
    amrex::Real dx = geom.CellSize()[0]; //* pow(2.0,max_coarsening_level);

    amrex::Print() << "** Initializing drum geometry... \n";
    
    amrex::ParmParse pp("EB");
    
    bool has_fluid_inside=true;
    int direction=2;
    RealArray center;
    pp.get("cylinder_center", center);
    Real radius;
    pp.get("cylinder_radius", radius);
    Real height = -1.0;
    pp.queryAdd("cylinder_height", height);

    EB2::PlaneIF cyl_wall({AMREX_D_DECL(radius, 0., 0.)},
                          {AMREX_D_DECL(1.,0.,0.)});
    //EB2::CylinderIF outer_cyl(radius, height, direction, center, has_fluid_inside);
    
    has_fluid_inside=false;
    pp.get("small_cylinder_center", center);
    pp.get("small_cylinder_radius", radius);
    pp.queryAdd("small_cylinder_height", height);
    EB2::CylinderIF small_cyl(radius, height, direction, center, has_fluid_inside);

    //auto polys = EB2::makeUnion(outer_cyl);

    auto pr = EB2::lathe(cyl_wall);
    auto polys = EB2::makeUnion(pr,small_cyl);

    auto shop   = EB2::makeShop(polys);

    // EB2::Build(shop, geom, max_coarsening_level, max_coarsening_level);
    EB2::Build(shop, geom, required_coarsening_level,max_coarsening_level);

}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
setEBState(const amrex::Real xEBface[AMREX_SPACEDIM],
           amrex::Real s_ext[NVAR],
           const amrex::Real time,
           amrex::GeometryData const& /*geomdata*/,
           ProbParm const& /*prob_parm*/)
{
  // Provide the state on the EB centroid provided
  // by xEBface
  // Not the entire has to be specified. For instance, if using isothermal
  // boundary, only s_ext[TEMP] will be used.
  // Aborting here (message will not show on GPUs)
  //Abort("Default setEBState called. Copy and edit this file for your needs");
   
   amrex::Real r = sqrt(pow((0.0-xEBface[0]),2)+pow((0.0-xEBface[1]),2));
   if(r <= 12.5e-3 and xEBface[2] > -1.e-4){
     if(time < 9.e-3){
       s_ext[TEMP] = 2000.;
     }
     else{
       s_ext[TEMP] = 500.;
     }
   }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
setEBType(const amrex::Real * /*xEBface[AMREX_SPACEDIM]*/,
          amrex::Real & EBflagType,
          amrex::GeometryData const& /*geomdata*/,
          ProbParm const& /*prob_parm*/)
{
  // Provide a factor applied on the EB face diffusion
  // coefficient: 0.0 will turn off the diffusion (homogeneous Neumann BC)
  // while 1.0 will let it unchanged. Intermediate values can be specified
  // for partially Dirichlet BCs
  //Abort("Default setEBType called. Copy and edit this file for your needs");
   EBflagType = 1.0;
}
#endif
#endif
