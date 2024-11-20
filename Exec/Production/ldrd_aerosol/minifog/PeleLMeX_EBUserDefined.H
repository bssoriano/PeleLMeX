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

    amrex::Print() << "** Initialising MiniFog geometry... \n";
    
    amrex::ParmParse pp("EB");
    amrex::Real x_lower_point = 0.0;
    pp.query("x_lower_point",x_lower_point);
    amrex::Real y_lower_point = 0.0;
    pp.query("y_lower_point",y_lower_point);
    amrex::Real x_higher_point = 0.0;
    pp.query("x_higher_point",x_higher_point);
    amrex::Real y_higher_point = 0.0;
    pp.query("y_higher_point",y_higher_point);
    
    amrex::Array<amrex::Real,AMREX_SPACEDIM> point0_cat;
    amrex::Array<amrex::Real,AMREX_SPACEDIM> point1_cat;

    point1_cat[0] =  x_lower_point; //x-coordinate
    point1_cat[1] =  y_lower_point; //y-coordinate

    point0_cat[0] =  x_higher_point; //x-coordinate
    point0_cat[1] =  y_higher_point; //y-coordinate

    amrex::Array<amrex::Real,AMREX_SPACEDIM> norm0;

    norm0[0] =  (point0_cat[0]-point1_cat[0]);
    norm0[1] =  (point0_cat[1]-point1_cat[1]); 
    norm0[2] = 0.0;

    amrex::Real norm = sqrt(norm0[0]*norm0[0]+norm0[1]*norm0[1]);
    norm0[0] = norm0[0]/norm;
    norm0[1] = norm0[1]/norm;
    
    amrex::Real loc_x;
    pp.query("loc_x", loc_x );
    amrex::Real loc_y;
    pp.query("loc_y", loc_y );    
    EB2::PlaneIF wall1({AMREX_D_DECL(loc_x,0.0,loc_y)},{AMREX_D_DECL(-norm0[0],0,norm0[1])});


    point1_cat[0] = -x_lower_point; //x-coordinate
    point1_cat[1] =  y_lower_point; //y-coordinate

    point0_cat[0] = -x_higher_point; //x-coordinate
    point0_cat[1] =  y_higher_point; //y-coordinate
    
    norm0[0] =  (point0_cat[0]-point1_cat[0]);
    norm0[1] =  (point0_cat[1]-point1_cat[1]); 
    norm0[2] = 0.0;

    norm = sqrt(norm0[0]*norm0[0]+norm0[1]*norm0[1]);
    norm0[0] = norm0[0]/norm;
    norm0[1] = norm0[1]/norm;
    EB2::PlaneIF wall2({AMREX_D_DECL(-loc_x,0.0,loc_y)},{AMREX_D_DECL(-norm0[0],0,norm0[1])});


    EB2::PlaneIF top_wall({AMREX_D_DECL(0.0,0.0,y_higher_point)},{AMREX_D_DECL(0,0,1)});

    bool has_fluid_inside=true;
    int direction=2;
    RealArray center;
    pp.get("cylinder_center", center);
    Real radius;
    pp.get("cylinder_radius", radius);
    Real height = -1.0;
    pp.queryAdd("cylinder_height", height);
    
    EB2::CylinderIF outlet(radius, height, direction, center, has_fluid_inside);

    //EB2::PlaneIF outlet({AMREX_D_DECL(0.0,radius,y_higher_point)},{AMREX_D_DECL(1,0,0)});

    auto top_intersection = EB2::makeIntersection(outlet,top_wall);
    //auto top_intersection = EB2::makeIntersection(top_wall);

    //auto polys = EB2::makeUnion(wall1,wall2);
    auto polys = EB2::makeUnion(wall1,wall2,top_intersection);
    //auto polys = EB2::makeUnion(outlet);

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
