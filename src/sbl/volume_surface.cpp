#ifdef WITH_SBL



#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Fixed_alpha_shape_vertex_base_3.h>
#include <CGAL/Fixed_alpha_shape_cell_base_3.h>
#include <CGAL/Fixed_alpha_shape_3.h>
#include <CGAL/Exact_spherical_kernel_3.h>
#include <SBL/GT/Union_of_balls_surface_volume_3.hpp>
#include <iostream>
#include <fstream>
typedef CGAL::Exact_predicates_inexact_constructions_kernel      K;
typedef CGAL::Regular_triangulation_euclidean_traits_3<K>        Traits;
typedef CGAL::Fixed_alpha_shape_vertex_base_3<Traits>            Vb;
typedef CGAL::Fixed_alpha_shape_cell_base_3<Traits>              Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb>             Tds;
typedef CGAL::Regular_triangulation_3<Traits, Tds>               Rt;
typedef CGAL::Fixed_alpha_shape_3<Rt>                            Alpha_complex;
#ifndef NDEBUG
typedef CGAL::Exact_spherical_kernel_3                           SK;
#else
typedef CGAL::Simple_cartesian<double>                           KE;
typedef CGAL::Algebraic_kernel_for_spheres_2_3<double>           AK;
typedef CGAL::Spherical_kernel_3<KE, AK>                         SK;
#endif
typedef SBL::GT::T_Union_of_balls_surface_volume_3<Alpha_complex, SK>    Union_of_balls_surface_volume_3;

#include "sblvolumesurface.hpp"

namespace MCAC{

pair<std::vector<double>, std::vector<double>> compute_volume_surface(const ListSphere& spherelist)
{
  // convert data for SBL
  std::vector<K::Sphere_3> spheres;
  for (const Sphere* sphere :spherelist)
    {
      std::array<double, 3> pos(sphere->Position());
      K::Point_3 p(pos[0],pos[1],pos[2]);
      // TODO temp workaround
      K::FT r(sphere->Radius() * (1 + 1e-10) );
      spheres.push_back(K::Sphere_3(p, r*r));
    }
  Union_of_balls_surface_volume_3 surface_volume(spheres.begin(), spheres.end());

  std::vector<double> volumes(spherelist.size());
  std::vector<double> areas(spherelist.size());

  for(unsigned int i=0; i<spherelist.size() ; ++i)
  {
      volumes[i] = to_double(surface_volume.volume(i));
      areas[i] = to_double(surface_volume.area(i));
  }

  return make_pair(volumes, areas);
}

}

#endif
