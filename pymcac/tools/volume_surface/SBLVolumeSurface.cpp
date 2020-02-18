/*
 * MCAC
 * Copyright (C) 2020 CORIA
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Fixed_alpha_shape_vertex_base_3.h>
#include <CGAL/Fixed_alpha_shape_cell_base_3.h>
#include <CGAL/Fixed_alpha_shape_3.h>
#include <CGAL/Exact_spherical_kernel_3.h>

#include <SBLVolumeSurface.hpp>
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

//#ifndef NDEBUG
typedef CGAL::Exact_spherical_kernel_3                           SK;
//#else // This does not work, crash execution if try to use it
//typedef CGAL::Simple_cartesian<double>                           KE;
//typedef CGAL::Algebraic_kernel_for_spheres_2_3<double>           AK;
//typedef CGAL::Spherical_kernel_3<KE, AK>                         SK;
//#endif

typedef SBL::GT::T_Union_of_balls_surface_volume_3<Alpha_complex, SK>    Union_of_balls_surface_volume_3;


void compute_volume_surface(const double spherelist[],
                            const int& nspheres,
                            double& volume,
                            double& area)
{
  // copy data into sbl format
  std::vector<K::Sphere_3> spheres;
  for (int i=0; i<nspheres*4;i+=4)
    {
      K::Point_3 p(spherelist[i],spherelist[i+1],spherelist[i+2]);
      K::FT r(spherelist[i+3]);
      K::Sphere_3 sphere(p, r*r);
      spheres.push_back(sphere);
    }

  Union_of_balls_surface_volume_3 surface_volume(spheres.begin(), spheres.end());

  CGAL::Interval_nt<false> cgalvolume = surface_volume.volume();
  CGAL::Interval_nt<false> cgalarea = surface_volume.area();

  volume = to_double(cgalvolume);
  area = to_double(cgalarea);
}
