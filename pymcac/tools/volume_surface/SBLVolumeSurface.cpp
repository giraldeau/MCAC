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
// CGAL VERSION less than 4.10 are note supported
#define CGAL_VERSION_GEQ_4_10

#include <SBLVolumeSurface.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#ifndef CGAL_VERSION_GEQ_4_10
#ifndef _CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_3_H_
#define _CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_3_H_
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#endif
#else
#ifndef _CGAL_REGULAR_TRIANGULATION_VERTEX_BASE_3_H_
#define _CGAL_REGULAR_TRIANGULATION_VERTEX_BASE_3_H_
#include <CGAL/Regular_triangulation_vertex_base_3.h>
#endif
#endif

#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Fixed_alpha_shape_vertex_base_3.h>
#include <CGAL/Fixed_alpha_shape_cell_base_3.h>
#include <CGAL/Fixed_alpha_shape_3.h>

#include <CGAL/Exact_spherical_kernel_3.h>
#include <SBL/GT/Union_of_balls_surface_volume_3.hpp>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel      K;

// https://doc.cgal.org/latest/Alpha_shapes_3/Alpha_shapes_3_2ex_fixed_weighted_alpha_shapes_3_8cpp-example.html#_a6
#ifdef CGAL_VERSION_GEQ_4_10
typedef CGAL::Regular_triangulation_vertex_base_3<K>      Vbb;
typedef CGAL::Fixed_alpha_shape_vertex_base_3<K,Vbb>      Vb;
typedef CGAL::Regular_triangulation_cell_base_3<K>        Rcb;
typedef CGAL::Fixed_alpha_shape_cell_base_3<K,Rcb>        Cb;
typedef CGAL::Triangulation_data_structure_3<Vb,Cb>       Tds;
typedef CGAL::Regular_triangulation_3<K,Tds>              Triangulation_3;
typedef CGAL::Fixed_alpha_shape_3<Triangulation_3>        Alpha_complex;
#else
typedef K                                                        Traits;
typedef CGAL::Fixed_alpha_shape_vertex_base_3<Traits>            Vb;
typedef CGAL::Fixed_alpha_shape_cell_base_3<Traits>              Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb>             Tds;
typedef CGAL::Regular_triangulation_3<Traits, Tds>               Rt;
typedef CGAL::Fixed_alpha_shape_3<Rt>                            Alpha_complex;
#endif

//typedef CGAL::Interval_nt<false>                             FT;
typedef double                                               FT;
typedef CGAL::Simple_cartesian<FT>                                  KE;
typedef CGAL::Algebraic_kernel_for_spheres_2_3<FT>           AK;
typedef CGAL::Spherical_kernel_3<KE, AK>                     SK;

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

  volume = to_double(surface_volume.volume());
  area = to_double(surface_volume.area());
}
