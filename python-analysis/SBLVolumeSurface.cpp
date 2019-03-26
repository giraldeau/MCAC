#include <SBLVolumeSurface.hpp>
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
//#ifndef NDEBUG
typedef CGAL::Exact_spherical_kernel_3                           SK;
//#else // This does not work, crash compilation if try to use it
//typedef CGAL::Simple_cartesian<double>                           KE;
//typedef CGAL::Algebraic_kernel_for_spheres_2_3<double>           AK;
//typedef CGAL::Spherical_kernel_3<KE, AK>                         SK;
//#endif
typedef SBL::GT::T_Union_of_balls_surface_volume_3<Alpha_complex, SK>    Union_of_balls_surface_volume_3;


void compute_volume_surface(const double spherelist[], const int& nspheres, double& volume, double& area)
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


/* EXEMPLE OF A STANDALONE CPP APP 


int main(int argc, char *argv[])
{
  if(argc < 2) return -1;
  std::ifstream in(argv[1]);

  std::cout << "Reading Files" << std::endl;
  std::vector< double > spheres;
  int nspheres(0);
  while(!in.eof())
  {
      std::array<double,4> sphere;
      in >> sphere[0] >> sphere[1] >> sphere[2] >> sphere[3];
      if(!in.eof())
        spheres.insert(std::end(spheres), std::begin(sphere), std::end(sphere));
        nspheres++;
  }
  in.close();

  double volume, area;

  std::cout << "Computing Volume" << std::endl;
  compute_volume_surface(spheres.data(), nspheres, volume, area);

  std::cout << "Total volume:" << volume << std::endl;
  std::cout << "Total area:  " << area << std::endl;
  return 0;
}
*/
