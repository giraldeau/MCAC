
/*

Sphere.h and Sphere.cpp defines the data storage.

 * Sphere *
 This is an object representing a sphere (!) with its volume and surface
 Data can be stored in an external Aggregat for vectorization purposes

 Beyond managing its properties, it can compute
 - its distance to a point or an another sphere
 - its intersection with an another sphere (volume and surface) (TODO : NOT CORRECTLY COMPUTED)
 - detect a collision with an another sphere

 * Aggregat *
 This is container for an aggregat which is an enhanced list of spheres
 Data can be shared between multiple Aggregat

*/



#include "Sphere.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define POW2(a) ((a)*(a))
#define POW3(a) ((a)*(a)*(a))


/* #############################################################################################################
 * #################################                                       #####################################
 * #################################              AGREGATE                 #####################################
 * #################################                                       #####################################
 * #############################################################################################################*/



/* #############################################################################################################
 * ########################################### Grow all spheres ################################################
 * #############################################################################################################*/
void ListSphere::CroissanceSurface(const double dt)
{
    const int listSize = size;
    #pragma omp for simd
    for (int i = 0; i < listSize; i++)
    {
        list[i]->CroissanceSurface(dt);
    }
}

