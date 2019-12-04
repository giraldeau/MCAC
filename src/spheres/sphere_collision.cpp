
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



#include "spheres/sphere_collision.hpp"
#include "spheres/sphere_distance.hpp"
#include "tools/tools.hpp"
#include <cmath>


using namespace std;
namespace MCAC {
[[gnu::pure]]  pair<bool, double> collision(const Sphere &sphere_1, const Sphere &sphere_2,
                                            const array<double, 3> &displacement) noexcept {
    array<array<double, 3>, 3> rot_mat = get_rot_mat(displacement);
    return collision_r(sphere_1, sphere_2, rot_mat);
}
[[gnu::pure]]  pair<bool, double> collision_r(const Sphere &sphere_1, const Sphere &sphere_2,
                                              const array<array<double, 3>, 3> &rot_mat) noexcept {
    /*
     * We use a change of axis system
     * the center is placed on the mobil sphere
     * the x axis is chosen to be the movement vector
     *
     * Thus we have two rotation to perform : one around z and one around y
     *
     * Then a collision is easy to detect :
     * The distance to the x axis must be less than the sum of the radius
     *
     * Moreover we are only interested in collision that can happend in the future,
     * wich is also easy to detect (x coordinate must be positive)
     *
     * Finally we consider 27 possible version of the other sphere in order to take into account
     * any periodicity effect.
    */

    double dist_contact = POW_2(sphere_1.get_radius() + sphere_2.get_radius());
    double dist = distance_2(sphere_1, sphere_2);
    double minval = 0.;
    bool collision = true;
    if (dist > dist_contact) {
        double L = sphere_1.physicalmodel->L;
        collision = false;
        minval = 10 * L;
        array<double, 3> oldpos1{sphere_1.get_position()};
        array<double, 3> oldpos2{sphere_2.get_position()};
        array<double, 3> pos{};
        for (size_t l = 0; l < 3; ++l) {
            pos[l] = sum(rot_mat[l] * (oldpos2 - oldpos1));
        }
        array<double, 3> perx{L * rot_mat[0][0],
                              L * rot_mat[1][0],
                              L * rot_mat[2][0],};
        array<double, 3> pery{L * rot_mat[0][1],
                              L * rot_mat[1][1],
                              L * rot_mat[2][1],};
        array<double, 3> perz{L * rot_mat[0][2],
                              L * rot_mat[1][2],
                              L * rot_mat[2][2],};
        for (int i = -1; i <= 1; i++) {
            for (int j = -1; j <= 1; j++) {
                for (int k = -1; k <= 1; k++) {
                    array<double, 3> tmp{pos + i * perx + j * pery + k * perz};

                    // in the future
                    if (tmp[0] < 0.) {
                        continue;
                    }
                    double dist1 = POW_2(tmp[1]) + POW_2(tmp[2]);

                    // collision is possible
                    if (dist1 > dist_contact) {
                        continue;
                    }
                    collision = true;
                    double sol = dist_contact - dist1;
                    sol = tmp[0] - sqrt(sol);
                    minval = MIN(minval, sol);
                }
            }
        }
    }
    pair<bool, double> result = {collision, minval};
    return result;
}
[[gnu::pure]]  vector<double> collisions(const Sphere &sphere_1,
                                         const SphereList &list,
                                         const array<double, 3> &displacement) noexcept {
    array<array<double, 3>, 3> rot_mat{get_rot_mat(displacement)};
    vector<double> dist_to_collision;
    for (const Sphere *sphere_2 : list) {
        pair<bool, double> suspect = collision_r(sphere_1, *sphere_2, rot_mat);
        if (suspect.first) {
            dist_to_collision.push_back(suspect.second);
        }
    }
    return dist_to_collision;
}
[[gnu::pure]] array<array<double, 3>, 3> get_rot_mat(const array<double, 3> &displacement) noexcept {
    double anglez = -atan2(displacement[1], displacement[0]);
    array<array<double, 3>, 3> rotz{};
    rotz[0] = {{cos(anglez), -sin(anglez), 0}};
    rotz[1] = {{sin(anglez), cos(anglez), 0}};
    rotz[2] = {{0, 0, 1}};
    array<double, 3> tmp{};
    for (size_t i = 0; i < 3; ++i) {
        tmp[i] = sum(rotz[i] * displacement);
    }
    double angley = atan2(tmp[2], tmp[0]);
    array<array<double, 3>, 3> roty{};
    roty[0] = {{cos(angley), 0, sin(angley)}};
    roty[1] = {{0, 1, 0}};
    roty[2] = {{-sin(angley), 0, cos(angley)}};
    array<array<double, 3>, 3> matrot{};
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            matrot[i][j] = 0;
            for (size_t k = 0; k < 3; ++k) {
                matrot[i][j] += roty[i][k] * rotz[k][j];
            }
        }
    }
    return matrot;
}
}  // namespace MCAC

