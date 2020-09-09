/*
 * arvo_c - a C version of program for calculation volume and surface of
 * molecule according to article Busa et. al.: ARVO: A Fortran package for
 * computing the solvent accessible surface area and the excluded volume of
 * overlapping spheres via analytic equations, Computer Physics Communications,
 * Vol. 165, Iss. 1, Pages 59-96.
 *
 * Author(s): Jan Busa Jr.
 * Institution: Academia Sinica
 * Version 2.0, 27.12.2011
 *
 * Version history:
 * V2.0 - program rewritten into C. Added dynamic allocation of memory, changed
 *        loading of molecules (using files generated by input_structure),
 *        changed NorthPoleTest to NorthPoleFix removing the necessity of using
 *        rotation of whole molecule (which is in some cases impossible).
 * V1.0 - initial version written in FORTRAN. Can by downloaded from
 *        http://cpc.cs.qub.ac.uk/summaries/ADUL
 *
 * [Jose] Adapted by Jose MORAN, PhD student at CORIA laboratory in France
 * for soot nanoparticles simulations.
 * Rouen, France - 20.04.2020
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "arvo/arvo.hpp"
#include <omp.h>


namespace mcac {
// Constants
const double PI = 3.14159265358979323846264;
const double EPS_NORTH_POLE = 0.0001; // accuracy in the subroutine NorthPoleFix
const double NORTH_POLE_REDUCE = 0.9999;  // reduction parameter in NorthPoleFix
const double EPS_DELTAT = 1e-12; // accuracy in the subroutine CirclesIntersection
const double EPS_ANGLE = 1e-12; // accuracy in the subroutine DeleteEqual (angles)
const double EPS_TWO_PI = 1e-12;
const double VOLUME_FACTOR = 4 * PI / 3;
const double SURFACE_FACTOR = 4 * PI;
// Variables
int maxNeighbors;
std::vector<int> ind;
int spheresNumber;    // number of spheres
std::vector<double> spheres;
std::vector<int> neighborsNumber;
std::vector<int> indexStart;
// rest will have dynamic size
std::vector<int> neighborsIndices;
std::vector<double> sphereLocal;
std::vector<double> circles;
std::vector<double> arcs;
std::vector<double> newArcs2;
std::vector<double> angles2;
std::vector<double> newAngles2;
/********************************************************************************
 * Function to fill vect with zeros
 ********************************************************************************/
void fill_zero_d(std::vector<double> &vect) {
    for (size_t it = 0; it < vect.size(); it++) {
        vect[it] = 0.0;
    }
}
void fill_zero_i(std::vector<int> &vect) {
    for (size_t it = 0; it < vect.size(); it++) {
        vect[it] = 0;
    }
}
/********************************************************************************
 * Function Initialization sets all pointers used to NULL and counters to zero
 * so we can be sure, where we start.
 ********************************************************************************/
void Initialization() {
    maxNeighbors = 0;
    ind.resize(0);
    spheres.resize(0);
    neighborsNumber.resize(0);
    indexStart.resize(0);
    neighborsIndices.resize(0);
    sphereLocal.resize(0);
    circles.resize(0);
    arcs.resize(0);
    newArcs2.resize(0);
    angles2.resize(0);
    newAngles2.resize(0);
}
/********************************************************************************
 * Function Cleanup() checks whether some memory was allocated and if yes, it
 * relases this memory and sets appropriate pointer to NULL.
 ********************************************************************************/
void Cleanup() {
    ind.resize(0);
    spheres.resize(0);
    neighborsNumber.resize(0);
    indexStart.resize(0);
    neighborsIndices.resize(0);
    sphereLocal.resize(0);
    circles.resize(0);
    arcs.resize(0);
    newArcs2.resize(0);
    angles2.resize(0);
    newAngles2.resize(0);
}
/********************************************************************************
 * [Jose] Function loads aggregate information fpAts and fills the structure spheres.
 * First it allocates fixed-size structures (spheres, indexStart,
 * neighborsNumber) and finaly it fills the structure spheres using data from
 * the external vector "sph" of length "size".
 ********************************************************************************/
void LoadProtein_external(const double sph[], const int size) {
    int idx;
    if (spheresNumber < 1) {
        return;
    }
    spheres.resize(spheresNumber * 4);
    fill_zero_d(spheres);
    indexStart.resize(spheresNumber + 1);
    fill_zero_i(indexStart);
    neighborsNumber.resize(spheresNumber);
    fill_zero_i(neighborsNumber);

    // parse aggregate array: sph
    for (idx = 0; idx < size; idx++) {
        spheres[idx] = sph[idx];
    };
}
// fraction evaluation for integral
double Fract(double A, double B, double C, double sinphi, double cosphi, double k) {
    return (-B * sinphi + C * cosphi) / pow((A + B * cosphi + C * sinphi), k);
}
/********************************************************************************
 * Computing integrals over arcs given in arc structure according to paper
 * Busa et al.
 ********************************************************************************/
void AvIntegral(const int nArcs, double &pVolume, double &pArea, const double r1, const double z1) {
    int idx;
    double t, s, r, A, B, C, S, rr, ca, sa, cb, sb, be, al;
    double vIone, vItwo, vIthree, vJone, vJtwo, vJthree, delta_vint, delta_aint;
    pVolume = 0.0;
    pArea = 0.0;
    for (idx = 0; idx < nArcs; idx++) { // cycle over all arcs
        t = circles[((int) (arcs[idx * 3])) * 4];
        s = circles[((int) (arcs[idx * 3])) * 4 + 1];
        r = circles[((int) (arcs[idx * 3])) * 4 + 2];
        A = (4.0 * r1 * r1 + t * t + s * s + r * r) / 2.0;
        B = t * r;
        C = s * r;
        S = sqrt(A * A - B * B - C * C);
        rr = r * r - A;
        if (fabs(fabs(arcs[idx * 3 + 2]) - 2.0 * PI) < EPS_TWO_PI) { // full circle arc
            vIone = 2.0 * PI / S;
            vItwo = 2.0 * PI * A / (pow(S, 3.));
            vIthree = PI * (2.0 * A * A + B * B + C * C) / (pow(S, 5.));
            vJone = PI + rr / 2.0 * vIone;
            vJtwo = (vIone + rr * vItwo) / 4.0;
            vJthree = (vItwo + rr * vIthree) / 8.0;
            delta_vint = (128.0 * vJthree * pow(r1, 7.) + 8.0 * vJtwo * pow(r1, 5.) + 2.0 * vJone * pow(r1, 3.)) / 3.0
                         - 8.0 * pow(r1, 4.) * vJtwo * (z1 + r1);
            delta_aint = 2.0 * vJone * r1 * r1;
            if (arcs[idx * 3 + 2] < 0) {
                delta_vint = -delta_vint;
                delta_aint = -delta_aint;
            }
            pVolume = pVolume + delta_vint;
            pArea = pArea + delta_aint;
        } else { // integration over arcs
            if (arcs[idx * 3 + 2] < 0) {
                al = arcs[idx * 3 + 1] + arcs[idx * 3 + 2];
                be = arcs[idx * 3 + 1];
            } else {
                be = arcs[idx * 3 + 1] + arcs[idx * 3 + 2];
                al = arcs[idx * 3 + 1];
            }
            vIone = 2.0 * (PI / 2.0 - atan(
                (A * cos((be - al) / 2.0) + B * cos((al + be) / 2.0) + C * sin((al + be) / 2.0))
                / (S * sin((be - al) / 2.0)))) / S;
            sb = sin(be);
            cb = cos(be);
            sa = sin(al);
            ca = cos(al);
            vItwo = (Fract(A, B, C, sb, cb, 1) - Fract(A, B, C, sa, ca, 1) + A * vIone) / (S * S);
            vIthree = (Fract(A, B, C, sb, cb, 2) - Fract(A, B, C, sa, ca, 2)
                       + (Fract(A, B, C, sb, cb, 1) - Fract(A, B, C, sa, ca, 1)) / A
                       + (2.0 * A * A + B * B + C * C) * vItwo / A) / (2.0 * S * S);
            vJone = ((be - al) + rr * vIone) / 2.0;
            vJtwo = (vIone + rr * vItwo) / 4.0;
            vJthree = (vItwo + rr * vIthree) / 8.0;
            delta_vint = (128.0 * vJthree * pow(r1, 7.) + 8.0 * vJtwo * pow(r1, 5.) + 2.0 * vJone * pow(r1, 3.)) / 3.0
                         - 8.0 * pow(r1, 4.) * vJtwo * (z1 + r1);
            delta_aint = 2.0 * vJone * r1 * r1;
            if (arcs[idx * 3 + 2] < 0) {
                delta_vint = -delta_vint;
                delta_aint = -delta_aint;
            }
            pVolume = pVolume + delta_vint;
            pArea = pArea + delta_aint;
        }
    }
}
// Deletion of "equal" (to some precision eps_angle) angles in sorted vector angles
int DeleteEqual(const int nAngles) {
    double angle;
    int UniqueAngles, idx;
    UniqueAngles = 0;
    newAngles2.resize(angles2.size());
    fill_zero_d(newAngles2);
    angle = angles2[0];
    newAngles2[UniqueAngles] = angle;
    UniqueAngles++;
    for (idx = 1; idx < nAngles; idx++) {
        if (fabs(angles2[idx] - angle) > EPS_ANGLE) {
            angle = angles2[idx];
            newAngles2[UniqueAngles] = angle;
            UniqueAngles++;
        }
    }
    for (idx = 0; idx < UniqueAngles; idx++) {
        angles2[idx] = newAngles2[idx];
    }
    return UniqueAngles;
}
// sorting array angles in decreasing order num_angle is the angles array length
void MyDSort(const int nAngles) {
    // here no functions SetAng are necessary
    int idx, iidx, jdx;
    double amin;
    for (idx = 0; idx < nAngles - 1; idx++) {
        iidx = idx;
        amin = angles2[idx];
        for (jdx = idx + 1; jdx < nAngles; jdx++) {
            if (amin < angles2[jdx]) {
                iidx = jdx;
                amin = angles2[jdx];
            }
        }
        if (iidx != idx) {
            angles2[iidx] = angles2[idx];
            angles2[idx] = amin;
        }
    }
}
// sorting array angles in increasing order num_angle is the angles array length
void MySort(const int nAngles) {
    // here no functions SetAng are necessary
    int idx, iidx, jdx;
    double amax;
    for (idx = 0; idx < nAngles - 1; idx++) {
        iidx = idx;
        amax = angles2[idx];
        for (jdx = idx + 1; jdx < nAngles; jdx++) {
            if (amax > angles2[jdx]) {
                iidx = jdx;
                amax = angles2[jdx];
            }
        }
        if (iidx != idx) {
            angles2[iidx] = angles2[idx];
            angles2[idx] = amax;
        }
    }
}
/********************************************************************************
 * Function PointInCircle returns
 * 1  if point (t,s) is inside k-th positive circle or outside k-th negative circle
 * 0  otherwise
 * WE KNOW, THAT POINT IS NOT ON THE CIRCLE
 ********************************************************************************/
int PointInCircle(const double t, const double s, const int Circle) {
    int PointInCircle;
    double d;
    PointInCircle = 0;
    d = sqrt((t - circles[Circle * 4]) * (t - circles[Circle * 4]) + (s - circles[Circle * 4 + 1])
                                                                     * (s - circles[Circle * 4 + 1]));
    if (d < circles[Circle * 4 + 2]) {
        PointInCircle = (circles[Circle * 4 + 3] > 0) ? 1 : 0;
    } else {
        PointInCircle = (circles[Circle * 4 + 3] > 0) ? 0 : 1;
    }
    return PointInCircle;
}
/********************************************************************************
 * Function CircleInCircle returns
 * 1  if Circ1-th circle is inside Circ2-th positive circle or outside Circ2-th negative circle
 * 0  otherwise
 * WE KNOW, THAT CIRCLES HAVE LESS THAN 2 INTERSECTION POINTS
 * i -> Circ1, k->Circ2
 ********************************************************************************/
int CircleInCircle(const int Circ1, const int Circ2) {
    int CirclesInCircle;
    double d;
    CirclesInCircle = 0;
    d = sqrt((circles[Circ1 * 4] + circles[Circ1 * 4 + 2] - circles[Circ2 * 4])
             * (circles[Circ1 * 4] + circles[Circ1 * 4 + 2] - circles[Circ2 * 4])
             + (circles[Circ1 * 4 + 1] - circles[Circ2 * 4 + 1]) * (circles[Circ1 * 4 + 1]
                                                                    - circles[Circ2 * 4 + 1]));
    if (d < circles[Circ2 * 4 + 2]) {
        CirclesInCircle = (circles[Circ2 * 4 + 3] > 0) ? 1 : 0;
    } else if (d > circles[Circ2 * 4 + 2]) {
        CirclesInCircle = (circles[Circ2 * 4 + 3] > 0) ? 0 : 1;
    } else {
        d = sqrt((circles[Circ1 * 4] - circles[Circ2 * 4]) * (circles[Circ1 * 4]
                                                              - circles[Circ2 * 4 + 0])
                 + (circles[Circ1 * 4 + 1] - circles[Circ2 * 4 + 1])
                   * (circles[Circ1 * 4 + 1] - circles[Circ2 * 4 + 1]));
        if (d < circles[Circ2 * 4 + 2]) {
            CirclesInCircle = (circles[Circ2 * 4 + 3] > 0) ? 1 : 0;
        } else {
            CirclesInCircle = (circles[Circ2 * 4 + 3] > 0) ? 0 : 1;
        }
    }
    return CirclesInCircle;
}
/********************************************************************************
 * Function CirclesIntersection returns angles of two intersection points of
 * circles with indices ic1 and ic2 in circles structure circles (we will use
 * it ONLY IN CASE, WHEN 2 INTERSECTION POINTS EXIST)
 * a1 and a2 are corresponding angles with respect to the center of 1st circ
 * b1 and b2 are corresponding angles with respect to the center of 2nd circ
 ********************************************************************************/
void CirclesIntersection(const int Circ1, const int Circ2, double &a1, double &a2) {
    //     (t,s) - circle center, r - circle radius
    double t1, s1, r1, t2, s2, r2, b1, b2, A, B, C, D;
    a1 = 0;
    a2 = 0;
    b1 = 0;
    b2 = 0;
    t1 = circles[Circ1 * 4 + 0];
    s1 = circles[Circ1 * 4 + 1];
    r1 = circles[Circ1 * 4 + 2];
    t2 = circles[Circ2 * 4 + 0];
    s2 = circles[Circ2 * 4 + 1];
    r2 = circles[Circ2 * 4 + 2];
    if (fabs(t2 - t1) < EPS_DELTAT) { // t2 == t1
        B = ((r1 * r1 - r2 * r2) / (s2 - s1) - (s2 - s1)) / 2.0;
        A = sqrt(r2 * r2 - B * B);
        if (B == 0) {
            b1 = 0.0;
            b2 = PI;
        } else if (B > 0) {
            b1 = atan(fabs(B / A));
            b2 = PI - b1;
        } else {
            b1 = PI + atan(fabs(B / A));
            b2 = 3.0 * PI - b1;
        }
        B = B + s2 - s1;
        if (B == 0) {
            a1 = 0.0;
            a2 = PI;
        } else if (B > 0) {
            a1 = atan(fabs(B / A));
            a2 = PI - a1;
        } else {
            a1 = PI + atan(fabs(B / A));
            a2 = 3.0 * PI - a1;
        }
    } else { // t2 != t1
        C = ((r1 * r1 - r2 * r2 - (s2 - s1) * (s2 - s1)) / (t2 - t1) - (t2 - t1)) * 0.5;
        D = (s1 - s2) / (t2 - t1);
        B = (-C * D + sqrt((D * D + 1.0) * r2 * r2 - C * C)) / (D * D + 1.0);
        A = C + D * B;
        if (A == 0) {
            b1 = (B > 0) ? PI * 0.5 : -PI * 0.5;
        } else if (A > 0) {
            b1 = atan(B / A);
        } else {
            b1 = PI + atan(B / A);
        }
        B = B + s2 - s1;
        A = A + t2 - t1;
        if (A == 0) {
            a1 = (B > 0) ? PI * 0.5 : -PI * 0.5;
        } else if (A > 0) {
            a1 = atan(B / A);
        } else {
            a1 = PI + atan(B / A);
        }
        B = (-C * D - sqrt((D * D + 1.0) * r2 * r2 - C * C)) / (D * D + 1.0);
        A = C + D * B;
        if (A == 0) {
            b2 = (B > 0) ? PI * 0.5 : -PI * 0.5;
        } else if (A > 0) {
            b2 = atan(B / A);
        } else {
            b2 = PI + atan(B / A);
        }
        B = B + s2 - s1;
        A = A + t2 - t1;
        if (A == 0) {
            a2 = (B > 0) ? PI * 0.5 : -PI * 0.5;
        } else if (A > 0) {
            a2 = atan(B / A);
        } else {
            a2 = PI + atan(B / A);
        }
    }
    if (a1 < 0) {
        a1 = a1 + 2.0 * PI;
    }
    if (a2 < 0) {
        a2 = a2 + 2.0 * PI;
    }
    if (b1 < 0) {
        b1 = b1 + 2.0 * PI;
    }
    if (b2 < 0) {
        b2 = b2 + 2.0 * PI;
    }
}
/********************************************************************************
 * Function NewArcs prepares arcs, which are part of i-th circle in circle
 * structure circles. Interesting are these arcs, which are inside other
 * positive circles or outside other negative circles
 *
 * Matrix arcsnew in each row has elements
 *
 * arcsnew(i,1)=ic - ic is the index of arc-circle in circle
 * arcsnew(i,2)=sigma - sigma is the starting angle of arc
 * arcsnew(i,3)=delta - delta is oriented arc angle
 ********************************************************************************/
int NewArcs(const int Circle, const int NumAtoms) {
    int numArc, numAngle, numCond, idx, jdx;
    double ti, si, ri, t, s, r, d, a1, a2;
    numArc = 0;
    numAngle = 0;
    numCond = 0;
    angles2.resize(0);
    ti = circles[Circle * 4 + 0];
    si = circles[Circle * 4 + 1];
    ri = circles[Circle * 4 + 2];
    for (idx = 0; idx < NumAtoms; idx++) { // composition of angles vector, consisting of intersection points
        if (idx != Circle) {
            t = circles[idx * 4];
            s = circles[idx * 4 + 1];
            r = circles[idx * 4 + 2];
            d = sqrt((ti - t) * (ti - t) + (si - s) * (si - s));
            if ((d < r + ri) && (fabs(r - ri) < d)) { // 2 intersection points exist
                CirclesIntersection(Circle, idx, a1, a2);
                angles2.push_back(a1);
                angles2.push_back(a2);
                numAngle += 2;
            }
        }
    }
    if (!numAngle) { // there are no double intersections of idx-th circles with
        // others
        numCond = 0;
        for (idx = 0; idx < NumAtoms; idx++) { // if i-th circle is inside of all other
            // positive and outside of all other
            // negative circles, it will be new arc
            if (idx != Circle) {
                numCond = numCond + CircleInCircle(Circle, idx);
            }
        }
        if (numCond == (NumAtoms - 1)) { // all conditions hold
            newArcs2[3 * numArc + 0] = Circle;
            newArcs2[3 * numArc + 1] = 0.0;
            newArcs2[3 * numArc + 2] = 2.0 * PI * circles[Circle * 4 + 3];
            numArc++;
        }
    } else { // there are double intersection points
        if (circles[Circle * 4 + 3] > 0) {
            MySort(numAngle);
        } else {
            MyDSort(numAngle);
        }
        numAngle = DeleteEqual(numAngle);
        for (idx = 0; idx < numAngle - 1; idx++) {
            numCond = 0;
            for (jdx = 0; jdx < NumAtoms; jdx++) {
                if (jdx != Circle) {
                    t = ti + ri * cos((angles2[idx] + angles2[idx + 1]) / 2.0);
                    s = si + ri * sin((angles2[idx] + angles2[idx + 1]) / 2.0);
                    numCond = numCond + PointInCircle(t, s, jdx);
                }
            }
            if (numCond == (NumAtoms - 1)) { // all conditions hold
                newArcs2[3 * numArc + 0] = Circle;
                newArcs2[3 * numArc + 1] = angles2[idx];
                newArcs2[3 * numArc + 2] = angles2[idx + 1] - angles2[idx];
                numArc++; // zero based indices
            }
        }
        numCond = 0;
        for (idx = 0; idx < NumAtoms; idx++) {
            if (idx != Circle) {
                t = ti + ri * cos((angles2[0] + 2.0 * PI + angles2[numAngle - 1]) / 2.0);
                s = si + ri * sin((angles2[0] + 2.0 * PI + angles2[numAngle - 1]) / 2.0);
                numCond = numCond + PointInCircle(t, s, idx);
            }
        }
        if (numCond == (NumAtoms - 1)) { // all conditions hold
            newArcs2[3 * numArc + 0] = Circle;
            newArcs2[3 * numArc + 1] = angles2[numAngle - 1];
            newArcs2[3 * numArc + 2] = angles2[0] + circles[Circle * 4 + 3] * 2.0 * PI - angles2[numAngle - 1];
            numArc++;
        }
    }
    return numArc;
}
/********************************************************************************
 * Function CirclesToArcs computes integration arcs
 *
 * arcs(i,1)=ci     - corresponding circle index
 * arcs(i,2)=sigma  - starting arc angle
 * arcs(i,3)=delta  - oriented arc angle
 * Arcs (with their orientation) are parts of circles, which bounds are
 * circles intersection points. If the center of arc lies inside all other
 * positive and outside all other negative circles, then we will put it
 * inside arcs structure
 ********************************************************************************/
int CirclesToArcs(const int NumAtoms) {
    int nArcs(0), idx(0), jdx(0), kdx(0), nna(0);
    arcs.resize(NumAtoms * 30);
    fill_zero_d(arcs);
    newArcs2.resize(NumAtoms * 30);
    fill_zero_d(newArcs2);
    if (NumAtoms == 1) { // we have only 1 circle
        nArcs = 1;
        arcs[0 * 3 + 0] = 0;
        arcs[0 * 3 + 1] = 0;
        arcs[0 * 3 + 2] = 2.0 * PI * circles[3];
    } else { // more than 1 circle
        for (idx = 0; idx < NumAtoms; idx++) {
            nna = NewArcs(idx, NumAtoms);
            if (nna) {
                for (jdx = 0; jdx < nna; jdx++) {
                    for (kdx = 0; kdx < 3; kdx++) {
                        arcs[(nArcs + jdx) * 3 + kdx] = newArcs2[jdx * 3 + kdx];
                    }
                }
                nArcs += nna;
            }
        }
    }
    return nArcs;
}
/********************************************************************************
 * Function MakeTsCircles prepares circles structure for 1st sphere
 * in array circles according to the paper Busa et al.
 *
 *   circles(i,1)=ti
 *   circles(i,2)=si    - ith circle's center coordinates
 *   circles(i,3)=ri    - ith circle's radius
 *   circles(i,4)=+1/-1 - circle orientation
 ********************************************************************************/
void MakeTsCircles(const int NumAtoms) {
    double r1, dx, dy, a, b, c, d;
    int idx;
    circles.resize(NumAtoms * 4);
    fill_zero_d(circles);
    r1 = sphereLocal[3];
    for (idx = 0; idx < NumAtoms; idx++) {
        dx = sphereLocal[0] - sphereLocal[(idx + 1) * 4];
        dy = sphereLocal[1] - sphereLocal[(idx + 1) * 4 + 1];
        a = dx * dx + dy * dy + (sphereLocal[2] + r1 - sphereLocal[(idx + 1) * 4 + 2])
                                * (sphereLocal[2] + r1 - sphereLocal[(idx + 1) * 4 + 2])
            - sphereLocal[(idx + 1) * 4 + 3]
              * sphereLocal[(idx + 1) * 4 + 3];
        b = 8.0 * r1 * r1 * dx;
        c = 8.0 * r1 * r1 * dy;
        d = 4.0 * r1 * r1 * (dx * dx + dy * dy + (sphereLocal[2] - r1 - sphereLocal[(idx + 1) * 4 + 2])
                                                 * (sphereLocal[2] - r1 - sphereLocal[(idx + 1) * 4 + 2])
                             - sphereLocal[(idx + 1) * 4 + 3]
                               * sphereLocal[(idx + 1) * 4 + 3]);
        circles[idx * 4 + 0] = -b / (2.0 * a);
        circles[idx * 4 + 1] = -c / (2.0 * a);
        circles[idx * 4 + 2] = sqrt((b * b + c * c - 4.0 * a * d) / (4.0 * a * a));
        if (a > 0) {
            circles[idx * 4 + 3] = -1;
        } else {
            circles[idx * 4 + 3] = 1;
        }
    }
}
/********************************************************************************
 * Function NorthPoleFix checks on local spheres if north pole of sphere
 * doesn't lie on the surface of other sphere. If it does, the radius of
 * sphere is multiplied by factor of NORTH_POLE_REDUCE. This function replaces
 * the NorthPoleTest and SpheresRotation from original algorithm.
 ********************************************************************************/
void NorthPoleFix(const int numAtoms) {
    double d;
    int idx;
    for (idx = 1; idx < numAtoms; idx++) { // all except atom 0
        d = fabs(
            sqrt(
                (sphereLocal[0] - sphereLocal[idx * 4])
                * (sphereLocal[0] - sphereLocal[idx * 4])
                + (sphereLocal[1] - sphereLocal[idx * 4 + 1])
                  * (sphereLocal[1] - sphereLocal[idx * 4 + 1])
                + (sphereLocal[2] + sphereLocal[3] - sphereLocal[idx * 4 + 2])
                  * (sphereLocal[2] + sphereLocal[3]
                     - sphereLocal[idx * 4 + 2]))
            - sphereLocal[idx * 4 + 3]);
        if (d < EPS_NORTH_POLE) {
            sphereLocal[idx * 4 + 3] = sphereLocal[idx * 4 + 3] * NORTH_POLE_REDUCE;
            idx--;
            //printf("Reduced sphere radius\n"); // print information that radius was
            // reduced
        }
    }
}
/********************************************************************************
 * Function LocalSpheres copies spheres whose indices are in the array ind
 * (actual atom and its neighbors) into the array sphereLocal and calls the
 * function NorthPoleFix to avoid that the north pole of actual atom lies on
 * some sphere.
 ********************************************************************************/
void LocalSpheres(const int numAtoms) {
    sphereLocal.resize(numAtoms * 4);
    fill_zero_d(sphereLocal);
    int idx, jdx, kdx = 0;
    for (idx = 0; idx < numAtoms; idx++) {
        for (jdx = 0; jdx < 4; jdx++) {
            sphereLocal[kdx++] = spheres[ind[idx] * 4 + jdx];
        }
    }
    NorthPoleFix(numAtoms);
}
/********************************************************************************
 * Function AreaVolume computes Atom-th part of the whole volume - the volume
 * of domain inside Atom-th and outside of all other spheres
 ********************************************************************************/
void AreaVolume(const int Atom, double &pVolume, double &pArea) {
    // nls was originaly neighborsNumber[Atom]+1, shifted to neighborsNumber[Atom]!!!
    pVolume = 0;
    pArea = 0;

    // Determination of i-th sphere's neighbors (row indices in matrix spheres)
    if (neighborsNumber[Atom] < 0) {
        // ith sphere is subset of other sphere, sphere(i,4) will be done negative
        pVolume = 0.0;
        pArea = 0.0;
    } else if (neighborsNumber[Atom] == 0) {
        // there are no neighbors (nls - number of local spheres = ith sphere + neigh
        pVolume = VOLUME_FACTOR * spheres[Atom * 4 + 3] * spheres[Atom * 4 + 3] * spheres[Atom * 4 + 3];
        pArea = SURFACE_FACTOR * spheres[Atom * 4 + 3] * spheres[Atom * 4 + 3];
    } else {
        int idx, nArcs, nPos;
        double z1, r1, partVol, partArea;
        // there are neighbors
        ind.resize(neighborsNumber[Atom] + 1);
        fill_zero_i(ind);
        ind[0] = Atom;
        for (idx = 0; idx < neighborsNumber[Atom]; idx++) {
            ind[idx + 1] = neighborsIndices[indexStart[Atom] + idx];
        }

        // we will work only with ith and neighbors spheres
        LocalSpheres(neighborsNumber[Atom] + 1);
        pVolume = 0.0;
        pArea = 0.0;
        MakeTsCircles(neighborsNumber[Atom]);
        nArcs = CirclesToArcs(neighborsNumber[Atom]);
        nPos = 0;
        for (idx = 0; idx < neighborsNumber[Atom]; idx++) {
            if (circles[idx * 4 + 3] > 0) {
                nPos++;
            }
        }
        z1 = sphereLocal[2];
        r1 = sphereLocal[3];
        AvIntegral(nArcs, partVol, partArea, r1, z1);
        if (nPos) { // there exists positive oriented circle
            pVolume = pVolume + partVol;
            pArea = pArea + partArea;
        } else { // all circles are negative oriented - we have computed complement
            pVolume = pVolume + partVol + VOLUME_FACTOR * sphereLocal[3] * sphereLocal[3] * sphereLocal[3];
            pArea = pArea + partArea + SURFACE_FACTOR * sphereLocal[3] * sphereLocal[3];
        }
    }
}
/********************************************************************************
 * If ith sphere is a subset of other sphere, index_number(i)=-1 and we change
 * radius in matrix spheres to -radius.
 * If some other sphere is subset of ith sphere, than we change its radius to
 * -radius.
 ********************************************************************************/
int Neighbors(const int Atom, const std::vector<std::unordered_map<size_t, double>> distances) {
    int NeighborsNum(0), idx(0); // number of neighbors
    NeighborsNum = distances[Atom].size();
    ind.resize(NeighborsNum);
    fill_zero_i(ind);
    for (const std::pair<size_t, double> &it : distances[Atom]) {
        ind[idx] = it.first;
        idx += 1;
    }
    return NeighborsNum;
}
/********************************************************************************
 * Determination of neighbors for all atoms. We construct next structure:
 * neighborsNumber(i) = neighbors number for ith atom
 * indexStart(i)      = start of neighbors indices for ith atom in array
 *                      neighbors_indices
 * neighborsIndices   = array of neighbors indices for each atom
 * neighborsIndices(indexStart(i)):neighborsInd(indexStart(i)+neighborsNumber(i)
 *
 * For example: 1. atom has neighbors with indices 2, 4, 7
 *              2. atom has neighbors with indices 1, 3
 *              3. atom has neighbors with indices 2, 4
 *              4. atom has neighbors with indices 1, 3
 *              5. atom is subset of some atom
 *              6. atom has no neighbors
 *              7. atom has neighbors with index 1
 * then we have:
 *               neighborsNumber = (3,2,2,2,-1,0,1)
 *                    indexStart = (0,3,5,7,9,9,9)
 *              neighborsIndices = (1,3,6,0,2,1,3,0,2,0)
 ********************************************************************************/
#ifdef FULL_INTERNAL_DISTANCES
void MakeNeighbors(const std::vector<std::unordered_map<size_t, double>> &distances) {
    int idx, jdx, total_neigh(0);
    indexStart[0] = 0;
    maxNeighbors = 0;
    neighborsIndices.resize(0);
    for (idx = 0; idx < spheresNumber; idx++) {
        neighborsNumber[idx] = Neighbors(idx, distances);
        if (neighborsNumber[idx] > maxNeighbors) {
            maxNeighbors = neighborsNumber[idx];
        }
        if (neighborsNumber[idx] <= 0) { // sphere is subset or there are no neighbors
            indexStart[idx + 1] = indexStart[idx];
        } else { //  there are neighbors
            indexStart[idx + 1] = indexStart[idx] + neighborsNumber[idx];
            for (jdx = 0; jdx < neighborsNumber[idx]; jdx++) {
                total_neigh += 1;
                neighborsIndices.push_back(ind[jdx]);
            }
        }
    }
}
#else
void MakeNeighbors(const std::vector<std::unordered_map<size_t, double>> &distances) {
    indexStart[0] = 0;
    maxNeighbors = 0;
    neighborsIndices.resize(0);
    for (size_t idx = 0; idx < distances.size(); idx++) {
        neighborsNumber[idx] = distances[idx].size();
        if (neighborsNumber[idx] > maxNeighbors) {
            maxNeighbors = neighborsNumber[idx];
        }
        indexStart[idx + 1] = indexStart[idx] + neighborsNumber[idx];
        for (const auto& neighbor: distances[idx]) {
            neighborsIndices.push_back(neighbor.first);
        }
    }
}
#endif
/********************************************************************************
 * Function PrintUsage provides user with short information on how to use
 * the program arvo.
 ********************************************************************************/
void PrintUsage() {
    printf("ARVO_C 2.0.\n\nUsage: arvo protein=file1 [log=file2]\n\n");
    printf("Mandatory:\n  protein - name of input file as created by input_structure\n");
    printf("            (you can also use short form 'p=')\n");
    printf("Optional:\n");
    printf("      log - name of file for results and log messages\n");
}
/********************************************************************************
 * [Jose] Main body of the program ARVO. It parses the parameters, loads agg. info,
 * creates the neighbors structure, calls function for calculating the
 * contribution to the area and volume of every atom and collects the results.
*********************************************************************************/
void main_arvo(const double sph[],
               double vol[],
               double surf[],
               const size_t nspheres,
               const std::vector<std::unordered_map<size_t, double>> &distances) {
    Initialization(); // initializate the memory

    LoadProtein_external(sph, nspheres * 4); // load protein and fill structure spheres
    spheresNumber = nspheres;

    if (spheresNumber <= 0) {
        printf("No spheres to calculate volume. Exiting.\n");
        Cleanup();
        return;
    }

    // Study the neighborhood relations [to be improved]
    MakeNeighbors(distances);

    // Computation of area and volume as a sum of surface integrals
    for (int idx = 0; idx < spheresNumber; idx++) {
        AreaVolume(idx, vol[idx], surf[idx]);
    }
    Cleanup();
    return;
}
}