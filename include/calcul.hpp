#ifndef CALCUL_H
#define CALCUL_H 1

#include "physical_model.hpp"
#include "aggregatList.hpp"
#include "statistics.hpp"
#include <string>

namespace MCAC{

void InitRandom();
double Random();
void Calcul(PhysicalModel&);
void Init(PhysicalModel&,
//    StatisticStorage&,
    ListAggregat&);
//int LectureSuiviTempo();
//int rechercheValTab();
//void Fermeture();
void print(std::string str);

}  // namespace MCAC


#endif // CALCUL_H
