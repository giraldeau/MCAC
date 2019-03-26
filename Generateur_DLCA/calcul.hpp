#ifndef CALCUL_H
#define CALCUL_H

#include "physical_model.hpp"
#include "aggregatList.hpp"
#include "statistics.hpp"
#include <string>

namespace DLCA{

void InitRandom();
double Random();
void Calcul(PhysicalModel&);
void Init(PhysicalModel&,StatisticStorage&, ListAggregat&);
//int LectureSuiviTempo();
//int rechercheValTab();
//void Fermeture();
void print(std::string str);

}  // namespace DLCA


#endif // CALCUL_H
