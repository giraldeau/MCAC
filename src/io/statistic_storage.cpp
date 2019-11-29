#include "statistics.hpp"


#ifdef WITH_HDF5
//
//#include "aggregatList.hpp"
//
//namespace MCAC{
//
//void StatisticStorage::Save()
//{
//  Save(false);
//}
//
//void StatisticStorage::Save(bool finish)
//{
//    auto aggregats_data = SavedAggregates->get_data();
//    aggregats_data->setName("StatsAggregats");
//    aggregats_data->insert(scalar("Time", FormatTimeData()));
//    aggregats_data->insert(scalar("Index", SavedAggregates->FormatIndexData()));
//    WriterAgg->write(physicalmodel->CheminSauve / "StatsAggregats", aggregats_data, finish);
//
//    SavedAggregates->lastSaved=SavedAggregates->size();
//
//    auto spheres_data = SavedAggregates->spheres.get_data();
//    spheres_data->setName("StatsSpheres");
//    WriterSph->write(physicalmodel->CheminSauve / "StatsSpheres", spheres_data, finish);
//
//    SavedAggregates->spheres.lastSaved=SavedAggregates->spheres.size();
//}
//
//std::vector<double> StatisticStorage::FormatTimeData() const
//{
//  const size_t N = times.size();
//
//  std::vector<double> time_data(N, std::allocator<double>());
//  for (size_t i=SavedAggregates->lastSaved; i<N;i++)
//  {
//    time_data[i] = times[i];
//  }
//  return time_data;
//}

//}  // namespace MCAC


#else

#define UNUSED(expr) do { (void)(expr); } while (0)

namespace MCAC{

void StatisticStorage::save()
{
    save(false);
}

void StatisticStorage::save(bool finish)
{
    UNUSED(finish);
}

}  // namespace MCAC


#endif
