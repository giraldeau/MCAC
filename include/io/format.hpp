#ifndef INCLUDE_IO_FORMAT_HPP_
#define INCLUDE_IO_FORMAT_HPP_ 1
#ifdef WITH_HDF5
#include "io/xmf_includes.hpp"


#define DEF_FORMATER(obj, varname, type) \
  std::vector<type> obj::Format_ ## varname () const \
{ \
  const size_t list_size = size(); \
  std::vector<type> data(list_size); \
  for (size_t i=last_saved; i<list_size;i++) \
  { \
    data[i] = list[i]->varname;\
  }\
  return data;\
}
#define DEF_FORMATER_PTR(obj, varname, type) \
  std::vector<type> obj::Format_ ## varname () const \
{ \
  const size_t list_size = size(); \
  std::vector<type> data(list_size); \
  for (size_t i=last_saved; i<list_size;i++) \
  { \
    data[i] = *list[i]->varname;\
  }\
  return data;\
}
#define DEF_FORMATER_FUNC(obj, varname, type) \
  std::vector<type> obj::Format_ ## varname () const \
{ \
  const size_t list_size = size(); \
  std::vector<type> data(list_size); \
  for (size_t i=last_saved; i<list_size;i++) \
  { \
    data[i] = list[i]->varname();\
  }\
  return data;\
}
#define POSITION_FORMATER(obj) \
  std::vector<double> obj::Format_Position () const \
{ \
  const size_t list_size = size(); \
  std::vector<double> data(list_size * 3); \
  for (size_t i=last_saved; i<list_size;i++) \
  { \
    auto pos = list[i]->GetPosition(); \
    data[3*i] = pos[0];    \
    data[3*i+1] = pos[1];  \
    data[3*i+2] = pos[2];  \
  }\
  return data;\
}
#define DEF_FORMATER_POSITION(obj) \
  std::vector<double> obj::Format_Position () const \
{ \
  const size_t list_size = size(); \
  std::vector<double> data(list_size * 3); \
  for (size_t i=last_saved; i<list_size;i++) \
  { \
    auto pos = list[i]->get_position(); \
    data[3*i] = pos[0];    \
    data[3*i+1] = pos[1];  \
    data[3*i+2] = pos[2];  \
  }\
  return data;\
}
namespace MCAC {
std::string filename(int step, size_t n);
std::string to_string(const double &value);
shared_ptr<XdmfInformation> xmf_format_double(const std::string &name, const double &number);
shared_ptr<XdmfInformation> xmf_format_integer(const std::string &name, const long &number);
shared_ptr<XdmfInformation> xmf_format_bool(const std::string &name, const bool &active);
shared_ptr<XdmfTime> format_time(const double &value);
}  // namespace MCAC

#endif // WITH_HDF5
#endif //INCLUDE_IO_FORMAT_HPP_

