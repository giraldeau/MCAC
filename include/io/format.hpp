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
#ifndef INCLUDE_IO_FORMAT_HPP
#define INCLUDE_IO_FORMAT_HPP 1
#ifdef WITH_HDF5
#include "io/xmf_includes.hpp"


#define DEF_FORMATER(obj, varname, type) \
  std::vector<type> obj::format_ ## varname () const \
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
  std::vector<type> obj::format_ ## varname () const \
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
  std::vector<type> obj::format_ ## varname () const \
{ \
  const size_t list_size = size(); \
  std::vector<type> data(list_size); \
  for (size_t i=last_saved; i<list_size;i++) \
  { \
    data[i] = list[i]->varname();\
  }\
  return data;\
}
#define DEF_FORMATER_POSITION(obj) \
  std::vector<double> obj::format_position () const \
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
namespace mcac {
std::string filename(int step, size_t n);
template<class T>
std::string to_string(const T &value);
template<class T>
boost::shared_ptr<XdmfInformation> xmf_format(const std::string &name, const T &number);
boost::shared_ptr<XdmfInformation> xmf_format(const std::string &name, const bool &active);
boost::shared_ptr<XdmfTime> format_time(const double &value);
}  // namespace mcac

#endif // WITH_HDF5
#endif //INCLUDE_IO_FORMAT_HPP

