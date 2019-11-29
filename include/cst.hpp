//
// Created by pouxa on 27/11/2019.
//
#ifndef INCLUDE_CST_HPP_
#define INCLUDE_CST_HPP_ 1

#define AGG_N_FIELD 15
#define SPH_N_FIELD 9


#ifdef WITH_HDF5
#define XMF_OUTPUT shared_ptr<XdmfInformation>

#else
#define XMF_OUTPUT void
#endif



namespace MCAC {

enum class Writer_status { Idle, Appending, Writing };


}  //namespace MCAC
#endif //INCLUDE_CST_HPP_
