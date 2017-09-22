#ifndef IO_H
#define IO_H

#include "physical_model.hpp"

#include <array>
#include <string>
#include <thread>

//#define WITH_HDF5

#ifdef WITH_HDF5
#include "XdmfDomain.hpp"
#include "XdmfGridCollectionType.hpp"
#include "XdmfUnstructuredGrid.hpp"
namespace DLCA{


class ThreadedIO
{
private:
    std::thread* writer;

    //0-> nothing
    //1-> file created, appending data
    //2-> writing data
    std::array<int,2> status;

    std::array<shared_ptr<XdmfDomain>,2> xmfFile;
    std::array<shared_ptr<XdmfGridCollection>,2> TimeCollection;

    PhysicalModel* physicalmodel;

    size_t current_thread;

    // Total number of time steps written (or ready to be saved)
    size_t step;

    // Number of time steps per file
    const size_t NTimePerFile;

    // data size
    const size_t N;

    // next file number
    int NumFile;

//    int padding;

public:
    void CreateFile();
    void Write(const std::string& prefix, shared_ptr<XdmfUnstructuredGrid>& data, bool all);

    ThreadedIO(PhysicalModel& _physicalmodel, size_t size);
    ~ThreadedIO() noexcept;

    /** Copy constructor */
    ThreadedIO(const ThreadedIO& other);

    /** Move constructor */
    ThreadedIO (ThreadedIO&&) noexcept; /* noexcept needed to enable optimizations in containers */

    /** Copy assignment operator */
    ThreadedIO& operator= (const ThreadedIO& other);

    /** Move assignment operator */
    ThreadedIO& operator= (ThreadedIO&& other) noexcept;

private:
    friend void WriteTask(std::string fileName, shared_ptr<XdmfDomain>* data);

};

void WriteTask(std::string fileName, shared_ptr<XdmfDomain>* data);
std::string to_string(double& value);
shared_ptr<XdmfInformation> xmfFormatDouble(std::string& name,double& number);
shared_ptr<XdmfInformation> xmfFormatInteger(std::string& name,int& number);
shared_ptr<XdmfInformation> xmfFormatBool(std::string& name, bool& active);

shared_ptr<XdmfTopology> theTopology();
shared_ptr<XdmfTime> FormatTime(const double& value);
shared_ptr<XdmfGeometry> thePositions(std::vector<double>& formatedPositions);
}  // namespace DLCA

#else //WITHOUT_HDF5
namespace DLCA{

class ThreadedIO{
public:
    ThreadedIO(PhysicalModel& _physicalmodel, size_t size){}
};
}  // namespace DLCA

#endif

#endif // IO_H

