#ifndef IO_H
#define IO_H

#include "physical_model.h"

#include <string>
#include <thread>
#include <array>

#include "XdmfDomain.hpp"
#include "XdmfGridCollectionType.hpp"
#include "XdmfUnstructuredGrid.hpp"


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

    int current_thread;

    // data size
    const int N;

    // Number of time steps per file
    const int NTimePerFile;

    // next file number
    int NumFile;

    // Total number of time steps written (or ready to be saved)
    int step;

    PhysicalModel* physicalmodel;


public:
    void CreateFile();
    void Write(const std::string& prefix, shared_ptr<XdmfUnstructuredGrid>& data, const bool all);

    ThreadedIO(PhysicalModel& _physicalmodel, const int size);
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
    friend void WriteTask(const std::string fileName, shared_ptr<XdmfDomain>* data);

};

void WriteTask(const std::string fileName, shared_ptr<XdmfDomain>* data);
std::string to_string(const double& value);
shared_ptr<XdmfInformation> xmfFormatDouble(const std::string& name,const double& number);
shared_ptr<XdmfInformation> xmfFormatInteger(const std::string& name,const int& number);
shared_ptr<XdmfInformation> xmfFormatBool(const std::string& name, const bool& active);

shared_ptr<XdmfTopology> theTopology(void);
shared_ptr<XdmfTime> settime(const double& value);
shared_ptr<XdmfGeometry> thePositions(const std::vector<double>& formatedPositions);

#endif // IO_H

