#include "Spherelist.h"
#include "physical_model.h"
#include "aggregatList.h"
#include "IO.h"

#include <iostream>
#include <iomanip>
#include <sstream>

#include "XdmfInformation.hpp"
#include "XdmfWriter.hpp"
#include "XdmfHDF5Writer.hpp"
#include "XdmfGridCollection.hpp"

/*
#include "XdmfSystemUtils.hpp"
#include "XdmfAttribute.hpp"
#include "XdmfCurvilinearGrid.hpp"
#include "XdmfRectilinearGrid.hpp"
#include "XdmfRegularGrid.hpp"
#include "XdmfMap.hpp"
#include "XdmfAttributeType.hpp"
#include "XdmfAttributeCenter.hpp"
#include "XdmfSet.hpp"
#include "XdmfArray.hpp"
#include "XdmfGeometry.hpp"
#include "XdmfTopology.hpp"
*/


// Usefull tool
std::string to_string(const double& value) {
    std::stringstream sstr;
    sstr << value;
    return sstr.str();
}


// Format for Physical Model
shared_ptr<XdmfInformation> xmfFormatDouble(const std::string& name,const double& number)
{
    shared_ptr<XdmfInformation> info = XdmfInformation::New(name, to_string(number));
    return info;
}

shared_ptr<XdmfInformation> xmfFormatInteger(const std::string& name,const int& number)
{
    shared_ptr<XdmfInformation> info = XdmfInformation::New(name, std::to_string(number));
    return info;
}

shared_ptr<XdmfInformation> xmfFormatBool(const std::string& name, const bool& active)
{
    shared_ptr<XdmfInformation> info;
    if (active)
        info = XdmfInformation::New(name, "True");
    else
        info = XdmfInformation::New(name, "False");
    return info;
}


// Format for ListSpheres
std::vector<double> ListSphere::FormatPositionData() const
{
    const int N = size();

    std::vector<double> PositionData(3*N);
    for (int i=0; i<N;i++)
    {
        std::array<double,3> pos = list[i]->Position();
        PositionData[3*i]=pos[0];
        PositionData[3*i+1]=pos[1];
        PositionData[3*i+2]=pos[2];
    }
    return PositionData;
}

std::vector<double> ListSphere::FormatRadiusData() const
{
    return (*Storage)[3];
}


std::vector<int> ListSphere::FormatLabelData() const
{
    const int N = size();

    std::vector<int> LabelData(N);
    for (int i=0; i<N;i++)
    {
        LabelData[i] = list[i]->AggLabel;
    }
    return LabelData;
}


// Format for ListAggretates
std::vector<double> ListAggregat::FormatPositionData() const
{
    const int N = size();

    std::vector<double> PositionData(3*N);
    for (int i=0; i<N;i++)
    {
        std::array<double,3> pos = list[i]->GetPosition();
        PositionData[3*i]=pos[0];
        PositionData[3*i+1]=pos[1];
        PositionData[3*i+2]=pos[2];
    }
    return PositionData;

}
std::vector<double> ListAggregat::FormatRgData() const
{
    return (*Storage)[0];
}

std::vector<int>    ListAggregat::FormatNpData() const
{
    const int N = size();

    std::vector<int> NpData(N);
    for (int i=0; i<N;i++)
    {
        NpData[i] = list[i]->Np;
    }
    return NpData;
}
std::vector<double> ListAggregat::FormatDmData() const
{
    return (*Storage)[1];
}
std::vector<double> ListAggregat::FormatlpmData() const
{
    return (*Storage)[2];
}
std::vector<double> ListAggregat::FormatdeltatData() const
{
    return (*Storage)[3];
}
std::vector<double> ListAggregat::FormatRmaxData() const
{
    return (*Storage)[4];
}
std::vector<double> ListAggregat::FormatVoumeData() const
{
    return (*Storage)[5];
}
std::vector<double> ListAggregat::FormatSurfaceData() const
{
    return (*Storage)[6];
}
std::vector<int>    ListAggregat::FormatLabelData() const
{
    const int N = size();

    std::vector<int> LabelData(N);
    for (int i=0; i<N;i++)
    {
        LabelData[i] = list[i]->indexInStorage;
    }
    return LabelData;

}

// Format filename
std::string filename(const int step, const int N)
{
    int witdh=int(ceil(log10(double(N))))+4;

    std::ostringstream fileNameStream;
    fileNameStream << "_" << std::setfill('0') << std::setw(witdh) << step;

    return fileNameStream.str();
}


// Actual write functions
shared_ptr<XdmfTopology> theTopology(void)
{
    shared_ptr<XdmfTopology> particules = XdmfTopology::New();
    particules->setType(XdmfTopologyType::Polyvertex());
    return particules;
}

shared_ptr<XdmfTime> time(const double& value)
{
    shared_ptr<XdmfTime> thetime = XdmfTime::New(value);
    return thetime;
}

shared_ptr<XdmfGeometry> thePositions(const std::vector<double>& formatedPositions)
{
    const int N = int(formatedPositions.size());
    shared_ptr<XdmfGeometry> Positions = XdmfGeometry::New();
    Positions->setType(XdmfGeometryType::XYZ());
    Positions->insert(0, formatedPositions.data(), N, 1, 1);
    return Positions;
}

template <class T>
shared_ptr<XdmfAttribute> Scalar(const std::string& name, const std::vector<T> formatedField)
{
    const int N = int(formatedField.size());

    shared_ptr<XdmfAttribute> xdmfField = XdmfAttribute::New();
    xdmfField->setName(name);
    xdmfField->setType(XdmfAttributeType::Scalar());
    xdmfField->setCenter(XdmfAttributeCenter::Node());
    xdmfField->insert(0, formatedField.data(), N, 1, 1);
    return xdmfField;
}

auto PhysicalModel::xmfWrite() const
{

    shared_ptr<XdmfInformation> info = XdmfInformation::New("Physics", "Physical properties of the simulation");

    info->insert(xmfFormatDouble("Asurfgrowth",Asurfgrowth));
    info->insert(xmfFormatDouble("dfe",dfe));
    info->insert(xmfFormatDouble("kfe",kfe));
    info->insert(xmfFormatDouble("xsurfgrowth",xsurfgrowth));
    info->insert(xmfFormatDouble("coeffB",coeffB));
    info->insert(xmfFormatDouble("lambda",lambda));
    info->insert(xmfFormatDouble("Dpeqmass",Dpeqmass));
    info->insert(xmfFormatDouble("rpeqmass",rpeqmass));
    info->insert(xmfFormatDouble("gamma_",gamma_));
    info->insert(xmfFormatDouble("P [Pa]",P));
    info->insert(xmfFormatDouble("T [K]",T));
    info->insert(xmfFormatDouble("Mu",Mu));
    info->insert(xmfFormatDouble("K",K));
    info->insert(xmfFormatDouble("Rho [kg/m3]",Rho));
    info->insert(xmfFormatDouble("Dpm [nm]",Dpm));
    info->insert(xmfFormatDouble("sigmaDpm [nm]",sigmaDpm));
    info->insert(xmfFormatDouble("X []",X));
    info->insert(xmfFormatDouble("FV [ppt]",FV));
    info->insert(xmfFormatDouble("L",L));
    info->insert(xmfFormatInteger("N []",N));
    info->insert(xmfFormatBool("ActiveModulephysique",ActiveModulephysique));
    info->insert(xmfFormatBool("ActiveVariationTempo",ActiveVariationTempo));

    return info;
}

void ListSphere::save(void) const
{
    save(false);
}


void ListSphere::save(const bool finish) const
{
    // Set geometry
    shared_ptr<XdmfUnstructuredGrid> SpheresData = XdmfUnstructuredGrid::New();
    SpheresData->setName("Spheres");
    SpheresData->setTopology(theTopology());

    // Set time
    SpheresData->setTime(time(physicalmodel->time));

    // Set Positions
    SpheresData->setGeometry(thePositions(FormatPositionData()));

    // Set Radius
    SpheresData->insert(Scalar("Radius",FormatRadiusData()));
    SpheresData->insert(Scalar("Label",FormatLabelData()));

    // Write data
    Writer->Write("Spheres",SpheresData, finish);
}

void ListAggregat::save(void) const
{
    save(false);
}


void ListAggregat::save(const bool finish) const
{

    // Set geometry
    shared_ptr<XdmfUnstructuredGrid> AggregatsData = XdmfUnstructuredGrid::New();
    AggregatsData->setName("Aggregats");
    shared_ptr<XdmfTopology> particules = XdmfTopology::New();
    particules->setType(XdmfTopologyType::Polyvertex());
    AggregatsData->setTopology(particules);

    // Set time
    shared_ptr<XdmfTime> time = XdmfTime::New(physicalmodel->time);
    AggregatsData->setTime(time);

    // Set Positions
    AggregatsData->setGeometry(thePositions(FormatPositionData()));

    // Set Radius
    AggregatsData->insert(Scalar("Rg",FormatRgData()));
    AggregatsData->insert(Scalar("Np",FormatNpData()));
    AggregatsData->insert(Scalar("Dm",FormatDmData()));
    AggregatsData->insert(Scalar("lpm",FormatlpmData()));
    AggregatsData->insert(Scalar("Deltat",FormatdeltatData()));
    AggregatsData->insert(Scalar("Rmax",FormatRmaxData()));
    AggregatsData->insert(Scalar("Volume",FormatVoumeData()));
    AggregatsData->insert(Scalar("Surface",FormatSurfaceData()));
    AggregatsData->insert(Scalar("Label",FormatLabelData()));

    // Write data
    Writer->Write("Aggregats",AggregatsData, finish);
}


void ThreadedIO::CreateFile()
{

        if (status[current_thread]==2)
        {
            writer->join();
            status[current_thread] = 0;
        }
        if (status[current_thread]==1)
        {
            std::cout << "And what should happen to the old file ?" << std::endl;
            exit(5);
        }
        if (status[current_thread]==0)
        {
            // Prepare Xmf file
            xmfFile[current_thread] = XdmfDomain::New();
            shared_ptr<XdmfInformation> xmfInfo = XdmfInformation::New("Copyright", "Produced by DLCA");
            xmfFile[current_thread]->insert(xmfInfo);

            // Save simulation properties in every xmf file
            xmfFile[current_thread]->insert(physicalmodel->xmfWrite());

            TimeCollection[current_thread] = XdmfGridCollection::New();
            TimeCollection[current_thread]->setType(XdmfGridCollectionType::Temporal());

            status[current_thread] = 1;
        }

}

void ThreadedIO::Write(const std::string& prefix, shared_ptr<XdmfUnstructuredGrid>& data, bool all)
{
    if (step%NTimePerFile == 0)
    {
        CreateFile();
    }
    if (status[current_thread] != 1)
    {
        std::cout << "File is not ready..."<< std::endl;
        exit(5);
    }
    TimeCollection[current_thread]->insert(data);

    step++;

    if(step%NTimePerFile == 0 || all)
    {
        xmfFile[current_thread]->insert(TimeCollection[current_thread]);

        std::string fileName = prefix + filename(NumFile, N);

        if (status[!current_thread] == 2)
        {
            writer->join();
            status[!current_thread] = 0;
        }


        writer = new std::thread(WriteTask, fileName, &xmfFile[current_thread]);
        status[current_thread] = 2;

        current_thread = !current_thread;

        NumFile++;

        if (all)
        {
            step=0;
            writer->join();
        }
    }

}

void WriteTask(const std::string fileName, shared_ptr<XdmfDomain>* data)
{
    shared_ptr<XdmfHDF5Writer> HDF5File = XdmfHDF5Writer::New(fileName+".h5");
    shared_ptr<XdmfWriter> XMFFile = XdmfWriter::New(fileName+".xmf", HDF5File);

    // Write data
    (*data)->accept(HDF5File);
    HDF5File->setMode(XdmfHeavyDataWriter::Overwrite);
    (*data)->accept(XMFFile);
}



ThreadedIO::ThreadedIO(PhysicalModel& _physicalmodel, int size):
    writer(),
    status({0,0}),
    xmfFile(),
    TimeCollection(),
    current_thread(0),
    N(size),
    NTimePerFile(size/10),
    NumFile(0),
    step(0),
    physicalmodel(&_physicalmodel)
{}


ThreadedIO::~ThreadedIO(void)
{
    if (status[current_thread]==1 || status[!current_thread]==1)
    {
        std::cout << "And what should happen to the data ?" << std::endl;
        exit(5);
    }
    if (status[current_thread]==2 || status[!current_thread]==2)
    {
        writer->join();
    }
}
