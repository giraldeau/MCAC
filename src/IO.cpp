#include "IO.hpp"
#include "Spherelist.hpp"
#include "aggregatList.hpp"
#include "physical_model.hpp"


#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

#define UNUSED(expr) do { (void)(expr); } while (0)

#ifdef WITH_HDF5

#include "XdmfGridCollection.hpp"
#include "XdmfHDF5Writer.hpp"
#include "XdmfInformation.hpp"
#include "XdmfWriter.hpp"

namespace fs = std::experimental::filesystem;

namespace DLCA{


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
    const size_t N = size();

    std::vector<double> PositionData(3*N);
    for (size_t i=lastSaved; i<N;i++)
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
    const size_t N = size();

    std::vector<double> RadiusData(N);
    for (size_t i=lastSaved; i<N;i++)
    {
        RadiusData[i] = list[i]->Radius();
    }
    return RadiusData;
}


std::vector<long> ListSphere::FormatLabelData() const
{
    const size_t N = size();

    std::vector<long> LabelData(N);
    for (size_t i=lastSaved; i<N;i++)
    {
        LabelData[i] = list[i]->AggLabel;
    }
    return LabelData;
}


// Format for ListAggretates
std::vector<double> ListAggregat::FormatPositionData() const
{
    const size_t N = size();

    std::vector<double> PositionData(3*N);
    for (size_t i=lastSaved; i<N;i++)
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
    const size_t N = size();

    std::vector<double> RgData(N);
    for (size_t i=lastSaved; i<N;i++)
    {
        RgData[i] = *list[i]->rg;
    }
    return RgData;
}

std::vector<int>    ListAggregat::FormatNpData() const
{
    const size_t N = size();

    std::vector<int> NpData(N);
    for (size_t i=lastSaved; i<N;i++)
    {
        NpData[i] = int(list[i]->Np);
    }
    return NpData;
}
std::vector<double> ListAggregat::FormatFaggData() const
{
    const size_t N = size();

    std::vector<double> FaggData(N);
    for (size_t i=lastSaved; i<N;i++)
    {
        FaggData[i] = *list[i]->f_agg;
    }
    return FaggData;
}
std::vector<double> ListAggregat::FormatlpmData() const
{
    const size_t N = size();

    std::vector<double> lpmData(N);
    for (size_t i=lastSaved; i<N;i++)
    {
        lpmData[i] = *list[i]->lpm;
    }
    return lpmData;
}
std::vector<double> ListAggregat::FormatdeltatData() const
{
    const size_t N = size();

    std::vector<double> deltatData(N);
    for (size_t i=lastSaved; i<N;i++)
    {
        deltatData[i] = *list[i]->time_step;
    }
    return deltatData;
}
std::vector<double> ListAggregat::FormatRmaxData() const
{
    const size_t N = size();

    std::vector<double> RmaxData(N);
    for (size_t i=lastSaved; i<N;i++)
    {
        RmaxData[i] = *list[i]->rmax;
    }
    return RmaxData;
}
std::vector<double> ListAggregat::FormatVolumeData() const
{
    const size_t N = size();

    std::vector<double> VolumeData(N);
    for (size_t i=lastSaved; i<N;i++)
    {
        VolumeData[i] = *list[i]->volAgregat;
    }
    return VolumeData;
}
std::vector<double> ListAggregat::FormatSurfaceData() const
{
    const size_t N = size();

    std::vector<double> SurfaceData(N);
    for (size_t i=lastSaved; i<N;i++)
    {
        SurfaceData[i] = *list[i]->surfAgregat;
    }
    return SurfaceData;
}
std::vector<int>    ListAggregat::FormatLabelData() const
{
    const size_t N = size();

    std::vector<int> LabelData(N);
    for (size_t i=lastSaved; i<N;i++)
    {
        LabelData[i] = int(list[i]->Label);
    }
    return LabelData;

}
std::vector<int>    ListAggregat::FormatIndexData() const
{
    const size_t N = size();

    std::vector<int> IndexData(N);
    for (size_t i=lastSaved; i<N;i++)
    {
        IndexData[i] = int(list[i]->indexInStorage);
    }
    return IndexData;

}

std::vector<double> StatisticStorage::FormatTimeData() const
{
    const size_t N = times.size();

    std::vector<double> TimeData(N);
    for (size_t i=SavedAggregates->lastSaved; i<N;i++)
    {
        TimeData[i] = times[i];
    }
    return TimeData;
}


// Format filename
std::string filename(const int step, const size_t N)
{
    int witdh=int(ceil(log10(double(N))))+4;

    std::ostringstream fileNameStream;
    fileNameStream << "_" << std::setfill('0') << std::setw(witdh) << step;

    return fileNameStream.str();
}


// Actual write functions
shared_ptr<XdmfTopology> theTopology()
{
    shared_ptr<XdmfTopology> particules = XdmfTopology::New();
    particules->setType(XdmfTopologyType::Polyvertex());
    return particules;
}

shared_ptr<XdmfTime> FormatTime(const double& value)
{
    shared_ptr<XdmfTime> thetime = XdmfTime::New(value);
    return thetime;
}

shared_ptr<XdmfGeometry> thePositions(const std::vector<double>& formatedPositions)
{
    const unsigned int N= unsigned(int(formatedPositions.size()));
    shared_ptr<XdmfGeometry> Positions = XdmfGeometry::New();
    Positions->setType(XdmfGeometryType::XYZ());
    Positions->insert(0, formatedPositions.data(), N, 1, 1);
    return Positions;
}

template <class T>
shared_ptr<XdmfAttribute> Scalar(const std::string& name, const std::vector<T>& formatedField)
{
    const unsigned int N= unsigned(int(formatedField.size()));

    shared_ptr<XdmfAttribute> xdmfField = XdmfAttribute::New();
    xdmfField->setName(name);
    xdmfField->setType(XdmfAttributeType::Scalar());
    xdmfField->setCenter(XdmfAttributeCenter::Node());
    xdmfField->insert(0, formatedField.data(), N, 1, 1);
    return xdmfField;
}

template <class T>
shared_ptr<XdmfAttribute> Attribute(const std::string& name, const T& value)
{
    shared_ptr<XdmfAttribute> xdmfField = XdmfAttribute::New();
    xdmfField->setName(name);
    xdmfField->setType(XdmfAttributeType::NoAttributeType());
    xdmfField->insert(0, value);
    return xdmfField;
}

auto PhysicalModel::xmfWrite() const
{

    shared_ptr<XdmfInformation> info = XdmfInformation::New("Physics", "Physical properties of the simulation");

    const int n = int(N);

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
    info->insert(xmfFormatInteger("N []",n));
    info->insert(xmfFormatBool("ActiveModulephysique",ActiveModulephysique));
    info->insert(xmfFormatBool("ActiveVariationTempo",ActiveVariationTempo));

    return info;
}

auto ListSphere::GetData() const
{
    // Set geometry
    shared_ptr<XdmfUnstructuredGrid> SpheresData = XdmfUnstructuredGrid::New();
    SpheresData->setName("Spheres");
    SpheresData->setTopology(theTopology());

    // Set time
    SpheresData->setTime(FormatTime(physicalmodel->Time));
    SpheresData->insert(Attribute("Time",physicalmodel->Time));
    SpheresData->insert(Attribute("BoxSize",physicalmodel->L));

    // Set Positions
    SpheresData->setGeometry(thePositions(FormatPositionData()));

    // Set Radius
    SpheresData->insert(Scalar("Radius",FormatRadiusData()));
    SpheresData->insert(Scalar("Label",FormatLabelData()));

    return SpheresData;
}


void ListSphere::save()
{
    save(false);
}


void ListSphere::save(const bool finish)
{
    //get everything
    lastSaved = 0;

    auto data = GetData();

    // Write data
    Writer->Write(physicalmodel->CheminSauve / "Spheres",data, finish);
}


auto ListAggregat::GetData() const
{
    // Set geometry
    shared_ptr<XdmfUnstructuredGrid> AggregatsData = XdmfUnstructuredGrid::New();
    AggregatsData->setName("Aggregats");
    shared_ptr<XdmfTopology> particules = XdmfTopology::New();
    particules->setType(XdmfTopologyType::Polyvertex());
    AggregatsData->setTopology(particules);

    // Set time
    shared_ptr<XdmfTime> time = XdmfTime::New(physicalmodel->Time);
    AggregatsData->setTime(time);
    AggregatsData->insert(Attribute("Time",physicalmodel->Time));
    AggregatsData->insert(Attribute("BoxSize",physicalmodel->L));

    // Set Positions
    AggregatsData->setGeometry(thePositions(FormatPositionData()));

    // Set Radius
    AggregatsData->insert(Scalar("Rg",FormatRgData()));
    AggregatsData->insert(Scalar("Np",FormatNpData()));
    AggregatsData->insert(Scalar("f_agg",FormatFaggData()));
    AggregatsData->insert(Scalar("lpm",FormatlpmData()));
    AggregatsData->insert(Scalar("Deltat",FormatdeltatData()));
    AggregatsData->insert(Scalar("Rmax",FormatRmaxData()));
    AggregatsData->insert(Scalar("Volume",FormatVolumeData()));
    AggregatsData->insert(Scalar("Surface",FormatSurfaceData()));
    AggregatsData->insert(Scalar("Label",FormatLabelData()));

    return AggregatsData;
}


void ListAggregat::save()
{
    save(false);
}


void ListAggregat::save(const bool finish)
{
    //get everything
    lastSaved = 0;

    auto data = GetData();

    // Write data
    Writer->Write(physicalmodel->CheminSauve / "Aggregats",data, finish);
}


void StatisticStorage::save()
{
    save(false);
}

void StatisticStorage::save(bool finish)
{
    auto AggregatsData = SavedAggregates->GetData();
    AggregatsData->setName("StatsAggregats");
    AggregatsData->insert(Scalar("Time",FormatTimeData()));
    AggregatsData->insert(Scalar("Index",SavedAggregates->FormatIndexData()));
    WriterAgg->Write(physicalmodel->CheminSauve / "StatsAggregats",AggregatsData, finish);

    SavedAggregates->lastSaved=SavedAggregates->size();

    auto SpheresData = SavedAggregates->spheres.GetData();
    SpheresData->setName("StatsSpheres");
    WriterSph->Write(physicalmodel->CheminSauve / "StatsSpheres",SpheresData, finish);

    SavedAggregates->spheres.lastSaved=SavedAggregates->spheres.size();
}


void ThreadedIO::CreateFile()
{

        if (status[current_thread]==2)
        {
            writer->join();
            delete writer;
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

void ThreadedIO::Write(const fs::path& prefix, shared_ptr<XdmfUnstructuredGrid>& data, bool all)
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

        std::string fileName = std::string(prefix) + filename(NumFile, N);

        if (status[!current_thread] == 2)
        {
            writer->join();
            delete writer;
            status[!current_thread] = 0;
        }
        //std::cout << "Writing " << fileName << std::endl;

        if (false)
        { // Multithread NOT POSSIBLE WITH ALL VERSION OF LIBHDF5
            writer = new std::thread(WriteTask, fileName, &xmfFile[current_thread]);
            status[current_thread] = 2;
        }
        else
        { // Sequential
            WriteTask(fileName, &xmfFile[current_thread]);
            status[current_thread] = 0;
        }
        current_thread = !current_thread;

        NumFile++;

        if (all)
        {
            if (status[!current_thread] == 2)
            {
                writer->join();
                delete writer;
                status[!current_thread] = 0;
            }
            step=0;
        }
    }

}

void WriteTask(const std::string fileName, shared_ptr<XdmfDomain>* data)
{
    shared_ptr<XdmfHDF5Writer> HDF5File = XdmfHDF5Writer::New(fileName+".h5");
    shared_ptr<XdmfWriter> XMFFile = XdmfWriter::New(fileName+".xmf", HDF5File);

    HDF5File->setUseDeflate(false); // do not use compression (too slow)
    HDF5File->setDeflateFactor(0);  // 0 to 6, 6 being the most compressed

    XMFFile->setLightDataLimit(0); // everything go to the hdf5

    // Write data
    (*data)->accept(XMFFile);
}



ThreadedIO::ThreadedIO(PhysicalModel& _physicalmodel, size_t size):
    writer(),
    status({0,0}),
    xmfFile(),
    TimeCollection(),
    physicalmodel(&_physicalmodel),
    current_thread(0),
    step(0),
    NTimePerFile(_physicalmodel.DeltaSauve),
    N(size),
    NumFile(0)
{}


ThreadedIO::~ThreadedIO() noexcept
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

}  // namespace DLCA


#else

namespace DLCA{

__attribute((const)) auto ListSphere::GetData() const
{
    return false;
}
__attribute((const)) auto ListAggregat::GetData() const
{
    return false;
}

void ListSphere::save()
{
    save(false);
}


void ListSphere::save(const bool finish)
{
    UNUSED(finish);
}

void ListAggregat::save()
{
    save(false);
}


void ListAggregat::save(const bool finish)
{
    UNUSED(finish);
}

void StatisticStorage::save()
{
    save(false);
}

void StatisticStorage::save(bool finish)
{
    UNUSED(finish);
}

}  // namespace DLCA


#endif
