#include "workflow/MPCNS_Pre_App.h"

#include <cmath>
#include <iostream>
#include <memory>

#include "workflow/MPCNS_Pre_CubicGridGenerator.h"
#include "workflow/MPCNS_Pre_CubicSphereGridGenerator.h"
#include "core/MPCNS_Pre_Data.h"
#include "partition/MPCNS_Pre_Datatrans.h"
#include "partition/MPCNS_Pre_DecOrientationTest.h"
#include "partition/MPCNS_Pre_Group.h"
#include "workflow/MPCNS_Pre_InputRefinement.h"
#include "io/MPCNS_Pre_IOinfo.h"
#include "core/MPCNS_Pre_Parameter.h"
#include "partition/MPCNS_Pre_Split.h"
#include "report/MPCNS_Pre_TopologyReport.h"
#include "io/Pre_ACANSOutput.h"
#include "io/Pre_MPCNSOutput.h"

#if MPCNS_Para_COUT == 1
std::mutex g_Cs;
#endif

namespace
{
class MpiSession
{
public:
    MpiSession(int argc, char **argv, int &rank)
    {
        rank = 0;
#if MPCNS_IF_MPI == 1
        DATATRANS::mpi_initial(argc, argv);
        DATATRANS::mpi_rank(&rank);
#else
        (void)argc;
        (void)argv;
#endif
    }

    ~MpiSession()
    {
#if MPCNS_IF_MPI == 1
        DATATRANS::mpi_barrier();
        DATATRANS::mpi_finalize();
#endif
    }

    void Barrier() const
    {
#if MPCNS_IF_MPI == 1
        DATATRANS::mpi_barrier();
#endif
    }

    int Size() const
    {
#if MPCNS_IF_MPI == 1
        int process_number = 0;
        DATATRANS::mpi_size(&process_number);
        return process_number;
#else
        return 1;
#endif
    }
};

class SplitGroupBuilder
{
public:
    SplitGroupBuilder(int rank, Param &param) : rank_(rank), param_(param) {}

    void BuildIfNeeded() const
    {
        if (rank_ != 0 || param_.GetBoo("if_split_group_info"))
            return;

        PrintInputSummary();

        std::cout << "---->(2) Reading the inp file...\n";
        std::unique_ptr<Preprocess_Data_Structure> data(new Preprocess_Data_Structure);
        data->read_inp(&param_);
        std::cout << "\t-->Grid and inp files have been processed ! ! !\n";

        std::cout << "---->(3) Splitting the large grid blocks...\n";
        Split(data.get(), &param_);
        std::cout << "\t-->The blocks have been split ! ! !\n";

        DecOrientationTest::ApplyIfEnabled(data.get(), &param_);

        std::cout << "---->(4) Settling down the grid blocks...\n";
        Preprocess_Group group(data.get(), &param_);
        group.Metis_allocate_group(data.get(), &param_);
        group.Load_Balance_Output();
        std::cout << "\t-->The blocks have been grouped ! ! !\n";

        std::cout << "---->(5) Output the grid and link information...\n";
        output_split_group_info(data.get(), &group);
        TopologyReport::Generate(data.get(), &group, &param_);

        std::cout << "\t-->The grid and link information have been output ^_^ \n";
    }

private:
    void PrintInputSummary() const
    {
        std::cout << "\t\tThe grids file is #" << param_.GetStr("gfilename") << "#\n";
        std::cout << "\t\tThe inp file is #" << param_.GetStr("bfilename") << "#\n";
        std::cout << "\t\tThe fvbnd file is #" << param_.GetStr("ffilename") << "#\n";
        std::cout << "\t-->The input setup file has been read successfully ! ! !\n";
    }

    int rank_;
    Param &param_;
};

class GeometryInfo
{
public:
    void Load()
    {
        data_.reset(new Preprocess_Data_Structure);
        group_.reset(new Preprocess_Group);
        input_split_group_info(data_.get(), group_.get());
    }

    Preprocess_Data_Structure *Data() const { return data_.get(); }
    Preprocess_Group *Group() const { return group_.get(); }

private:
    std::unique_ptr<Preprocess_Data_Structure> data_;
    std::unique_ptr<Preprocess_Group> group_;
};

class PeriodicBoundaryGuard
{
public:
    static void WarnForAcansIfNeeded(Preprocess_Data_Structure &data)
    {
        for (int i = 0; i < data.phy_name.size(); i++)
        {
            if (data.phy_name[i].size() > 7)
            {
                std::string str = data.phy_name[i].substr(0, 7);
                if (str == "PERIOD-" || str == "period_")
                {
                    std::cout << std::endl;
                    std::cout << "#########################################################################\n";
                    std::cout << "\t\t\t--->Warning ! ! !<---\n";
                    std::cout << "PERIOD- and period_ boundary condition exist ! ! !\n";
                    std::cout << "Errors will appear during the simulation process Later !\n";
                    std::cout << "#########################################################################\n\n";
                    std::cout << "-->if_preprocess_grid = 1 \t should be used...\n";
                    std::cout << "Press any key to CONTINUE the output process, but I do not recommend to do so\n";
                    std::cin >> str;
                    break;
                }
            }
        }
    }
};

class GeometryOutputWorkflow
{
public:
    GeometryOutputWorkflow(int rank, Param &param, const MpiSession &mpi)
        : rank_(rank), param_(param), mpi_(mpi) {}

    void Run() const
    {
        const int32_t preprocess_method = param_.GetInt("preprocess_method");
        if (preprocess_method == 0)
        {
            OutputAcansGeometry();
        }
        else if (preprocess_method == 1)
        {
            OutputMpcnsGeometry();
        }
    }

private:
    void OutputAcansGeometry() const
    {
#if MPCNS_IF_MPI == 0
        OutputAcansSerial();
#else
        WarnIfAcansLaunchedWithMultipleProcesses();
        if (rank_ == 0)
            OutputAcansSerial();
#endif
    }

    void OutputAcansSerial() const
    {
        std::cout << "---->(5) Reading the binary split_group file...\n";
        GeometryInfo geometry;
        geometry.Load();
        std::cout << "\t-->Successfully reading the binary split_group file.\n";

        PeriodicBoundaryGuard::WarnForAcansIfNeeded(*geometry.Data());

        std::cout << "---->(6) Starting to ouput the ACANS geometry txt information...\n";
        OUTPUT_ACANS output(geometry.Data(), &param_, geometry.Group());
        output.output();
        std::cout << "\t-->The ACANS geometry txt files have been ouput...\n";
    }

    void OutputMpcnsGeometry() const
    {
#if MPCNS_IF_MPI == 0
        OutputMpcnsSerial();
#else
        WarnIfTooManyProcesses();
        OutputMpcnsParallel();
#endif
    }

    void OutputMpcnsSerial() const
    {
        std::cout << "---->(5) Reading the binary split_group file...\n";
        GeometryInfo geometry;
        geometry.Load();
        std::cout << "\t-->Successfully reading the binary split_group file.\n";

        std::cout << "---->(6) Starting to ouput the MPCNS geometry txt information...\n";
        OUTPUT_MPCNS output(geometry.Data(), &param_, geometry.Group());
        for (int32_t output_rank = 0; output_rank < param_.GetInt("proc_num"); output_rank++)
            output.output(output_rank);

        std::cout << "\t-->The MPCNS geometry txt files have been ouput...\n";
    }

#if MPCNS_IF_MPI == 1
    void OutputMpcnsParallel() const
    {
        if (rank_ == 0)
            std::cout << "---->(5) Reading the binary split_group file...\n";

        GeometryInfo geometry;
        geometry.Load();
        mpi_.Barrier();

        if (rank_ == 0)
            std::cout << "\t-->Successfully reading the binary split_group file!\n";

        if (rank_ == 0)
            std::cout << "---->(6) Starting to ouput the MPCNS geometry txt information...\n";

        OUTPUT_MPCNS output(geometry.Data(), &param_, geometry.Group());
        const int process_number = mpi_.Size();
        const int32_t proc_num = param_.GetInt("proc_num");
        if (process_number > proc_num)
        {
            if (rank_ < proc_num)
                output.output(rank_);
        }
        else
        {
            for (int32_t output_rank = 0; output_rank < param_.GetInt("proc_num"); output_rank++)
            {
                if (std::fmod(output_rank + process_number - rank_, process_number) == 0)
                    output.output(output_rank);
            }
        }

        mpi_.Barrier();
        if (rank_ == 0)
            std::cout << "\t-->The MPCNS geometry txt files have been ouput...\n";
    }

    void WarnIfAcansLaunchedWithMultipleProcesses() const
    {
        if (mpi_.Size() != 1 && rank_ == 0)
        {
            std::cout << "#########################################################################\n";
            std::cout << "Serial program should be launched, but no influence will be on the output\n";
            std::cout << "#########################################################################\n";
        }
    }

    void WarnIfTooManyProcesses() const
    {
        if (mpi_.Size() > param_.GetInt("proc_num") && rank_ == 0)
        {
            std::cout << "###############################################################################\n";
            std::cout << "Too Many Parallel programs are launched, but no influence will be on the output\n";
            std::cout << "###############################################################################\n";
        }
    }
#endif

    int rank_;
    Param &param_;
    const MpiSession &mpi_;
};
} // namespace

int PreprocessApplication::Run(int argc, char **argv)
{
    int rank = 0;
    MpiSession mpi(argc, argv, rank);

    Param param;
    param.ReadParam(rank);

    if (rank == 0)
        CubicSphereGridGenerator::GenerateIfEnabled(param);
    if (rank == 0)
        CubicPeriodicGridGenerator::GenerateIfEnabled(param);
    mpi.Barrier();

    InputRefinement::RefineIfEnabled(param);
    mpi.Barrier();

    SplitGroupBuilder(rank, param).BuildIfNeeded();
    mpi.Barrier();

    GeometryOutputWorkflow(rank, param, mpi).Run();

    if (rank == 0)
        std::cout << "^_^ You can press any key to continue ^_^ \n";

    return 0;
}
