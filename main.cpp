#include <iostream>

#include "MPCNS_Pre_Datatrans.h"
#include "MPCNS_Pre_Parameter.h"
#include "MPCNS_Pre_Data.h"
#include "MPCNS_Pre_Split.h"
#include "MPCNS_Pre_Group.h"
#include "MPCNS_Pre_IOinfo.h"
#include "Pre_ACANSOutput.h"
#include "Pre_MPCNSOutput.h"
#include <math.h>

// #include "0_ARRAY_DEFINE.h"
// #include "1_Preprocess_Parameter.h"
// #include "1_Preprocess_Data_Structure.h"
// #include "preprocess_split.h"
// #include "preprocess_output.h"
// #include "preprocess_group.h"
// #include "preprocess_grid.h"

#if MPCNS_Para_COUT == 1
CRITICAL_SECTION g_Cs;
#endif

int main(int arg, char **argv)
{
    //=================================================================================================
    // MAIN函数通用开头
    int myid = 0;
//=================================================================================================

//=================================================================================================
// MPCNS_IF_MPI=1才会生效，mpi初始化
#if MPCNS_IF_MPI == 1
#if MPCNS_Para_COUT == 1
    InitializeCriticalSection(&g_Cs);
#endif
    DATATRANS::mpi_initial(arg, argv);
    DATATRANS::mpi_rank(&myid);
#endif
    //=================================================================================================

    //=================================================================================================
    // 所有进程都读入参数，存入Param的数据结构中
    Param *par = new Param;
    par->ReadParam(myid);
    //=================================================================================================

    //=================================================================================================
    // 首先采用串行的方式进行网格的剖分组合工作，输出包含网格剖分组合结果的二进制信息文件
    // 根据这些二进制信息文件可以从原始输入文件中获取全部几何信息，便于后面的处理

    if (myid == 0 && !par->GetBoo("if_split_group_info"))
    {
        std::cout << "\t\tThe grids file is #" << par->GetStr("gfilename") << "#\n";
        std::cout << "\t\tThe inp file is #" << par->GetStr("bfilename") << "#\n";
        std::cout << "\t\tThe fvbnd file is #" << par->GetStr("ffilename") << "#\n";
        std::cout << "\t-->The input setup file has been read successfully ! ! !\n";
        //-------------------------------------------------------------------------------
        std::cout << "---->(2) Reading the inp file...\n";
        Preprocess_Data_Structure *ptr = new Preprocess_Data_Structure;
        ptr->read_inp(par);
        std::cout << "\t-->Grid and inp files have been processed ! ! !\n";
        //-------------------------------------------------------------------------------
        std::cout << "---->(3) Splitting the large grid blocks...\n";
        Split(ptr, par);
        std::cout << "\t-->The blocks have been split ! ! !\n";
        //-------------------------------------------------------------------------------
        std::cout << "---->(4) Settling down the grid blocks...\n";
        Preprocess_Group grp(ptr, par);
        grp.Metis_allocate_group(ptr, par);
        grp.Load_Balance_Output();
        std::cout << "\t-->The blocks have been grouped ! ! !\n";
        //-------------------------------------------------------------------------------
        std::cout << "---->(5) Output the grid and link information...\n";
        // 将剖分组合信息输出，从而能够并行读入，分别生成各个进程自己的几何信息
        output_split_group_info(ptr, &grp);
        //---------------------------------------------------------------------

        std::cout << "\t-->The grid and link information have been output ^_^ \n";
        delete ptr;
    }
//=================================================================================================

//=================================================================================================
// 将所有进程拦到此处
// MPCNS_IF_MPI=1才会生效，mpi同步
#if MPCNS_IF_MPI == 1
    DATATRANS::mpi_barrier();
#endif
    //=================================================================================================

    //=================================================================================================
    int32_t preprocess_method = par->GetInt("preprocess_method");
    if (preprocess_method == 0)
    {
        //---------------------------------------------------------------------
        // 由于输出的时候需要对每个块网格几乎所有点进行读入并输出,本程序针对ACANS采用串行，效率略低
        // 输出ACANS文本读入grd、boundary_condition的信息
        // 存在【周期边界条件】请勿使用！ ！ ！

#if MPCNS_IF_MPI == 0
        //---------------------------------------------------------------------
        // 将剖分组合信息读入，分别生成各个进程自己的几何信息
        std::cout << "---->(5) Reading the binary split_group file...\n";
        Preprocess_Group *grp = new Preprocess_Group;
        Preprocess_Data_Structure *ptr = new (Preprocess_Data_Structure);
        input_split_group_info(ptr, grp);
        std::cout << "\t-->Successfully reading the binary split_group file.\n";
        //---------------------------------------------------------------------
        // 为了防止错误，这里检测周期边界条件
        for (int i = 0; i < ptr->phy_name.size(); i++)
        {
            if (ptr->phy_name[i].size() > 7)
            {
                std::string str = ptr->phy_name[i].substr(0, 7);
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
        //---------------------------------------------------------------------
        // 开始生成各个进程的网格并输出到文件，形成原始的grd、boundary_condition文本文件
        std::cout << "---->(6) Starting to ouput the ACANS geometry txt information...\n";
        OUTPUT_ACANS *opt = new OUTPUT_ACANS(ptr, par, grp);
        opt->output();
        delete opt, ptr, grp;
        std::cout << "\t-->The ACANS geometry txt files have been ouput...\n";
#else
        int32_t process_number = 0;
        DATATRANS::mpi_size(&process_number);
        if (process_number != 1 && myid == 0)
        {
            std::cout << "#########################################################################\n";
            std::cout << "Serial program should be launched, but no influence will be on the output\n";
            std::cout << "#########################################################################\n";
        }

        if (myid == 0)
        {
            //---------------------------------------------------------------------
            // 将剖分组合信息读入，分别生成各个进程自己的几何信息
            std::cout << "---->(5) Reading the binary split_group file...\n";
            Preprocess_Group *grp = new Preprocess_Group;
            Preprocess_Data_Structure *ptr = new (Preprocess_Data_Structure);
            input_split_group_info(ptr, grp);
            std::cout << "\t-->Successfully reading the binary split_group file.\n";
            //---------------------------------------------------------------------
            // 为了防止错误，这里检测周期边界条件
            for (int i = 0; i < ptr->phy_name.size(); i++)
            {
                if (ptr->phy_name[i].size() > 7)
                {
                    std::string str = ptr->phy_name[i].substr(0, 7);
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
            //---------------------------------------------------------------------
            // 开始生成各个进程的网格并输出到文件，形成原始的grd、boundary_condition文本文件
            std::cout << "---->(6) Starting to ouput the ACANS geometry txt information...\n";
            OUTPUT_ACANS *opt = new OUTPUT_ACANS(ptr, par, grp);
            opt->output();
            delete opt, ptr, grp;
            std::cout << "\t-->The ACANS geometry txt files have been ouput...\n";
        }
#endif
    }
    else if (preprocess_method == 1)
    {
        // 存在【周期边界条件】使用preprocess_method=1输出MPCNS使用的文件（仍然会输出grd,boundary_condition）
        // 注意在MPCNS中需要正确处理周期边界条件（周期边界条件直接加入inp信息中，而虚网格不能传！！）
        // 存在【周期边界条件】请勿使用preprocess_method=0，否则虚网格坐标传递会出现错误
        // 输出MPCNS文本读入grd、boundary_condition的信息,是MPCNS形式！！确保求解程序能够正确读取
        //  MPCNS形式与ACANS的区别：
        /**
         * (1)修改了grd文件格式，在网格块的大小之前，总网格块块数之后有各个块的物理场名称（Fluid Solid或其他默认名称）
         * (2)修改了parallel   X.txt中的send_flag recv_flag信息，对于由周期边界条件产生的一律为负，这一点不影响ACANS计算
         * (3)修改了inner   X.txt中的id信息，块号不再从0开始，均从1开始编号，此外由周期边界条件产生的一律为负
         */
#if MPCNS_IF_MPI == 0
        //---------------------------------------------------------------------
        // 将剖分组合信息读入，生成几何信息
        std::cout << "---->(5) Reading the binary split_group file...\n";
        Preprocess_Group *grp = new Preprocess_Group;
        Preprocess_Data_Structure *ptr = new (Preprocess_Data_Structure);
        input_split_group_info(ptr, grp);
        std::cout << "\t-->Successfully reading the binary split_group file.\n";
        //---------------------------------------------------------------------
        // 开始生成各个进程的网格并输出到文件，形成原始的grd、boundary_condition文本文件
        std::cout << "---->(6) Starting to ouput the MPCNS geometry txt information...\n";
        OUTPUT_MPCNS *opt = new OUTPUT_MPCNS(ptr, par, grp);

        for (int32_t myid_num = 0; myid_num < par->GetInt("proc_num"); myid_num++)
            opt->output(myid_num);

        delete opt, ptr, grp;
        std::cout << "\t-->The MPCNS geometry txt files have been ouput...\n";
#else
        int32_t process_number = 0;
        DATATRANS::mpi_size(&process_number);
        int32_t proc_num = par->GetInt("proc_num");
        if (process_number > proc_num && myid == 0)
        {
            std::cout << "###############################################################################\n";
            std::cout << "Too Many Parallel programs are launched, but no influence will be on the output\n";
            std::cout << "###############################################################################\n";
        }
        //---------------------------------------------------------------------
        // 将剖分组合信息读入，分别生成各个进程自己的几何信息
        if (myid == 0)
            std::cout << "---->(5) Reading the binary split_group file...\n";
        Preprocess_Group *grp = new Preprocess_Group;
        Preprocess_Data_Structure *ptr = new (Preprocess_Data_Structure);
        input_split_group_info(ptr, grp);
        DATATRANS::mpi_barrier();
        if (myid == 0)
            std::cout << "\t-->Successfully reading the binary split_group file!\n";
        //---------------------------------------------------------------------
        // 开始【并行】生成各个进程的网格并输出到文件，形成原始的grd、boundary_condition文本文件
        if (myid == 0)
            std::cout << "---->(6) Starting to ouput the MPCNS geometry txt information...\n";
        OUTPUT_MPCNS *opt = new OUTPUT_MPCNS(ptr, par, grp);

        if (process_number > proc_num)
        {
            if (myid < proc_num)
                opt->output(myid);
        }
        else
        {
            for (int32_t myid_num = 0; myid_num < par->GetInt("proc_num"); myid_num++)
            {
                if (fmod(myid_num + process_number - myid, process_number) == 0)
                    opt->output(myid_num);
            }
        }

        delete ptr, grp, opt;
        DATATRANS::mpi_barrier();
        if (myid == 0)
            std::cout << "\t-->The MPCNS geometry txt files have been ouput...\n";
            //---------------------------------------------------------------------
#endif
    }
    //=================================================================================================

    // if (par->preprocess_method == 2)
    // {
    //     int32_t process_number = 0;
    //     DATATRANS::mpi_size(&process_number);
    //     if (process_number != par->proc_num)
    //     {
    //         std::cout << "The target process number is\t" << par->proc_num << ", but only " << process_number << " is launched.\n";
    //         exit(-1);
    //     }

    //     //---------------------------------------------------------------------
    //     //将剖分组合信息读入，分别生成各个进程自己的几何信息
    //     if (myid == 0)
    //         std::cout << "---->(5) Reading the binary split_group file...\n";
    //     preprocess_group *grp = new (preprocess_group);
    //     Preprocess_Data_Structure *ptr = new (Preprocess_Data_Structure);
    //     input_split_group_info(ptr, grp);
    //     DATATRANS::mpi_barrier();
    //     if (myid == 0)
    //         std::cout << "\t-->Successfully reading the binary split_group file!\n";
    //     //---------------------------------------------------------------------
    //     //开始【并行】生成各个进程的网格并输出到文件，形成原始的grd、boundary_condition文本文件
    //     if (myid == 0)
    //         std::cout << "---->(6) Starting to ouput the ACANS/MPCNS geometry txt information...\n";
    //     preprocess_output *opt = new OUTPUT_MPCNS;
    //     opt->build_parallel(ptr, par, grp, myid);
    //     opt->output_parallel();
    //     delete ptr, grp, opt;
    //     DATATRANS::mpi_barrier();
    //     if (myid == 0)
    //         std::cout << "\t-->The ACANS/MPCNS geometry txt files have been ouput...\n";
    //     //---------------------------------------------------------------------
    //     //各个进程开始读入几何信息，开展预处理过程
    //     preprocess_grid *grd = new (preprocess_grid);
    //     if (myid == 0)
    //         std::cout << "---->(7) Starting the geometry process...\n";
    //     grd->process_grid_from_ordinary_geometry_File(myid, par->ngg, par->scale);
    //     DATATRANS::mpi_barrier();
    //     delete grd;
    //     if (myid == 0)
    //         std::cout << "\t-->Successfully Preprocessing the Geometry Files!\n";
    // }
    // else if (par->preprocess_method == 1)
    // {
    //     //---------------------------------------------------------------------
    //     //由于输出的时候需要对每个块网格几乎所有点进行读入并输出，采用串行输出，速度较慢
    //     //输出MPCNS文本读入grd、boundary_condition的信息,是MPCNS形式！！确保求解程序能够正确读取
    //     // MPCNS形式与ACANS的区别：
    //     /**
    //      * (1)修改了grd文件格式，在网格块的大小之前，总网格块块数之后有各个块的物理场名称（Fluid Solid或其他默认名称）
    //      * (2)修改了parallel   X.txt中的send_flag recv_flag信息，对于由周期边界条件产生的一律为负，这一点不影响ACANS计算
    //      * (3)修改了inner   X.txt中的id信息，块号不再从0开始，均从1开始编号，此外由周期边界条件产生的一律为负
    //      */
    //     int32_t process_number = 0;
    //     DATATRANS::mpi_size(&process_number);
    //     if (process_number != 1 && myid == 0)
    //     {
    //         std::cout << "#########################################################################\n";
    //         std::cout << "Serial program should be launched, but no influence will be on the output\n";
    //         std::cout << "#########################################################################\n";
    //     }

    //     if (myid == 0)
    //     {
    //         //---------------------------------------------------------------------
    //         //将剖分组合信息读入，生成几何信息
    //         std::cout << "---->(5) Reading the binary split_group file...\n";
    //         preprocess_group *grp = new (preprocess_group);
    //         Preprocess_Data_Structure *ptr = new (Preprocess_Data_Structure);
    //         input_split_group_info(ptr, grp);
    //         std::cout << "\t-->Successfully reading the binary split_group file.\n";
    //         //---------------------------------------------------------------------
    //         //开始生成各个进程的网格并输出到文件，形成原始的grd、boundary_condition文本文件
    //         std::cout << "---->(6) Starting to ouput the ACANS/MPCNS geometry txt information...\n";
    //         preprocess_output *opt = new OUTPUT_MPCNS;
    //         opt->build(ptr, par, grp);
    //         opt->output();
    //         delete opt, ptr, grp;
    //         std::cout << "\t-->The ACANS/MPCNS geometry txt files have been ouput...\n";
    //     }
    // }
    // else if (par->preprocess_method == 0)
    // {
    //     //---------------------------------------------------------------------
    //     //由于输出的时候需要对每个块网格几乎所有点进行读入并输出，采用串行输出，速度较慢
    //     //输出ACANS文本读入grd、boundary_condition的信息
    //     //存在【周期边界条件】请勿使用！ ！ ！
    //     int32_t process_number = 0;
    //     DATATRANS::mpi_size(&process_number);
    //     if (process_number != 1 && myid == 0)
    //     {
    //         std::cout << "#########################################################################\n";
    //         std::cout << "Serial program should be launched, but no influence will be on the output\n";
    //         std::cout << "#########################################################################\n";
    //     }

    //     if (myid == 0)
    //     {
    //         //---------------------------------------------------------------------
    //         //将剖分组合信息读入，分别生成各个进程自己的几何信息
    //         std::cout << "---->(5) Reading the binary split_group file...\n";
    //         preprocess_group *grp = new (preprocess_group);
    //         Preprocess_Data_Structure *ptr = new (Preprocess_Data_Structure);
    //         input_split_group_info(ptr, grp);
    //         std::cout << "\t-->Successfully reading the binary split_group file.\n";
    //         //---------------------------------------------------------------------
    //         //为了防止错误，这里检测周期边界条件
    //         for (int i = 0; i < ptr->phy_name.size(); i++)
    //         {
    //             if (ptr->phy_name[i].size() > 7)
    //             {
    //                 std::string str = ptr->phy_name[i].substr(0, 7);
    //                 if (str == "PERIOD-" || str == "period_")
    //                 {
    //                     std::cout << std::endl;
    //                     std::cout << "#########################################################################\n";
    //                     std::cout << "\t\t\t--->Warning ! ! !<---\n";
    //                     std::cout << "PERIOD- and period_ boundary condition exist ! ! !\n";
    //                     std::cout << "Errors will appear during the simulation process Later !\n";
    //                     std::cout << "#########################################################################\n\n";
    //                     std::cout << "-->if_preprocess_grid = 1 \t should be used...\n";
    //                     std::cout << "Press any key to CONTINUE the output process, but I do not recommend to do so\n";
    //                     std::cin >> str;
    //                     break;
    //                 }
    //             }
    //         }
    //         //---------------------------------------------------------------------
    //         //开始生成各个进程的网格并输出到文件，形成原始的grd、boundary_condition文本文件
    //         std::cout << "---->(6) Starting to ouput the ACANS/MPCNS geometry txt information...\n";
    //         preprocess_output *opt = new OUTPUT_ACANS;
    //         opt->build(ptr, par, grp);
    //         opt->output();
    //         delete opt, ptr, grp;
    //         std::cout << "\t-->The ACANS/MPCNS geometry txt files have been ouput...\n";
    //     }
    // }

    delete par;
    if (myid == 0)
        std::cout << "^_^ You can press any key to continue ^_^ \n";

//=================================================================================================
// MPCNS_IF_MPI=1才会生效，mpi初始化
#if MPCNS_IF_MPI == 1
    DATATRANS::mpi_barrier();
    DATATRANS::mpi_finalize();
#endif
    //=================================================================================================

    return 0;
}
