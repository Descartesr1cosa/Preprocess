#pragma once
#include <fstream>
#include "core/MPCNS_Pre_Data.h"
#include "partition/MPCNS_Pre_Group.h"

class OUTPUT_MPCNS
{
public:
    Preprocess_Data_Structure *pt;
    Param *pr;
    Preprocess_Group *gp;
    Array2D<Inner> parallel_id;
    bool if_multiphy_facecouple;

public:
    void output(int32_t my_id)
    {
        output_boundary_file(my_id);

        output_inner_file(my_id);

        output_parallel_file(my_id);
        output_gridfile(my_id);
    }

public:
    void output_boundary_file(int32_t my_id);

    /**
     *  @brief 在ACANS输出中，inner文件内id是从零开始，不用正负号表示周期边界条件,MPCNS中id从1开始用正负号表示周期边界条件
     *
     *  @warning 如果存在周期边界条件而使用了ACANS模式输出，虚网格传坐标会出现问题！！！
     */
    void output_inner_file(int32_t my_id);

    /**
     *  @brief 如果是周期边界条件产生的内部边界条件，flag标记为奇数，其他为偶数
     *
     *  @warning 如果存在周期边界条件而使用了ACANS模式输出，虚网格并行传坐标会出现问题！！！
     */
    void process_para_id();

    void output_parallel_file(int32_t my_id);

    /**
     *  @brief 在ACANS模式的输出中，网格文件应严格按照Plot3D输出，这样会损失块的名称
     */
    void output_gridfile(int32_t my_id);

    void output_fvbnd_file();

public:
    OUTPUT_MPCNS(Preprocess_Data_Structure *ptr, Param *par, Preprocess_Group *group)
    {
        pt = ptr;
        pr = par;
        gp = group;
        if_multiphy_facecouple = par->GetBoo("if_multiphy_facecouple");
        if (pr->GetInt("myid") == 0)
            output_fvbnd_file();
        process_para_id();
    };
    ~OUTPUT_MPCNS() = default;
};
