/*****************************************************************************
@Copyright: NLCFD
@File name: preprocess_output.h
@Author: Descartes
@Version: 1.0
@Date: 2022年10月01日
@Description:	基于1_Preprocess_Data_Structure数据结构输出功能实现
@History:		（修改历史记录列表， 每条修改记录应包括修改日期、修改者及修改内容简述。）
                1、修改了文件名称，与具体问题有关的应用功能头文件名用采用全部小写的方式
*****************************************************************************/
#pragma once
#include <fstream>
#include "MPCNS_Pre_Data.h"
#include "MPCNS_Pre_Group.h"

class OUTPUT_ACANS
{
public:
    Preprocess_Data_Structure *pt;
    Param *pr;
    Preprocess_Group *gp;
    Array2D<Inner> parallel_id;

public:
    void output()
    {
        output_inner_file();
        output_boundary_file();
        output_parallel_file();
        output_gridfile();
    }

public:
    /**
     *  @brief 在ACANS模式的输出中，inner文件id是从零开始，不再用正负号表示周期边界条件
     *
     *  @warning 如果存在周期边界条件而使用了ACANS模式输出，虚网格传坐标会出现问题！！！
     */
    void output_inner_file();

    void output_boundary_file();

    /**
     *  @brief 在ACANS模式的输出中，parallel文件flag不再用正负号表示周期边界条件,全部为正
     *
     *  @warning 如果存在周期边界条件而使用了ACANS模式输出，虚网格并行传坐标会出现问题！！！
     */
    void process_para_id();

    void output_parallel_file();

    /**
     *  @brief 在ACANS模式的输出中，网格文件应严格按照Plot3D输出，这样会损失块的名称
     */
    void output_gridfile();

public:
    OUTPUT_ACANS(Preprocess_Data_Structure *ptr, Param *par, Preprocess_Group *group)
    {
        pt = ptr;
        pr = par;
        gp = group;
    };
    ~OUTPUT_ACANS() = default;
};