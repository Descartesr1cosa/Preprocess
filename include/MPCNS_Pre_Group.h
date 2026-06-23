/*****************************************************************************
@Copyright: NLCFD
@File name: preprocess_group.h
@Author: Descartes
@Version: 1.0
@Date: 2022年10月01日
@Description:	基于1_Preprocess_Data_Structure数据结构与开源的METIS实现图剖分方法，此头文件
                为应用层的功能。
@History:		（修改历史记录列表， 每条修改记录应包括修改日期、修改者及修改内容简述。）
                1、修改了文件名称，与具体问题有关的应用功能头文件名用采用全部小写的方式
*****************************************************************************/
#pragma once
#include "MPCNS_Pre_Data.h"
/*
 * @brief       用来存储网格块分配信息的数据结构
 * @param   block_proc_index    网格在其被分配进程中的序号
 * @param   block_proc_number   网格被分配到的进程号
 * @param   proc_num            进程中的网格块数量
 * @param   pro_grid_number    分配到进程中网格点的数量
 * @param   pro_block_index    分配到进程中网格块的序号
 */
class Preprocess_Group
{
public:
    // 网格在其被分配进程中的序号
    int1D block_proc_index;
    // 网格被分配到的进程号
    int1D block_proc_number;
    // 进程中的网格块数量
    int1D proc_num;
    // 分配到进程中网格点的数量
    int1D pro_grid_number;
    // 分配到进程中网格块的序号
    Array1D<std::vector<int32_t>> pro_block_index;

public:
    Preprocess_Group() {};
    Preprocess_Group(Preprocess_Data_Structure *ptr, Param *par)
    {
        int32_t grid_number = ptr->blk.Getsize1();
        block_proc_index.SetSize(grid_number);
        block_proc_number.SetSize(grid_number);

        int32_t proc_number = par->GetInt("proc_num");
        proc_num.SetSize(proc_number);
        pro_grid_number.SetSize(proc_number);

        pro_block_index.SetSize(proc_number);
    };
    ~Preprocess_Group() = default;
    void Metis_allocate_group(Preprocess_Data_Structure *ptr, Param *par);
    void Load_Balance_Output();
};