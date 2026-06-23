/*****************************************************************************
@Copyright: NLCFD
@File name: MPCNS_Pre_IOinfo.h
@Author: Descartes
@Version: 1.0
@Date: 2023年8月31日
@Description:	基于1_Preprocess_Data_Structure数据结构输出功能实现
@History:		（修改历史记录列表， 每条修改记录应包括修改日期、修改者及修改内容简述。）
                1、修改了文件名称，与具体问题有关的应用功能头文件名用采用全部小写的方式
*****************************************************************************/

#pragma once
#include <fstream>
#include "MPCNS_Pre_Data.h"
#include "MPCNS_Pre_Group.h"

/**
 *  @brief      通过将剖分以及分配的信息输出到二进制文件，从而实现可以根据剖分信息直接并行输出网格
 *              也可以根据现有信息继续开展其他操作（在有较大网格的时候，处理周期边界条件非常耗时）
 */
void output_split_group_info(Preprocess_Data_Structure *ptr, Preprocess_Group *group);

/**
 *  @brief      通过将剖分以及分配的信息输出得到的二进制文件，从而实现剖分信息的重建，从而能够
 *              根据现有信息继续开展其他操作（在有较大网格的时候，处理周期边界条件非常耗时）
 */
void input_split_group_info(Preprocess_Data_Structure *ptr, Preprocess_Group *group);

void output_split_group_info_DEBUG_only(Preprocess_Data_Structure *ptr, Preprocess_Group *group);