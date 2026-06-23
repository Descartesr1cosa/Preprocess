/*****************************************************************************
@Copyright: NLCFD
@File name: preprocess_split.h
@Author: Descartes
@Version: 1.0
@Date: 2022年10月01日
@Description:	基于1_Preprocess_Data_Structure数据结构的多维贪婪自动剖分算法，此头文件为应用
                层的功能，主要基于Preprocess_Data_Structure数据结构以及其内部提供的剖分功能
                split_block。通过调用该剖分功能，本程序实现自动剖分过程。只需要修改本文件及实现
                就能形成不同的自动剖分算法。
@History:		（修改历史记录列表， 每条修改记录应包括修改日期、修改者及修改内容简述。）
                1、修改了文件名称，与具体问题有关的应用功能头文件名用采用全部小写的方式
*****************************************************************************/
#pragma once
#include "MPCNS_Pre_Data.h"

void Split(Preprocess_Data_Structure *ptr, Param *par);

void split_mid_location(Preprocess_Data_Structure *ptr, Param *par, int32_t split_index);
void split_1D_location(Preprocess_Data_Structure *ptr, Param *par, int32_t split_index, int32_t averge_grid);
void split_2D_location(Preprocess_Data_Structure *ptr, Param *par, int32_t split_index, int32_t averge_grid);
void split_3D_location(Preprocess_Data_Structure *ptr, Param *par, int32_t split_index, int32_t averge_grid);

void output_split_info_inp(Preprocess_Data_Structure *ptr, Param *par);
//=========================================================================================================
