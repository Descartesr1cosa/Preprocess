/*****************************************************************************
@Copyright: NLCFD
@File name: MPCNS_Pre_Datatrans.h
@Author: Descartes
@Version: 1.0
@Date: 2023年08月31日
@Description:	封装形成MPCNS能够统一调用的内部/MPI并行传值的接口，不同的MPI环境只需要更新几个
                基本函数即可。/>对于MPCNS计算使用的数据结构提供了从三维甚至多维数组中抽取形成一维
                以及从一维抽取形成多维数组的接口，方便并行传输数据</
@History:		（修改历史记录列表， 每条修改记录应包括修改日期、修改者及修改内容简述。）
                1、修改了文件名称，基本的功能头文件名用数字1+全部首字母大写的方式
*****************************************************************************/
#pragma once
#include "0_CONST_DEFINE.h"
//=================================================================================================
// MPCNS_IF_MPI=1才会生效
#if MPCNS_IF_MPI == 1
#include <mpi.h>
namespace DATATRANS
{
    void mpi_initial(int arg, char **argv);
    void mpi_finalize();
    void mpi_barrier();
    void mpi_data_send(int tar_id, int send_flag, double *data, int length, MPI_Request *status);
    void mpi_data_recv(int tar_id, int recv_flag, double *data, int length, MPI_Request *status);
    void mpi_rank(int *rank);
    void mpi_size(int *size);
    void mpi_wait(int &count, MPI_Request *request_s, MPI_Status *status_s);
}
#endif
//=================================================================================================