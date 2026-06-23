/*****************************************************************************
@Copyright: NLCFD
@File name: 1_Preprocess_Data_Structure.h
@Author: Descartes
@Version: 1.0
@Date: 2022年10月01日
@Description:	基于map、vector、Array1D等实现的可变数组搭建了用于存储网格块连接、继承关系的数据
                存储方式，基于此数据存储能够通过程序实现网格块的剖分、组合甚至融合的过程，而且不用读
                入所有的网格点信息，防止处理网格量过大而超出内存的限制。根据这些数据结构，要能保证包
                含了正确输出剖分组合后网格、连接关系所需要的全部信息。
@History:		（修改历史记录列表， 每条修改记录应包括修改日期、修改者及修改内容简述。）
                1、修改了文件名称，基本功能的与具体问题有关的头文件名用1_+首字母全部大写的方式
*****************************************************************************/
#pragma once
#include "MPCNS_Pre_Parameter.h"

// Remark 除了ijk的坐标，Tranform的值是以1开始的，其他的数组、数值均从0开始

/**
 * @brief 用来存储网格块内部连接关系的数据结构
 * @param    sub         本块网格坐标下界
 * @param    sup         本块网格坐标上界
 * @param    tar_sub     本块网格坐标对应网格的下界
 * @param    tar_sup     本块网格坐标对应网格的上界
 * @param    volume      本块网格坐标对应网格的面积（长度）
 * @param    is_parallel 本内部关系是否从周期边界条件而来
 *                       若为1,表示周期边界条件，send/recv_flag为负，MPCNS中虚网格坐标不传
 * @remark   这里存储数据均按照inp文件的格式进行存储
 */
class coordinate_inp
{
public:
    int32_t sub[3], sup[3], tar_sub[3], tar_sup[3], volume, is_parallel;

public:
    coordinate_inp() {};
    ~coordinate_inp() = default;
};

/*
 * @brief 用来存储网格块与原始网格块父子关系的数据结构
 * @param    sub         本块网格坐标下界
 * @param    sup         本块网格坐标上界
 * @param    oring_sub   本块网格坐标对应初始网格的下界
 * @param    oring_sup   本块网格坐标对应初始网格的上界
 * @param    oring_num   本块网格坐标对应初始网格的编号
 * @param    Transform   表示本块网格变换到初始网格坐标的调整关系（如本块j为初始网格-i方向则Transform[1]=-1）
 * @param    blk_name    表示本块网格的物理场信息，一般默认为普通流体，可以作为多物理场计算的指示器
 */
class coordinate_blk
{
public:
    int32_t sub[3], sup[3], oring_sub[3], oring_sup[3], Transform[3], oring_num;
    std::string blk_name;

public:
    coordinate_blk() {};
    ~coordinate_blk() = default;
};

/*
 * @brief 用来存储网格块物理边界的数据结构
 * @param    sub     物理边界条件的网格坐标下界
 * @param    sup     物理边界条件的网格坐标上界
 * @param    phy_kind物理边界条件的类型
 */
class coordinate_phy
{
public:
    int32_t sub[3], sup[3], phy_kind;

public:
    coordinate_phy() {};
    ~coordinate_phy() = default;
};

/*
 * @brief 用来存储网格块连接关系的数据结构
 * @param    tar_num  my_num        本内部面是连接哪两个块的
 * @param    inner                  存储连接关系信息
 */
class Inner
{
public:
    int32_t tar_num, my_num; // 本内部面是连接哪两个块的
    std::vector<coordinate_inp> inner;

public:
    Inner() {};
    ~Inner() = default;
};

/*
 * @brief 用来存储网格块自身信息的数据结构
 * @param    hiera         存储剖分合成后父子关系的结构
 * @param    bound          存储物理边界条件的数据结构
 */
class block_info
{
public:
    std::vector<coordinate_blk> hiera; // 存储剖分合成后父子关系的结构
    std::vector<coordinate_phy> bound; // 存储物理边界条件的数据结构
    int32_t ijkmax[3];

public:
    block_info() {};
    ~block_info() = default;
};

/*
 * @brief 用来存储网格所有信息的数据结构，可以用于剖分、合成操作，最后也能根据本数据得到最终网格
 * @param    inp            以二维矩阵的形式存储，用来类似图论中的邻接矩阵的方法表达blk间的连接关系
 * @param    blk            以一维数组的形式存储块自身的信息
 * @param    phy_index      以map存储物理边界条件名称到索引号信息
 * @param    phy_name       以一维数组的形式存储索引号到边界条件名称信息
 */
class Preprocess_Data_Structure
{
public:
    Array2D<Inner> inp;
    Array1D<block_info> blk;

    std::map<std::string, int32_t> phy_index;
    std::vector<std::string> phy_name;

public:
    Preprocess_Data_Structure() {};
    ~Preprocess_Data_Structure() = default;

    /*
     *   @brief 读入inp文件，建立所需要的信息
     */
    void read_inp(Param *par);

    //---------------------------------------------------------------------------------------------
    // 处理周期边界条件所需要的模块

private:
    /*
     *   @brief 如果读入信息中发现了周期边界条件，需要额外读入部分网格点信息，重构出一个内部边界条件
     *          周期边界条件要求命名为一组的PERIOD-'name'和period_'name',对应的name应相同，可以加相同的连接符
     *          例如PERIOD-XUKE与period_XUKE为一组对应的边界条件
     */
    void process_period_condition(Param *par);

    /*
     *   @brief 根据输入块的ijk坐标信息读入坐标，获取对应关系Transform
     */
    void find_corresponding_para(int32_t *Transform, int32_t block_num, int32_t *sub, int32_t *sup, int32_t block_num2, int32_t *tar_sub, int32_t *tar_sup, Param *par);

    /*
     *   @brief 根据输入块的信息读入坐标
     *          为了避免重复多次读取较大的网格文件，需要读取一次后存储，从而能够方便获取a_index到对应坐标的数组xyz &xyz2
     *          对应方法为a_index每一个数映射为0,1（-1为0,1为1），组成的二进制数作为第一个输入e.g.
     *          a_index={-1,1,-1}==>{0 1 0}=2即xyz[2][0~2]
     */
    void read_grid_coordinate(Param *par, int32_t blk_n, int32_t *sub, int32_t *sup, double **xyz);

    /*
     *   @brief 根据周期边界条件的对应点向量平行特性判断对应关系是否正确
     *          begin->tar_begin的矢量与end->tar_end的矢量
     */
    bool check_coordinate_parallel(double *begin, double *tar_begin, double *end, double *tar_end);

    //---------------------------------------------------------------------------------------------

public:
    /*
     *   @brief 检查读入信息，主要检查二维数组inp是否具有对称性（内部关系是互相存在的）
     *          检查是否所有边界面都有边界条件，无论是内部还是物理边界条件
     *          检查是否还存在周期边界条件
     */
    void check_read_info(Param *par);

public:
    /*
     * @brief   网格剖分方法针对此数据结构中的数据，剖分其中blk_number的direction方向的location位置；
     *          blk_number从0开始；direction为1～3表示ijk；location表示位置从1开始
     *
     * @remarks 剖分程序按照如下的思路进行
     *          剖分一个网格块将产生一个新的网格块，我们默认direction方向小于location的仍然用blk_number标号记录，而
     *  direction方向大于location的放在blk的最后位置。这样能够尽可能减少修正连接关系。
     *       针对父子关系即ptr->blk[blk_number]->hiera中所有信息需要重构，其中blk_name，oring_num直接继承，剖分过程
     *  不会改变网格的方向和正负，故而Transform也可以直接继承，需要根据切割的位置和方向对坐标的范围对sub，sup, oring_sub
     *  以及oring_sup进行处理，主要是sub【direction】～sup【direction】范围中出现了location位置，就需要进行处理。例如，
     *  某一个块有多个hiera来源，某一个sub[direction]<location<sup[direction]则sub[direction]<location放在本块，而
     *  location<sup[direction]放在新加块中，其他未被分割的则只需要判断direction方向的坐标范围是在location之上还是之下来
     *  直接修正对应的sub sup坐标即可。
     *       针对物理边界条件ptr->blk[blk_number]->bound也是只需要改变坐标范围，同上；phy_kind直接继承即可
     *       针对最为复杂的内部边界条件需要慎重处理。ptr->inp(blk_number,:),还要新建一个扩充个数 tar_num my_num容易继承
     *  而剩下的inner需要仔细考量。主要有sub sup tar_sub tar_sup。首先根据location判断各个面放在哪一个block，一旦决定后进
     *  行处理即可。对于放在本块的不需要作任何处理，放在最后新加的一个网格块的需要修改direction方向的坐标上下限，同时修改目标块的
     *  对应信息。
     */
    void split_block(Param *par, int32_t blk_number, int32_t direction, int32_t location);
};

/*
 * @brief 根据sub sup的范围查找_inp->inner第几个面为所需要面的面号
 */
int32_t find_index_inner_face(Inner *_inp, int32_t *_sub, int32_t *_sup);