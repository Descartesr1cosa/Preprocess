/*****************************************************************************
@Copyright: NLCFD
@File name: MPCNSS_Pre_Data.cpp
@Author: Descartes
@Version: 1.0
@Date: 2025年1月8日
@Description:	实现inp读入处理程序
@History:		（修改历史记录列表， 每条修改记录应包括修改日期、修改者及修改内容简述。）
                1、修改了621~623行x,y,z Setsize()的第三个参数应为mz而不是my, 713~715同
*****************************************************************************/

#include "MPCNS_Pre_Data.h"
#include <algorithm>
#include <cctype>
#include <math.h>

namespace
{
const int32_t CANONICAL_POLE_KIND = 7;

bool IsPoleAliasKind(int32_t phy_kind)
{
    return phy_kind == 72 || phy_kind == 73;
}

void NormalizeInpPoleBoundaryKind(Preprocess_Data_Structure *ptr)
{
    int32_t changed_face_number = 0;
    for (int32_t iblock = 0; iblock < ptr->blk.Getsize1(); iblock++)
    {
        for (int32_t iface = 0; iface < ptr->blk(iblock).bound.size(); iface++)
        {
            coordinate_phy &boundary = ptr->blk(iblock).bound[iface];
            if (IsPoleAliasKind(boundary.phy_kind))
            {
                boundary.phy_kind = CANONICAL_POLE_KIND;
                changed_face_number++;
            }
        }
    }

    if (changed_face_number > 0)
    {
        std::cout << "\t-->Normalize " << changed_face_number
                  << " Pole boundary faces in inp from kind 72/73 to "
                  << CANONICAL_POLE_KIND << ".\n";
    }
}

std::string ToLower(std::string str)
{
    std::transform(str.begin(), str.end(), str.begin(), [](unsigned char c) { return std::tolower(c); });
    return str;
}

std::vector<std::string> SplitKeywords(const std::string &keyword_group)
{
    std::vector<std::string> keywords;
    std::string keyword;
    for (int32_t i = 0; i < keyword_group.size(); i++)
    {
        if (keyword_group[i] == '|')
        {
            if (!keyword.empty())
                keywords.push_back(keyword);
            keyword.clear();
        }
        else
        {
            keyword.push_back(keyword_group[i]);
        }
    }
    if (!keyword.empty())
        keywords.push_back(keyword);
    return keywords;
}

std::string NormalizeBlockName(const std::string &raw_block_name, Param *par)
{
    if (!par->HasBoo("if_multiphy_facecouple") || !par->GetBoo("if_multiphy_facecouple"))
        return raw_block_name;

    std::string matched_name;
    const std::string raw_lower = ToLower(raw_block_name);
    if (par->HasStr_List("canonical_block_name_keywords"))
    {
        List<std::string> keyword_groups = par->GetStr_List("canonical_block_name_keywords");
        for (std::map<std::string, std::string>::iterator iter = keyword_groups.data.begin(); iter != keyword_groups.data.end(); iter++)
        {
            const std::string canonical_name = iter->first;
            const std::vector<std::string> keywords = SplitKeywords(iter->second);
            bool group_matched = false;
            for (int32_t i = 0; i < keywords.size(); i++)
            {
                if (raw_lower.find(ToLower(keywords[i])) != std::string::npos)
                {
                    group_matched = true;
                    break;
                }
            }
            if (!group_matched)
                continue;

            if (!matched_name.empty() && matched_name != canonical_name)
            {
                std::cout << "#Fatal Error: Block name '" << raw_block_name
                          << "' matches more than one canonical block-name keyword group: "
                          << matched_name << " and " << canonical_name << std::endl;
                exit(-1);
            }
            matched_name = canonical_name;
        }
    }

    if (matched_name.empty())
        return "Fluid";
    return matched_name;
}
}
//=================================================================================================
//---------------------------------------读入inp文件信息---------------------------------------------

void Preprocess_Data_Structure::read_inp(Param *par)
{
    int32_t blk_num, dimension = par->GetInt("dimension");
    std::string line, filename = par->GetStr("bfilename");
    std::fstream file;
    file.open(filename, std::ios_base::in);
    std::getline(file, line);
    // 读入网格块的个数
    file >> blk_num;

    inp.SetSize(blk_num, blk_num);
    blk.SetSize(blk_num);

    for (int32_t i = 0; i < blk_num; i++)
        for (int32_t j = 0; j < blk_num; j++)
        {
            inp(i, j).my_num = i;
            inp(i, j).tar_num = j;
        }

    for (int32_t iblock = 0; iblock < blk_num; iblock++)
    {
        coordinate_blk temp;
        for (int32_t i = 0; i < 3; i++)
            temp.sub[i] = 1;
        if (dimension == 2)
        {
            file >> temp.sup[0] >> temp.sup[1];
            temp.sup[2] = 1;
        }
        else
        {
            file >> temp.sup[0] >> temp.sup[1] >> temp.sup[2];
        }

        for (int32_t i = 0; i < 3; i++)
        {
            temp.oring_sub[i] = temp.sub[i];
            temp.oring_sup[i] = temp.sup[i];
            temp.Transform[i] = i + 1;
        }
        temp.oring_num = iblock;
        file >> temp.blk_name;
        temp.blk_name = NormalizeBlockName(temp.blk_name, par);
        // 完成了coordinate_blk读入
        for (int32_t i = 0; i < 3; i++)
            blk(iblock).ijkmax[i] = temp.sup[i];

        blk(iblock).hiera.push_back(temp);
        //=====================================================

        int32_t face_num;
        file >> face_num;

        int sub[3], sup[3], tar_sub[3], tar_sup[3], index;
        for (int32_t i = 0; i < face_num; i++)
        {
            for (int j = 0; j < dimension; j++)
                file >> sub[j] >> sup[j];
            file >> index;
            if (dimension == 2)
            {
                sub[2] = 1;
                sup[2] = 1;
            }
            // file >> sub[0] >> sup[0] >> sub[1] >> sup[1] >> sub[2] >> sup[2] >> index;
            if (index > 0)
            {
                /* Physical boundary conditions */
                coordinate_phy temp_phy;
                temp_phy.phy_kind = index;
                for (int32_t i = 0; i < 3; i++)
                {
                    temp_phy.sub[i] = sub[i];
                    temp_phy.sup[i] = sup[i];
                }
                // 完成了一个coordinate_phy读入
                blk(iblock).bound.push_back(temp_phy);
            }
            else
            {
                /* Inner boundary conditions */
                coordinate_inp temp_inp;
                for (int j = 0; j < dimension; j++)
                    file >> tar_sub[j] >> tar_sup[j];
                file >> index;
                if (dimension == 2)
                {
                    tar_sub[2] = 1;
                    tar_sup[2] = 1;
                }
                // file >> tar_sub[0] >> tar_sup[0] >> tar_sub[1] >> tar_sup[1] >> tar_sub[2] >> tar_sup[2] >> index;
                if (index > 0)
                {
                    // 从inp文件读入的信息，真实的内部边界条件，虚网格和求解时边界数据均需传递
                    // 一般商用软件的网格Plot3D格式输出就是此形式，直接表示块号
                    temp_inp.is_parallel = 0;
                    index--;
                }
                else
                {
                    // 本MPCNS程序默认，目标块号为负的时候表示周期边界条件
                    temp_inp.is_parallel = 1;
                    index = -index;
                    index--;
                }

                for (int32_t i = 0; i < 3; i++)
                {
                    temp_inp.sub[i] = sub[i];
                    temp_inp.sup[i] = sup[i];
                    temp_inp.tar_sub[i] = tar_sub[i];
                    temp_inp.tar_sup[i] = tar_sup[i];
                }
                temp_inp.volume = (abs(abs(sub[0]) - abs(sup[0])) + 1) * (abs(abs(sub[1]) - abs(sup[1])) + 1) * (abs(abs(sub[2]) - abs(sup[2])) + 1);
                // 完成了一个coordinate_inp读入
                inp(iblock, index).inner.push_back(temp_inp);
            }
        }
    }
    file.close();

    std::cout << "\t-->Finish reading the inp file\n";
    //=============================================================================================
    filename = par->GetStr("ffilename");
    file.open(filename, std::ios_base::in);
    std::getline(file, line);                                   // FVBND 1 3
    file >> line;                                               //(0)No_boundary_Conditions,注意没有为1的类型
    phy_index.insert(std::pair<std::string, int32_t>(line, 0)); //(0)No_boundary_Conditions
    phy_name.push_back(line);                                   //(0)No_boundary_Conditions
    phy_name.push_back("");                                     //(1)""
    int bc_number = 1;
    while (true)
    {
        bc_number = bc_number + 1;
        file >> line; // 读入名称
        if (line == "BOUNDARIES")
        {
            break;
        }
        phy_index.insert(std::pair<std::string, int>(line, bc_number));
        phy_name.push_back(line);
    }
    file.close();

    std::cout << "\t-->Finish reading the fvbnd file\n";
    NormalizeInpPoleBoundaryKind(this);
    //=============================================================================================
    std::string str;
    for (int32_t i = 0; i < phy_name.size(); i++)
    {
        str = phy_name[i];
        if (str.size() < 7)
            continue;
        str = str.substr(0, 7);
        if (str == "PERIOD-" || str == "period_")
        {
            process_period_condition(par);
            break;
        }
    }
    std::cout << "\t-->Finish processing the period boundary condition\n";
    check_read_info(par);
}
//=================================================================================================

//=================================================================================================
//---------------------------------------处理周期边界条件---------------------------------------------

void Preprocess_Data_Structure::process_period_condition(Param *par)
{
    std::string str, str2, per_name;
    int32_t dimension = par->GetInt("dimension");
    int32_t index, index2, block_num, block_num2, block_num_index, block_num_index2;
    int32_t Transform[3];
    // index index2表示对应周期边界条件的编号
    // block_num block_num2表示周期边界条件所在的块的块号
    // block_num_index block_num_index2表示所在块的物理边界条件序号
    // Transform表示需要重构的对应关系即第一块j对应第二块-i 则Transform[1]=-1
    do
    {
        bool end_loop = true;
        // 首先找到PERIOD-开头的边界条件str，名称name部分为per_name，编号为block_num_index
        // 找到与周期边界条件对应的面所在位置
        block_num = -1;
        block_num2 = -1;
        block_num_index = -1;
        block_num_index2 = -1;
        index = -1;
        index2 = -1;
        for (int32_t iblock = 0; iblock < blk.Getsize1(); iblock++)
        {
            for (int32_t jface = 0; jface < blk(iblock).bound.size(); jface++)
            {
                int32_t i = blk(iblock).bound[jface].phy_kind;
                str = phy_name[i];
                if (str.size() < 7)
                    continue;
                str = str.substr(0, 7);
                if (str == "PERIOD-")
                {
                    // 获取该周期边界条件的名称，用于寻找对应的面
                    per_name = phy_name[i].substr(7, phy_name[i].size() - 1);
                    // 用str表示全称
                    str = phy_name[i];
                    // 记录周期边界条件的数字代号
                    index = i;
                    // 该周期边界条件在本块中的序号、块号
                    block_num = iblock;
                    block_num_index = jface;
                    // 阻止跳出循环，因为还有周期边界条件没有处理
                    end_loop = false;
                    break;
                }
            }
        }
        if (end_loop)
            break;
        // 然后找到与str对应的period边界条件
        for (int32_t iblock = 0; iblock < blk.Getsize1(); iblock++)
        {
            for (int32_t jface = 0; jface < blk(iblock).bound.size(); jface++)
            {
                int32_t i = blk(iblock).bound[jface].phy_kind;
                str2 = phy_name[i];
                if (str2.size() < 8)
                    continue;
                // 取出字符串period_后面的内容
                str2 = str2.substr(7, str2.size() - 1);
                // 若内容与该周期边界条件的名称per_name一致
                if (str2 == per_name && phy_name[i].substr(0, 7) == "period_")
                {
                    // 用str2表示全称
                    str2 = phy_name[i];
                    // 记录周期边界条件的数字代号
                    index2 = i;
                    // 该周期边界条件在本块中的序号、块号
                    block_num2 = iblock;
                    block_num_index2 = jface;
                    break;
                }
            }
        }
        if (block_num2 == -1)
        {
            std::cout << "Error occurs for the period boundary condition\t" << str << "\t" << str2 << std::endl;
            exit(-1);
        }

        std::cout << "\t\tNow processing the:\t" << str << "\t" << str2 << std::endl;
        //---------------------------------------------------------------------
        //---------------------------------------------------------------------
        // 在找到边界条件名称和编号后进行处理，主要是读入网格重建内部边界条件，删去物理边界条件
        // 首先建立temp_inp temp_inp_conj用于压入inner中
        coordinate_inp temp_inp, temp_inp_conj;
        coordinate_phy &bound = blk(block_num).bound[block_num_index];
        coordinate_phy &bound2 = blk(block_num2).bound[block_num_index2];
        int32_t sub[3], sup[3], tar_sup[3], tar_sub[3];
        for (int32_t i = 0; i < 3; i++)
        {
            sub[i] = bound.sub[i];
            sup[i] = bound.sup[i];
            tar_sub[i] = bound2.sub[i];
            tar_sup[i] = bound2.sup[i];
        }
        // 利用额外网格坐标点的xyz信息重建对应关系，存入Transform中
        find_corresponding_para(Transform, block_num, sub, sup, block_num2, tar_sub, tar_sup, par);
        //-------------------------------------------------------------------------------
        // 首先处理反向的边
        for (int32_t i = 0; i < 3; i++)
        {
            if (Transform[i] < 0)
            {
                int32_t temp;
                temp = tar_sub[abs(Transform[i]) - 1];
                tar_sub[abs(Transform[i]) - 1] = tar_sup[abs(Transform[i]) - 1];
                tar_sup[abs(Transform[i]) - 1] = temp;
            }
        }
        // 二维方向唯一确定（一个法向，另一个自然对应），而三维还需要一个负号来确定
        if (dimension == 3)
        {
            int32_t i, my_norm;
            // 找到法向
            for (int32_t ii = 0; ii < 3; ii++)
            {
                if (sub[ii] == sup[ii])
                {
                    my_norm = ii;
                    break;
                }
            }
            // 选择去掉法向外的随意一个方向，添加负号用于表征对应关系
            i = (my_norm == 0) ? 1 : 0;

            sub[i] = -sub[i];
            sup[i] = -sup[i];
            tar_sub[abs(Transform[i]) - 1] = -tar_sub[abs(Transform[i]) - 1];
            tar_sup[abs(Transform[i]) - 1] = -tar_sup[abs(Transform[i]) - 1];
        }
        // 至此sub sup tar_sub tar_sup已经完全调整好了，可以存入inp中
        //===============================================================================
        // 存入inp
        for (int32_t i = 0; i < 3; i++)
        {
            temp_inp.sub[i] = sub[i];
            temp_inp.sup[i] = sup[i];
            temp_inp.tar_sub[i] = tar_sub[i];
            temp_inp.tar_sup[i] = tar_sup[i];

            temp_inp_conj.sub[i] = tar_sub[i];
            temp_inp_conj.sup[i] = tar_sup[i];
            temp_inp_conj.tar_sub[i] = sub[i];
            temp_inp_conj.tar_sup[i] = sup[i];

            temp_inp.volume = (abs(sub[0] - sup[0]) + 1) * (abs(sub[1] - sup[1]) + 1) * (abs(sub[2] - sup[2]) + 1);
            temp_inp_conj.volume = temp_inp.volume;

            temp_inp.is_parallel = 1;      // 由周期边界条件产生的内部条件，虚网格不能传递，用1标记
            temp_inp_conj.is_parallel = 1; // 由周期边界条件产生的内部条件，虚网格不能传递，用1标记
        }
        inp(block_num, block_num2).inner.push_back(temp_inp);
        inp(block_num2, block_num).inner.push_back(temp_inp_conj);
        //================================================================
        // 删除原有的物理信息
        std::vector<coordinate_phy>::iterator iter = blk(block_num).bound.begin();
        iter = iter + block_num_index;
        // 移除第一个的周期边界条件
        blk(block_num).bound.erase(iter);
        // erase不会删除指针的对象或内存空间，如下释放内存
        // {
        //     std::vector<coordinate_phy> temp(blk(block_num).bound);
        //     temp.swap(blk(block_num).bound);
        // }
        // std::vector(blk(block_num).bound).swap(blk(block_num).bound);

        // 移除第二个的周期边界条件，erase不会删除指针的对象或内存空间，如下释放内存
        if (block_num2 == block_num && block_num_index < block_num_index2)
            block_num_index2--;
        std::vector<coordinate_phy>::iterator iter2 = blk(block_num2).bound.begin();
        iter2 = iter2 + block_num_index2;
        blk(block_num2).bound.erase(iter2);
        // {
        //     std::vector<coordinate_phy> temp(blk(block_num2).bound);
        //     temp.swap(blk(block_num2).bound);
        // }
        // std::vector(blk(block_num2).bound).swap(blk(block_num2).bound);

        // 移除物理边界名称到编号的map中元素
        //  Remark phy_kind其他不能变，因为其他序号表示其他的物理边界条件
        //       同理phy_name也不能变
        phy_index.erase(str);
        phy_index.erase(str2);

        // for (auto i = phy_index.begin(); i != phy_index.end(); i++)
        //     std::cout << i->first << "\t" << i->second << "\t\n";
        // std::cout << std::endl;

        std::cout << "\t\tWE process and delete the period boundary:\t" << str << "\t" << str2 << std::endl;

    } while (true);
    check_read_info(par);
}

void Preprocess_Data_Structure::find_corresponding_para(int32_t *Transform, int32_t block_num, int32_t *sub, int32_t *sup, int32_t block_num2, int32_t *tar_sub, int32_t *tar_sup, Param *par)
{
    int32_t dimension = par->GetInt("dimension");
    // 新建数组用于存储网格坐标信息
    double begin[3], end[3], mid[3], tar_begin[3], tar_end[3], tar_mid[3], tar_extra[3];
    // double xyz[8][3], xyz2[8][3];
    double **xyz, **xyz2;
    xyz = new double *[8];
    xyz2 = new double *[8];
    for (int i = 0; i < 8; i++)
    {
        xyz[i] = new double[3];
        xyz2[i] = new double[3];
    }
    // Remark 这里的begin表示block_num块的sub[0]sub[1]sub[2]的xyz坐标
    //               end表示block_num块的sup[0]sup[1]sup[2]的xyz坐标
    //               mid表示block_num块的sup[0]sub[1]sub[2]的xyz坐标(若0为法向则为sub[0]sup[1]sub[2])
    int32_t a_index[3], b_index[3];
    int32_t my_norm, tar_norm, temp = 0;
    read_grid_coordinate(par, block_num, sub, sup, xyz);
    read_grid_coordinate(par, block_num2, tar_sub, tar_sup, xyz2);

    //--------------------------------------------------------------------------------------------
    // 搜寻重建对应关系的思路如下
    // 存在法向，故而Transform只剩下两个自由度，分别可取-3 -2 -1 1 2 3，同时需要排除重复数字的情况，一共8种
    // 前面的begin end mid可以用数组a_index b_index表示。如a_index={-1,-1,-1}表示subsubsub即begin
    // a_index仅仅由正负号表示sub sup。
    // 在有Transform后由a_index到b_index的计算关系为b(abs(Transform[i])-1)=a(i)sgn(Transform[i])
    // 为了避免重复多次读取较大的网格文件，需要读取一次后存储，从而能够方便获取a_index到对应坐标的数组xyz &xyz2
    // 对应方法为a_index每一个数映射为0,1（-1为0,1为1），组成的二进制数作为第一个输入e.g.
    // a_index={-1,1,-1}==>{0 1 0}=2即xyz[2][0~2]=*(*(xyz+2)+0~2)
    //--------------------------------------------------------------------------------------------

    // Transform初始化为0
    for (int32_t i = 0; i < 3; i++)
        Transform[i] = 0;
    // 找到对应的法向
    for (int32_t i = 0; i < dimension; i++)
    {
        if (sub[i] == sup[i])
            my_norm = i;
        if (tar_sub[i] == tar_sup[i])
            tar_norm = i;
    }
    Transform[my_norm] = tar_norm + 1;
    for (int i = 0; i < 3; i++)
    {
        begin[i] = *(*(xyz + 0) + i);
        end[i] = *(*(xyz + 7) + i);
        // 对于mid仅仅取除去norm外的第一个方向为sup,其他为sub即1 0 0或0 1 0
        if (my_norm == 0)
        {
            mid[i] = *(*(xyz + 2) + i);
        }
        else
        {
            mid[i] = *(*(xyz + 4) + i);
        }
    }

    if (dimension == 3)
    {
        bool stop_loop = false;
        for (int32_t i = -3; i <= 3; i++)
        {
            for (int32_t j = -3; j <= 3; j++)
            {
                if (abs(i) - 1 == my_norm || abs(j) - 1 == my_norm || j == 0 || i == 0)
                    continue;
                if (abs(i) == abs(j))
                    continue;
                if (my_norm == 0)
                {
                    Transform[1] = i;
                    Transform[2] = j;
                }
                else
                {
                    Transform[0] = i;
                    Transform[3 - my_norm] = j;
                }
                // 至此实现了所有Transform的循环，下面检测即可
                std::cout << "\t\t  testing Transform " << Transform[0] << "\t" << Transform[1] << "\t" << Transform[2] << "\n";
                //-------------------------------------------------------------------------------
                // 验证向量是否平行，首先计算向量
                for (int i = 0; i < 3; i++)
                    a_index[i] = -1;
                for (int i = 0; i < 3; i++)
                    b_index[abs(Transform[i]) - 1] = (Transform[i] > 0) ? a_index[i] : -a_index[i];
                // 将b_index转化为所谓二进制，大小存入temp
                temp = (1 + b_index[0]) * 2 + (1 + b_index[1]) + (1 + b_index[2]) / 2;
                for (int i = 0; i < 3; i++)
                    tar_begin[i] = *(*(xyz2 + temp) + i);
                //-----------------------------------------------
                for (int i = 0; i < 3; i++)
                    a_index[i] = 1;
                for (int i = 0; i < 3; i++)
                    b_index[abs(Transform[i]) - 1] = (Transform[i] > 0) ? a_index[i] : -a_index[i];
                // 将b_index转化为所谓二进制，大小存入temp
                temp = (1 + b_index[0]) * 2 + (1 + b_index[1]) + (1 + b_index[2]) / 2;
                for (int i = 0; i < 3; i++)
                    tar_end[i] = *(*(xyz2 + temp) + i);
                //-----------------------------------------------
                for (int i = 0; i < 3; i++)
                    a_index[i] = -1;
                a_index[0] = 1;
                if (my_norm == 0)
                    a_index[1] = 1;
                // 若法向为0方向，需要1方向变为sup
                for (int i = 0; i < 3; i++)
                    b_index[abs(Transform[i]) - 1] = (Transform[i] > 0) ? a_index[i] : -a_index[i];
                // 将b_index转化为所谓二进制，大小存入temp
                temp = (1 + b_index[0]) * 2 + (1 + b_index[1]) + (1 + b_index[2]) / 2;
                for (int i = 0; i < 3; i++)
                    tar_mid[i] = *(*(xyz2 + temp) + i);
                //-----------------------------------------------
                if (check_coordinate_parallel(begin, tar_begin, end, tar_end) && check_coordinate_parallel(mid, tar_mid, end, tar_end))
                {
                    stop_loop = true;
                    break;
                }
            }
            if (stop_loop)
                break;
        }
    }
    else
    {
        // 2D只有一个自由度
        Transform[2] = 3;

        Transform[my_norm] = tar_norm + 1;
        // 剩下最后一个方向,但是有可能是反向的，需要用坐标验证
        Transform[1 - my_norm] = 3 - Transform[my_norm];
        //-------------------------------------------------------------------------------
        // 假设现在是对的，验证向量是否平行
        for (int i = 0; i < 3; i++)
            a_index[i] = -1;
        for (int i = 0; i < 3; i++)
            b_index[abs(Transform[i]) - 1] = (Transform[i] > 0) ? a_index[i] : -a_index[i];
        // 将b_index转化为所谓二进制，大小存入temp
        temp = (1 + b_index[0]) * 2 + (1 + b_index[1]) + (1 + b_index[2]) / 2;
        for (int i = 0; i < 3; i++)
            tar_begin[i] = *(*(xyz2 + temp) + i);
        //-----------------------------------------------
        for (int i = 0; i < 3; i++)
            a_index[i] = 1;
        for (int i = 0; i < 3; i++)
            b_index[abs(Transform[i]) - 1] = (Transform[i] > 0) ? a_index[i] : -a_index[i];
        // 将b_index转化为所谓二进制，大小存入temp
        temp = (1 + b_index[0]) * 2 + (1 + b_index[1]) + (1 + b_index[2]) / 2;
        for (int i = 0; i < 3; i++)
            tar_end[i] = *(*(xyz2 + temp) + i);
        //-----------------------------------------------

        if (!check_coordinate_parallel(begin, tar_begin, end, tar_end))
        {
            // 若不平行，反向即可
            Transform[1 - my_norm] = Transform[my_norm] - 3;

            // 最后检验，防止输入网格以及周期边界条件本身就存在问题
            for (int i = 0; i < 3; i++)
                a_index[i] = -1;
            for (int i = 0; i < 3; i++)
                b_index[abs(Transform[i]) - 1] = (Transform[i] > 0) ? a_index[i] : -a_index[i];
            // 将b_index转化为所谓二进制，大小存入temp
            temp = (1 + b_index[0]) * 2 + (1 + b_index[1]) + (1 + b_index[2]) / 2;
            for (int i = 0; i < 3; i++)
                tar_begin[i] = *(*(xyz2 + temp) + i);
            //-----------------------------------------------
            for (int i = 0; i < 3; i++)
                a_index[i] = 1;
            for (int i = 0; i < 3; i++)
                b_index[abs(Transform[i]) - 1] = (Transform[i] > 0) ? a_index[i] : -a_index[i];
            // 将b_index转化为所谓二进制，大小存入temp
            temp = (1 + b_index[0]) * 2 + (1 + b_index[1]) + (1 + b_index[2]) / 2;
            for (int i = 0; i < 3; i++)
                tar_end[i] = *(*(xyz2 + temp) + i);

            if (!check_coordinate_parallel(begin, tar_begin, end, tar_end))
            {
                std::cout << "Error grids for parallel boundary_condition ! !" << std::endl;
                std::cout << "The Error blocks are\t" << block_num << "and\t" << block_num2 << std::endl;
                exit(-1);
            }
            return;
        }
    }
    for (int i = 0; i < 8; i++)
    {
        delete[] xyz[i];
        delete[] xyz2[i];
    }
    delete[] xyz;
    delete[] xyz2;
}

void Preprocess_Data_Structure::read_grid_coordinate(Param *par, int32_t blk_n, int32_t *sub, int32_t *sup, double **xyz)
{
    std::string filename, input;
    std::fstream file;
    filename = par->GetStr("gfilename");
    int32_t grd_readtype = par->GetInt("grd_readtype");

    double3D x, y, z;
    if (grd_readtype == 1 || grd_readtype == 2)
    {
        //==================================================
        // 读入原始网格文件,尝试读入gridgen的unformatted格式
        file.open(par->GetStr("gfilename"), std::ios::in | std::ios::binary);
        int32_t temp;
        if (grd_readtype == 1)
        {
            file.read((char *)&temp, sizeof(temp));
        }
        file.read((char *)&temp, sizeof(temp));
        int32_t numofblock = temp;
        if (grd_readtype == 1)
        {
            file.read((char *)&temp, sizeof(temp));
            file.read((char *)&temp, sizeof(temp));
        }

        int32_t *mx = new int32_t[numofblock];
        int32_t *my = new int32_t[numofblock];
        int32_t *mz = new int32_t[numofblock];
        for (int iiblock = 0; iiblock < numofblock; iiblock++)
        {
            file.read((char *)&temp, sizeof(temp));
            mx[iiblock] = temp;
            file.read((char *)&temp, sizeof(temp));
            my[iiblock] = temp;
            file.read((char *)&temp, sizeof(temp));
            mz[iiblock] = temp;
            // std::cout << " "<<mx[iiblock] << " " << my[iiblock] << " " << mz[iiblock]<<std::endl;
        }
        // 注意这里的mx my mz是临时的网格点数，且是从1开始的
        if (grd_readtype == 1)
        {
            file.read((char *)&temp, sizeof(temp));
        }

        x.SetSize(mx[blk_n], my[blk_n], mz[blk_n], 0);
        y.SetSize(mx[blk_n], my[blk_n], mz[blk_n], 0);
        z.SetSize(mx[blk_n], my[blk_n], mz[blk_n], 0);

        for (int iiblock = 0; iiblock < numofblock; iiblock++)
        {
            double temp_double;
            double temp_binary;
            if (grd_readtype == 1)
            {
                file.read((char *)&temp, sizeof(temp));
            }
            if (iiblock == blk_n)
            {
                for (int k = 0; k < mz[iiblock]; ++k)
                    for (int j = 0; j < my[iiblock]; ++j)
                        for (int i = 0; i < mx[iiblock]; ++i)
                        {
                            file.read((char *)&temp_binary, sizeof(temp_binary));
                            x(i, j, k) = temp_binary;
                            // std::cout << " " << temp_binary << std::endl;
                        }
                for (int k = 0; k < mz[iiblock]; ++k)
                    for (int j = 0; j < my[iiblock]; ++j)
                        for (int i = 0; i < mx[iiblock]; ++i)
                        {
                            file.read((char *)&temp_binary, sizeof(temp_binary));
                            y(i, j, k) = temp_binary;
                            // std::cout << " " << temp_binary << std::endl;
                        }
                for (int k = 0; k < mz[iiblock]; ++k)
                    for (int j = 0; j < my[iiblock]; ++j)
                        for (int i = 0; i < mx[iiblock]; ++i)
                        {
                            file.read((char *)&temp_binary, sizeof(temp_binary));
                            z(i, j, k) = temp_binary;
                            // std::cout << " " << temp_binary << std::endl;
                        }
                break;
            }
            else
            {
                for (int k = 0; k < mz[iiblock]; ++k)
                    for (int j = 0; j < my[iiblock]; ++j)
                        for (int i = 0; i < mx[iiblock]; ++i)
                        {
                            file.read((char *)&temp_binary, sizeof(temp_binary));
                            temp_double = temp_binary;
                            // std::cout << " " << temp_binary << std::endl;
                        }
                for (int k = 0; k < mz[iiblock]; ++k)
                    for (int j = 0; j < my[iiblock]; ++j)
                        for (int i = 0; i < mx[iiblock]; ++i)
                        {
                            file.read((char *)&temp_binary, sizeof(temp_binary));
                            temp_double = temp_binary;
                            // std::cout << " " << temp_binary << std::endl;
                        }
                for (int k = 0; k < mz[iiblock]; ++k)
                    for (int j = 0; j < my[iiblock]; ++j)
                        for (int i = 0; i < mx[iiblock]; ++i)
                        {
                            file.read((char *)&temp_binary, sizeof(temp_binary));
                            temp_double = temp_binary;
                            // std::cout << " " << temp_binary << std::endl;
                        }
            }
            if (grd_readtype == 1)
            {
                file.read((char *)&temp, sizeof(temp));
            }
        }
        delete[] mx;
        delete[] my;
        delete[] mz;
    }
    else
    {
        //==================================================
        // 读入原始网格文件,尝试读入gridgen的文本格式
        file.open(filename, std::ios_base::in);

        /****************读取网格块数*****************/
        int32_t numofblock = 0;
        file >> numofblock;
        /***************读取网格点数****************/
        int32_t *mx = new int32_t[numofblock];
        int32_t *my = new int32_t[numofblock];
        int32_t *mz = new int32_t[numofblock];
        int32_t dimension = 3;
        for (int iblock = 0; iblock < numofblock; iblock++)
            file >> mx[iblock] >> my[iblock] >> mz[iblock];
        // 注意这里的mx my mz是临时的网格点数，且是从1开始的

        x.SetSize(mx[blk_n], my[blk_n], mz[blk_n], 0);
        y.SetSize(mx[blk_n], my[blk_n], mz[blk_n], 0);
        z.SetSize(mx[blk_n], my[blk_n], mz[blk_n], 0);
        double temp = 0.0;

        for (int izone = 0; izone < numofblock; izone++)
        {
            if (izone != blk_n)
            {
                /*读取grd坐标*/
                for (int k = 0; k < mz[izone]; ++k)
                    for (int j = 0; j < my[izone]; ++j)
                        for (int i = 0; i < mx[izone]; ++i)
                            file >> temp;
                for (int k = 0; k < mz[izone]; ++k)
                    for (int j = 0; j < my[izone]; ++j)
                        for (int i = 0; i < mx[izone]; ++i)
                            file >> temp;
                for (int k = 0; k < mz[izone]; ++k)
                    for (int j = 0; j < my[izone]; ++j)
                        for (int i = 0; i < mx[izone]; ++i)
                            file >> temp;
            }
            else
            {
                /*读取grd坐标*/
                for (int k = 0; k < mz[izone]; ++k)
                    for (int j = 0; j < my[izone]; ++j)
                        for (int i = 0; i < mx[izone]; ++i)
                            file >> x(i, j, k);

                for (int k = 0; k < mz[izone]; ++k)
                    for (int j = 0; j < my[izone]; ++j)
                        for (int i = 0; i < mx[izone]; ++i)
                            file >> y(i, j, k);

                for (int k = 0; k < mz[izone]; ++k)
                    for (int j = 0; j < my[izone]; ++j)
                        for (int i = 0; i < mx[izone]; ++i)
                            file >> z(i, j, k);
                break;
            }
        }
        delete[] mx;
        delete[] my;
        delete[] mz;
    }
    int32_t ii[2], jj[2], kk[2];
    ii[0] = abs(sub[0]) - 1;
    ii[1] = abs(sup[0]) - 1;

    jj[0] = abs(sub[1]) - 1;
    jj[1] = abs(sup[1]) - 1;

    kk[0] = abs(sub[2]) - 1;
    kk[1] = abs(sup[2]) - 1;
    // 二维数组指针取值xyz[a][b]=*(*(xyz+a)+b)
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
                *(*(xyz + i * 4 + j * 2 + k) + 0) = x(ii[i], jj[j], kk[k]);
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
                *(*(xyz + i * 4 + j * 2 + k) + 1) = y(ii[i], jj[j], kk[k]);
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
                *(*(xyz + i * 4 + j * 2 + k) + 2) = z(ii[i], jj[j], kk[k]);

    file.close();
}

bool Preprocess_Data_Structure::check_coordinate_parallel(double *begin, double *tar_begin, double *end, double *tar_end)
{
    double r1[3], r2[3];
    for (int32_t i = 0; i < 3; i++)
    {
        r1[i] = tar_begin[i] - begin[i];
        r2[i] = tar_end[i] - end[i];
    }
    double temp = 0.0, rr1 = 0.0, rr2 = 0.0;
    for (int32_t i = 0; i < 3; i++)
    {
        temp += r1[i] * r2[i];
        rr1 += r1[i] * r1[i];
        rr2 += r2[i] * r2[i];
    }
    temp = temp / sqrt(rr1 * rr2);
    return (abs(1 - abs(temp)) < SMALL);
}

//---------------------------------------处理周期边界条件---------------------------------------------
//=================================================================================================

void Preprocess_Data_Structure::check_read_info(Param *par)
{
    int32_t blk_num = blk.Getsize1();
    int32_t dimension = par->GetInt("dimension");
    for (int32_t iblock = 0; iblock < blk_num; iblock++)
    {

        int32_t boundary_num = blk(iblock).bound.size();
        for (int32_t ibound = 0; ibound < boundary_num; ibound++)
        {
            coordinate_phy &boundary = blk(iblock).bound[ibound];
            if (phy_name[boundary.phy_kind].substr(0, 6) == "PERIOD-" || phy_name[boundary.phy_kind].substr(0, 6) == "period_")
            {
                std::cout << "Error appears for parallel boundary_condition processing! !" << std::endl;
                exit(-1);
            }
        }
    }
    //===================================================================================

    for (int32_t iblk = 0; iblk < blk_num; iblk++)
    {
        for (int32_t jblk = 0; jblk < blk_num; jblk++)
        {
            int32_t face_num = inp(iblk, jblk).inner.size();
            for (int32_t kface = 0; kface < face_num; kface++)
            {
                bool find_it = false;
                coordinate_inp &AA = inp(iblk, jblk).inner[kface];
                for (int32_t kkface = 0; kkface < inp(jblk, iblk).inner.size(); kkface++)
                {
                    coordinate_inp &BB = inp(jblk, iblk).inner[kkface];
                    if (AA.tar_sub[0] == BB.sub[0] && AA.tar_sub[2] == BB.sub[2] && AA.tar_sub[1] == BB.sub[1] && AA.tar_sup[0] == BB.sup[0] && AA.tar_sup[2] == BB.sup[2] && AA.tar_sup[1] == BB.sup[1])
                    {
                        find_it = true;
                        break;
                    }
                }
                if (!find_it)
                {
                    for (int32_t kkface = 0; kkface < inp(jblk, iblk).inner.size(); kkface++)
                    {
                        coordinate_inp &BB = inp(jblk, iblk).inner[kkface];
                        if (abs(abs(AA.tar_sub[0]) - abs(AA.tar_sup[0])) == abs(abs(BB.sub[0]) - abs(BB.sup[0])) && abs(abs(AA.tar_sub[1]) - abs(AA.tar_sup[1])) == abs(abs(BB.sub[1]) - abs(BB.sup[1])) && abs(abs(AA.tar_sub[2]) - abs(AA.tar_sup[2])) == abs(abs(BB.sub[2]) - abs(BB.sup[2])) && fmin(abs(AA.tar_sub[0]), abs(AA.tar_sup[0])) == fmin(abs(BB.sub[0]), abs(BB.sup[0])) && fmin(abs(AA.tar_sub[1]), abs(AA.tar_sup[1])) == fmin(abs(BB.sub[1]), abs(BB.sup[1])) && fmin(abs(AA.tar_sub[2]), abs(AA.tar_sup[2])) == fmin(abs(BB.sub[2]), abs(BB.sup[2])))
                        {
                            // 存在对应的块，应进一步检查AA.sub sup与BB.tar_sub tar_sup对应关系是否一致
                            int32_t TransA[3], TransB[3];
                            if (dimension == 2)
                            {
                                // 对于2维的情况
                                TransA[2] = 3;
                                TransB[2] = 3;
                                // 获得AA的转化关系
                                for (int32_t iii = 0; iii < 2; iii++)
                                {
                                    if (AA.tar_sub[iii] == AA.tar_sup[iii])
                                    {
                                        for (int32_t jjj = 0; jjj < 2; jjj++)
                                        {
                                            if (AA.sub[jjj] == AA.sup[jjj])
                                            {
                                                TransA[iii] = jjj + 1;
                                                if ((abs(AA.tar_sub[1 - iii]) - abs(AA.tar_sup[1 - iii])) * (abs(AA.sub[1 - jjj]) - abs(AA.sup[1 - jjj])) > 0)
                                                    TransA[1 - iii] = 3 - TransA[iii];
                                                else
                                                    TransA[1 - iii] = TransA[iii] - 3;

                                                break;
                                            }
                                        }
                                        break;
                                    }
                                }
                                // 获得BB中与AA对应的转换关系
                                for (int32_t iii = 0; iii < 2; iii++)
                                {
                                    if (BB.sub[iii] == BB.sup[iii])
                                    {
                                        for (int32_t jjj = 0; jjj < 2; jjj++)
                                        {
                                            if (BB.tar_sub[jjj] == BB.tar_sup[jjj])
                                            {
                                                TransB[iii] = jjj + 1;
                                                if ((abs(BB.sub[1 - iii]) - abs(BB.sup[1 - iii])) * (abs(BB.tar_sub[1 - jjj]) - abs(BB.tar_sup[1 - jjj])) > 0)
                                                    TransB[1 - iii] = 3 - TransB[iii];
                                                else
                                                    TransB[1 - iii] = TransB[iii] - 3;

                                                break;
                                            }
                                        }
                                        break;
                                    }
                                }
                                // 检测是否匹配
                                if (TransA[0] == TransB[0] && TransA[1] == TransB[1] && TransA[2] == TransB[2])
                                {
                                    // 调整BB形成对称
                                    for (int32_t iii = 0; iii < 2; iii++)
                                    {
                                        BB.sub[iii] = AA.tar_sub[iii];
                                        BB.sup[iii] = AA.tar_sup[iii];
                                        BB.tar_sub[iii] = AA.sub[iii];
                                        BB.tar_sup[iii] = AA.sup[iii];
                                    }
                                    find_it = true;
                                }
                            }
                            else
                            {
                                TransA[0] = 0;
                                TransA[1] = 0;
                                TransA[2] = 0;
                                TransB[0] = 0;
                                TransB[1] = 0;
                                TransB[2] = 0;

                                // 找到AA的法向
                                for (int32_t iii = 0; iii < 3; iii++)
                                {
                                    if (AA.tar_sub[iii] == AA.tar_sup[iii])
                                    {
                                        for (int32_t jjj = 0; jjj < 3; jjj++)
                                        {
                                            if (AA.sub[jjj] == AA.sup[jjj])
                                            {
                                                TransA[iii] = jjj + 1;
                                                break;
                                            }
                                        }
                                        break;
                                    }
                                }
                                // 找到AA用负号标识的方向
                                for (int32_t iii = 0; iii < 3; iii++)
                                {
                                    if (AA.tar_sub[iii] < 0)
                                    {
                                        for (int32_t jjj = 0; jjj < 3; jjj++)
                                        {
                                            if (AA.sub[jjj] < 0)
                                            {
                                                if ((abs(AA.tar_sub[iii]) - abs(AA.tar_sup[iii])) * (abs(AA.sub[jjj]) - abs(AA.sup[jjj])) > 0)
                                                    TransA[iii] = jjj + 1;
                                                else
                                                    TransA[iii] = -jjj - 1;
                                                break;
                                            }
                                        }
                                        break;
                                    }
                                }
                                // 获得AA剩下的方向
                                for (int32_t iii = 0; iii < 3; iii++)
                                {
                                    if (TransA[iii] == 0)
                                    {
                                        int32_t jjj = 5 - abs(TransA[0]) - abs(TransA[1]) - abs(TransA[2]);
                                        if ((abs(AA.tar_sub[iii]) - abs(AA.tar_sup[iii])) * (abs(AA.sub[jjj]) - abs(AA.sup[jjj])) > 0)
                                            TransA[iii] = jjj + 1;
                                        else
                                            TransA[iii] = -jjj - 1;
                                        break;
                                    }
                                }
                                // 找到BB的法向
                                for (int32_t iii = 0; iii < 3; iii++)
                                {
                                    if (BB.sub[iii] == BB.sup[iii])
                                    {
                                        for (int32_t jjj = 0; jjj < 3; jjj++)
                                        {
                                            if (BB.tar_sub[jjj] == BB.tar_sup[jjj])
                                            {
                                                TransB[iii] = jjj + 1;
                                                break;
                                            }
                                        }
                                        break;
                                    }
                                }
                                // 找到BB用负号标识的方向
                                for (int32_t iii = 0; iii < 3; iii++)
                                {
                                    if (BB.sub[iii] < 0)
                                    {
                                        for (int32_t jjj = 0; jjj < 3; jjj++)
                                        {
                                            if (BB.tar_sub[jjj] < 0)
                                            {
                                                if ((abs(BB.sub[iii]) - abs(BB.sup[iii])) * (abs(BB.tar_sub[jjj]) - abs(BB.tar_sup[jjj])) > 0)
                                                    TransB[iii] = jjj + 1;
                                                else
                                                    TransB[iii] = -jjj - 1;
                                                break;
                                            }
                                        }
                                        break;
                                    }
                                }
                                // 获得BB剩下的方向
                                for (int32_t iii = 0; iii < 3; iii++)
                                {
                                    if (TransB[iii] == 0)
                                    {
                                        int32_t jjj = 5 - abs(TransB[0]) - abs(TransB[1]) - abs(TransB[2]);
                                        if ((abs(BB.sub[iii]) - abs(BB.sup[iii])) * (abs(BB.tar_sub[jjj]) - abs(BB.tar_sup[jjj])) > 0)
                                            TransB[iii] = jjj + 1;
                                        else
                                            TransB[iii] = -jjj - 1;
                                        break;
                                    }
                                }

                                // 检测是否匹配
                                if (TransA[0] == TransB[0] && TransA[1] == TransB[1] && TransA[2] == TransB[2])
                                {
                                    // 调整BB形成对称
                                    std::cout << "\t\t#We Adjust the inp information:\t->Block\t\t" << jblk << "\tto\t" << iblk << "\tFace num is\t" << kkface << std::endl;
                                    std::cout << "\t\t" << BB.sub[0] << "\t" << BB.sup[0] << "\t" << BB.sub[1] << "\t" << BB.sup[1] << "\t" << BB.sub[2] << "\t" << BB.sup[2] << std::endl;
                                    std::cout << "\t\t" << BB.tar_sub[0] << "\t" << BB.tar_sup[0] << "\t" << BB.tar_sub[1] << "\t" << BB.tar_sup[1] << "\t" << BB.tar_sub[2] << "\t" << BB.tar_sup[2] << std::endl;
                                    std::cout << "\t\tinto the inp information:\n";
                                    for (int32_t iii = 0; iii < 3; iii++)
                                    {
                                        BB.sub[iii] = AA.tar_sub[iii];
                                        BB.sup[iii] = AA.tar_sup[iii];
                                        BB.tar_sub[iii] = AA.sub[iii];
                                        BB.tar_sup[iii] = AA.sup[iii];
                                    }
                                    std::cout << "\t\t" << BB.sub[0] << "\t" << BB.sup[0] << "\t" << BB.sub[1] << "\t" << BB.sup[1] << "\t" << BB.sub[2] << "\t" << BB.sup[2] << std::endl;
                                    std::cout << "\t\t" << BB.tar_sub[0] << "\t" << BB.tar_sup[0] << "\t" << BB.tar_sub[1] << "\t" << BB.tar_sup[1] << "\t" << BB.tar_sub[2] << "\t" << BB.tar_sup[2] << std::endl
                                              << std::endl;
                                    find_it = true;
                                }
                            }
                            break;
                        }
                    }
                    if (find_it)
                        continue;
                    std::cout << "Error occurs in INP data structure, it is not symmetry!!\n";
                    std::cout << "Block " << iblk << "to " << jblk << "Face num is " << kface << std::endl;
                    std::cout << AA.sub[0] << "\t" << AA.sup[0] << "\t" << AA.sub[1] << "\t" << AA.sup[1] << "\t" << AA.sub[2] << "\t" << AA.sup[2] << std::endl;
                    std::cout << AA.tar_sub[0] << "\t" << AA.tar_sup[0] << "\t" << AA.tar_sub[1] << "\t" << AA.tar_sup[1] << "\t" << AA.tar_sub[2] << "\t" << AA.tar_sup[2] << std::endl;

                    exit(-1);
                }
            }
        }
    }
    //===================================================================================
    // #################NOT FINISHED#####################
    // 检查是否所有的边界都有边界条件
    // #################NOT FINISHED#####################
}

/*
 * @brief 对本数据结构的blk_number块按照direction方向的location位置进行剖分自动剖分，更新连接关系
 */
void Preprocess_Data_Structure::split_block(Param *par, int32_t blk_number, int32_t direction, int32_t location)
{
    // Remark
    //       剖分程序按照如下的思路进行
    //       剖分一个网格块将产生一个新的网格块，我们默认direction方向小于location的仍然用blk_number标号记录，而
    // direction方向大于location的放在blk的最后位置。这样能够尽可能减少修正连接关系。
    //       针对父子关系即ptr->blk[blk_number]->hiera中所有信息需要重构，其中blk_name，oring_num直接继承，剖分过程
    // 不会改变网格的方向和正负，故而Transform也可以直接继承，需要根据切割的位置和方向对坐标的范围对sub，sup, oring_sub
    // 以及oring_sup进行处理，主要是sub【direction】～sup【direction】范围中出现了location位置，就需要进行处理。例如，
    // 某一个块有多个hiera来源，某一个sub[direction]<location<sup[direction]则sub[direction]<location放在本块，而
    // location<sup[direction]放在新加块中，其他未被分割的则只需要判断direction方向的坐标范围是在location之上还是之下来
    // 直接修正对应的sub sup坐标即可。
    //       针对物理边界条件ptr->blk[blk_number]->bound也是只需要改变坐标范围，同上；phy_kind直接继承即可
    //       针对最为复杂的内部边界条件需要慎重处理。ptr->inp(blk_number,:),还要新建一个扩充个数 tar_num my_num容易继承
    // 而剩下的inner需要仔细考量。主要有sub sup tar_sub tar_sup。首先根据location判断各个面放在哪一个block，一旦决定后进
    // 行处理即可。对于放在本块的不需要作任何处理，放在最后新加的一个网格块的需要修改direction方向的坐标上下限，同时修改目标块的
    // 对应信息。
    //=============================================================================================
    //   剖分开始
    //=============================================================================================
    Preprocess_Data_Structure *ptr = this;
    int32_t dimension = par->GetInt("dimension");
    // （0）扩大空间
    Array2D<Inner> _inp_temp = ptr->inp;
    Array1D<block_info> _blk_temp = ptr->blk;
    int32_t number_of_block = _blk_temp.Getsize1();
    ptr->inp.SetSize(number_of_block + 1, number_of_block + 1);
    ptr->blk.SetSize(number_of_block + 1);

    // for (int32_t i = 0; i <= number_of_block; i++)
    // {
    //     ptr->blk(i).bound.clear();
    //     ptr->blk(i).hiera.clear();
    //     for (int32_t j = 0; j <= number_of_block; j++)
    //         ptr->inp(i, j).inner.clear();
    // }

    for (int32_t i = 0; i < number_of_block; i++)
    {
        ptr->blk(i) = _blk_temp(i);
        for (int32_t j = 0; j < number_of_block; j++)
        {
            ptr->inp(i, j) = _inp_temp(i, j);
        }
    }
    for (int32_t i = 0; i <= number_of_block; i++)
    {
        ptr->inp(number_of_block, i).my_num = number_of_block;
        ptr->inp(number_of_block, i).tar_num = i;
        ptr->inp(i, number_of_block).my_num = i;
        ptr->inp(i, number_of_block).tar_num = number_of_block;
    }
    //=============================================================================================
    // （1）针对父子关系进行处理
    std::vector<coordinate_blk> &hier = ptr->blk(blk_number).hiera;
    std::vector<coordinate_blk> &hier_new = ptr->blk(number_of_block).hiera;
    int32_t sub_block_number = hier.size();
    // 第一次循环只会把大于location的加入hier_new中，把要分割的加入hier_new中，hier需要删除的用remove_index记录
    std::vector<int32_t> remove_index;
    for (int32_t i_subblock = 0; i_subblock < sub_block_number; i_subblock++)
    {
        if (hier[i_subblock].sub[direction - 1] <= location && hier[i_subblock].sup[direction - 1] <= location)
        {
            // 该块的对应方向均小于切割位置location，仍然存在本块中，故而不变
            continue;
        }
        else if (hier[i_subblock].sub[direction - 1] >= location && hier[i_subblock].sup[direction - 1] >= location)
        {
            // 该块的对应方向均大于切割位置location，存在新的网格块中
            //-------------------------------------------------------------------------------------
            //  hier应当删除，计入remove_index
            int32_t index = i_subblock;
            remove_index.push_back(index);
            //-------------------------------------------------------------------------------------
            // 处理hier_new
            coordinate_blk new_subblock(hier[i_subblock]);
            // 由于放入新的网格块中，sub sup需要变化
            // 更新sub sup
            new_subblock.sup[direction - 1] = new_subblock.sup[direction - 1] - location + 1;
            new_subblock.sub[direction - 1] = new_subblock.sub[direction - 1] - location + 1;
            // 推入hier_new中
            hier_new.push_back(new_subblock);
        }
        else if (hier[i_subblock].sub[direction - 1] < location && hier[i_subblock].sup[direction - 1] > location)
        {
            // 该块的对应方向大于切割位置location的存在新的网格块中
            coordinate_blk new_subblock(hier[i_subblock]);
            // 由于放入新的网格块中，sub sup oring_sub oring_sup需要变化
            int32_t origin_dir = new_subblock.Transform[direction - 1];
            if (origin_dir > 0)
            {
                // 本块与原始块同向，因此sub增大即可
                new_subblock.oring_sub[abs(origin_dir) - 1] += (location - new_subblock.sub[direction - 1]);
            }
            else
            {
                // 本块与原始块异向，因此sup增大即可
                // 注意本块ijk认为是正向的，即sub<sup严格成立，若Transform为负表示目标sup<sub,在纯分割过程中不会出现
                new_subblock.oring_sup[abs(origin_dir) - 1] += (location - new_subblock.sub[direction - 1]);
            }
            // 更新sub sup
            new_subblock.sup[direction - 1] = new_subblock.sup[direction - 1] - location + 1;
            new_subblock.sub[direction - 1] = 1;
            // 推入hier_new中
            hier_new.push_back(new_subblock);
            //-------------------------------------------------------------------------------------
            // hier也要修正
            if (origin_dir > 0)
            {
                // 本块与原始块同向，因此sup减小即可
                hier[i_subblock].oring_sup[abs(origin_dir) - 1] -= (hier[i_subblock].sup[direction - 1] - location);
            }
            else
            {
                // 本块与原始块异向，因此sub减小即可
                // 注意本块ijk认为是正向的，即sub<sup严格成立，若Transform为负表示目标sup<sub,在纯分割过程中不会出现
                hier[i_subblock].oring_sub[abs(origin_dir) - 1] -= (hier[i_subblock].sup[direction - 1] - location);
            }
            // 更新sub sup
            hier[i_subblock].sup[direction - 1] = location;
        }
    }
    // 第二次循环会把全部大于location的hier删除,标号记录在remove_index中
    for (int32_t i_subblock = 0; i_subblock < remove_index.size(); i_subblock++)
    {
        // 该块的对应方向均大于切割位置location，存在新的网格块中，本块直接删除
        //  erase不会删除指针的对象或内存空间，如下释放内存
        std::vector<coordinate_blk>::iterator iter2 = hier.begin();
        iter2 = iter2 + remove_index[i_subblock];
        hier.erase(iter2);
        ////如果存在带指针的vector可以考虑如下形式
        // {
        //     std::vector<coordinate_blk> temp(hier);
        //     temp.swap(hier);
        // }
        // 由于删除了一个元素，remove_index中在i_subblock之后的标号需要减去1
        for (int32_t ii = i_subblock + 1; ii < remove_index.size(); ii++)
            remove_index[ii]--;
    }
    // 修改blk的尺寸信息

    for (int32_t i = 0; i < 3; i++)
        ptr->blk(number_of_block).ijkmax[i] = ptr->blk(blk_number).ijkmax[i];
    ptr->blk(blk_number).ijkmax[direction - 1] = location;
    ptr->blk(number_of_block).ijkmax[direction - 1] -= (location - 1);
    //=============================================================================================
    // （2）针对物理边界条件进行处理
    std::vector<coordinate_phy> &phy = ptr->blk(blk_number).bound;
    std::vector<coordinate_phy> &phy2 = ptr->blk(number_of_block).bound;
    for (int32_t iface = 0; iface < phy.size(); iface++)
    {
        if (phy[iface].sub[direction - 1] <= location && phy[iface].sup[direction - 1] <= location)
        {
            // 该块的对应方向均小于切割位置location，仍然存在本块中，故而不变
            continue;
        }
        else if (phy[iface].sub[direction - 1] >= location && phy[iface].sup[direction - 1] >= location)
        {
            // 该块的对应方向均大于切割位置location，存在新的网格块中
            //-------------------------------------------------------------------------------------
            //  处理phy2
            coordinate_phy phy_new(phy[iface]);
            phy_new.sub[direction - 1] = phy_new.sub[direction - 1] - location + 1;
            phy_new.sup[direction - 1] = phy_new.sup[direction - 1] - location + 1;
            phy2.push_back(phy_new);
            //-------------------------------------------------------------------------------------
            // phy删除本个元素
            std::vector<coordinate_phy>::iterator iter2 = phy.begin();
            iter2 = iter2 + iface;
            phy.erase(iter2);
            // 由于删去了一个变量，需要将标号前移一个，继续循环
            iface--;
        }
        else if ((phy[iface].sub[direction - 1] - location) * (phy[iface].sup[direction - 1] - location) < 0)
        {
            // 该块的对应方向跨越了切割位置location
            //-------------------------------------------------------------------------------------
            //  处理phy2
            coordinate_phy phy_new(phy[iface]);
            phy_new.sup[direction - 1] = (phy_new.sup[direction - 1] > phy_new.sub[direction - 1]) ? phy_new.sup[direction - 1] - location + 1 : phy_new.sub[direction - 1] - location + 1;
            phy_new.sub[direction - 1] = 1;
            phy2.push_back(phy_new);
            //-------------------------------------------------------------------------------------
            // 处理phy
            if (phy[iface].sup[direction - 1] > phy[iface].sub[direction - 1])
            {
                phy[iface].sup[direction - 1] = location;
            }
            else
            {
                phy[iface].sub[direction - 1] = location;
            }
        }
    }
    //=============================================================================================
    // （3）针对内部条件进行处理
    // Remark:在增加新块前，一共有number_of_block个块，0～number_of_block-1，增加后为0～number_of_block
    // 这里循环到number_of_block
    // 每次处理两件事情，第一本块的信息，第二目标块的信息；最后增加额外的切面
    for (int32_t iblock = 0; iblock <= number_of_block; iblock++)
    {
        std::vector<coordinate_inp> &inp = ptr->inp(blk_number, iblock).inner;
        std::vector<coordinate_inp> &inp2 = ptr->inp(number_of_block, iblock).inner;
        for (int32_t iface = 0; iface < inp.size(); iface++)
        {
            if (abs(inp[iface].sub[direction - 1]) <= location && abs(inp[iface].sup[direction - 1]) <= location)
            {
                // 该块的对应方向均小于切割位置location，仍然存在本块中，故而不变
                continue;
            }
            else if (abs(inp[iface].sub[direction - 1]) >= location && abs(inp[iface].sup[direction - 1]) >= location)
            {
                // 该块的对应方向均大于切割位置location，存在新块中，本块删除
                //---------------------------------------------------------------------------------
                // 首先处理inp2,以及目标块与本块的关系，故而一共有四个需要处理
                coordinate_inp inp_temp(inp[iface]);
                coordinate_inp inp_temp_conj(inp[iface]);

                inp_temp.sub[direction - 1] = (inp_temp.sub[direction - 1] > 0) ? abs(inp_temp.sub[direction - 1]) - location + 1 : -abs(inp_temp.sub[direction - 1]) + location - 1;
                inp_temp.sup[direction - 1] = (inp_temp.sup[direction - 1] > 0) ? abs(inp_temp.sup[direction - 1]) - location + 1 : -abs(inp_temp.sup[direction - 1]) + location - 1;
                for (int32_t i = 0; i < 3; i++)
                {
                    inp_temp_conj.sub[i] = inp_temp.tar_sub[i];
                    inp_temp_conj.sup[i] = inp_temp.tar_sup[i];
                    inp_temp_conj.tar_sub[i] = inp_temp.sub[i];
                    inp_temp_conj.tar_sup[i] = inp_temp.sup[i];
                }
                inp2.push_back(inp_temp);
                ptr->inp(iblock, number_of_block).inner.push_back(inp_temp_conj);
                //---------------------------------------------------------------------------------
                // 1、删除本块此边界条件的信息
                // 删除前记录tar_信息
                int32_t record_tar_sup[3], record_tar_sub[3];
                for (int iii = 0; iii < 3; iii++)
                {
                    record_tar_sub[iii] = inp[iface].tar_sub[iii];
                    record_tar_sup[iii] = inp[iface].tar_sup[iii];
                }
                // 开始删除本块此边界条件的信息
                std::vector<coordinate_inp>::iterator iter2 = ptr->inp(blk_number, iblock).inner.begin();
                iter2 = iter2 + iface;
                ptr->inp(blk_number, iblock).inner.erase(iter2);
                // 由于删去了此面的一个变量，需要将标号前移一个，继续循环
                iface--;

                // 2、删除对称位置的信息
                iter2 = ptr->inp(iblock, blk_number).inner.begin();
                // 找到对应块中的界面序号
                int32_t face_index = find_index_inner_face(&ptr->inp(iblock, blk_number), record_tar_sub, record_tar_sup);
                // 删除对应块中的信息
                iter2 = iter2 + face_index;
                ptr->inp(iblock, blk_number).inner.erase(iter2);
            }
            else if ((abs(inp[iface].sub[direction - 1]) - location) * (abs(inp[iface].sup[direction - 1]) - location) < 0)
            {
                // 该块的对应方向跨越了切割位置location，存在新块中，本块修改
                //---------------------------------------------------------------------------------
                // 首先处理inp2,以及目标块与本块的关系，故而一共有四个需要处理
                coordinate_inp inp_temp(inp[iface]);
                coordinate_inp inp_temp_conj(inp[iface]);
                int32_t dir = 1; // 表示与direction对应块的方向以及正负
                // 找到对应关系
                if (dimension == 2)
                {
                    if (inp_temp.tar_sub[0] == inp_temp.tar_sup[0])
                    {
                        dir = 1;
                    }
                    else
                    {
                        dir = 0;
                    }
                    if ((abs(inp_temp.sub[direction - 1]) - abs(inp_temp.sup[direction - 1])) * (abs(inp_temp.tar_sub[dir]) - abs(inp_temp.tar_sup[dir])) < 0)
                    {
                        dir = -(dir + 1);
                    }
                    else
                    {
                        dir++;
                    }
                }
                else
                {
                    for (dir = 0; dir < dimension; dir++)
                    {
                        if (inp_temp.tar_sub[dir] == inp_temp.tar_sup[dir])
                            continue;
                        if (inp_temp.sub[direction - 1] * inp_temp.tar_sub[dir] > 0)
                        {
                            // 找到方向
                            if ((abs(inp_temp.sub[direction - 1]) - abs(inp_temp.sup[direction - 1])) * (abs(inp_temp.tar_sub[dir]) - abs(inp_temp.tar_sup[dir])) < 0)
                            {
                                dir = -(dir + 1);
                            }
                            else
                            {
                                dir++;
                            }
                            break;
                        }
                    }
                }

                int32_t _sub = 0, _sup = 0, ss_min, ss_max;
                ss_min = std::min(abs(inp_temp.sub[direction - 1]), abs(inp_temp.sup[direction - 1]));
                ss_max = std::max(abs(inp_temp.sub[direction - 1]), abs(inp_temp.sup[direction - 1]));
                _sub = (inp_temp.sub[direction - 1] > 0) ? 1 : -1;
                _sup = (inp_temp.sub[direction - 1] > 0) ? ss_max - location + 1 : -(ss_max - location + 1);
                inp_temp.sub[direction - 1] = (abs(inp_temp.sub[direction - 1]) == ss_min) ? _sub : _sup;
                inp_temp.sup[direction - 1] = (abs(inp_temp.sup[direction - 1]) == ss_min) ? _sub : _sup;

                int32_t tar_ss_min, tar_ss_max, tar_sup, tar_sub;
                tar_ss_min = std::min(abs(inp_temp.tar_sub[abs(dir) - 1]), abs(inp_temp.tar_sup[abs(dir) - 1]));
                tar_ss_max = std::max(abs(inp_temp.tar_sub[abs(dir) - 1]), abs(inp_temp.tar_sup[abs(dir) - 1]));
                if (dir > 0)
                {
                    tar_sub = (inp_temp.tar_sub[abs(dir) - 1] > 0) ? tar_ss_max - (ss_max - location) : -(tar_ss_max - (ss_max - location));
                    tar_sup = (inp_temp.tar_sub[abs(dir) - 1] > 0) ? tar_ss_max : -tar_ss_max;
                }
                else
                {
                    tar_sub = (inp_temp.tar_sub[abs(dir) - 1] > 0) ? tar_ss_min : -tar_ss_min;
                    tar_sup = (inp_temp.tar_sub[abs(dir) - 1] > 0) ? tar_ss_min + (ss_max - location) : -(tar_ss_min + (ss_max - location));
                }

                inp_temp.tar_sub[abs(dir) - 1] = (abs(inp_temp.tar_sub[abs(dir) - 1]) == tar_ss_min) ? tar_sub : tar_sup;
                inp_temp.tar_sup[abs(dir) - 1] = (abs(inp_temp.tar_sup[abs(dir) - 1]) == tar_ss_min) ? tar_sub : tar_sup;

                inp_temp.volume = (abs(abs(inp_temp.sub[0]) - abs(inp_temp.sup[0])) + 1) * (abs(abs(inp_temp.sub[1]) - abs(inp_temp.sup[1])) + 1) * (abs(abs(inp_temp.sub[2]) - abs(inp_temp.sup[2])) + 1);

                for (int32_t i = 0; i < 3; i++)
                {
                    inp_temp_conj.sub[i] = inp_temp.tar_sub[i];
                    inp_temp_conj.sup[i] = inp_temp.tar_sup[i];
                    inp_temp_conj.tar_sub[i] = inp_temp.sub[i];
                    inp_temp_conj.tar_sup[i] = inp_temp.sup[i];
                }
                inp_temp_conj.volume = inp_temp.volume;

                inp2.push_back(inp_temp);
                ptr->inp(iblock, number_of_block).inner.push_back(inp_temp_conj);
                //---------------------------------------------------------------------------------
                // 然后处理本块和本块的对称块
                _sub = (inp[iface].sub[direction - 1] > 0) ? ss_min : -ss_min;
                _sup = (inp[iface].sub[direction - 1] > 0) ? location : -location;
                if (dir > 0)
                {
                    tar_sub = (inp[iface].tar_sub[abs(dir) - 1] > 0) ? tar_ss_min : -tar_ss_min;
                    tar_sup = (inp[iface].tar_sub[abs(dir) - 1] > 0) ? tar_ss_min + location - ss_min : -(tar_ss_min + location - ss_min);
                }
                else
                {
                    tar_sub = (inp[iface].tar_sub[abs(dir) - 1] > 0) ? tar_ss_max - (location - ss_min) : -(tar_ss_max - (location - ss_min));
                    tar_sup = (inp[iface].tar_sub[abs(dir) - 1] > 0) ? tar_ss_max : -tar_ss_max;
                }
                int32_t tar_face = find_index_inner_face(&ptr->inp(iblock, blk_number), inp[iface].tar_sub, inp[iface].tar_sup);

                inp[iface].sub[direction - 1] = (abs(inp[iface].sub[direction - 1]) == ss_min) ? _sub : _sup;
                inp[iface].sup[direction - 1] = (abs(inp[iface].sup[direction - 1]) == ss_min) ? _sub : _sup;
                inp[iface].tar_sub[abs(dir) - 1] = (abs(inp[iface].tar_sub[abs(dir) - 1]) == tar_ss_min) ? tar_sub : tar_sup;
                inp[iface].tar_sup[abs(dir) - 1] = (abs(inp[iface].tar_sup[abs(dir) - 1]) == tar_ss_min) ? tar_sub : tar_sup;
                inp[iface].volume = (abs(abs(inp[iface].sub[0]) - abs(inp[iface].sup[0])) + 1) * (abs(abs(inp[iface].sub[1]) - abs(inp[iface].sup[1])) + 1) * (abs(abs(inp[iface].sub[2]) - abs(inp[iface].sup[2])) + 1);

                for (int32_t i = 0; i < 3; i++)
                {
                    ptr->inp(iblock, blk_number).inner[tar_face].sub[i] = inp[iface].tar_sub[i];
                    ptr->inp(iblock, blk_number).inner[tar_face].sup[i] = inp[iface].tar_sup[i];
                    ptr->inp(iblock, blk_number).inner[tar_face].tar_sub[i] = inp[iface].sub[i];
                    ptr->inp(iblock, blk_number).inner[tar_face].tar_sup[i] = inp[iface].sup[i];
                }
                ptr->inp(iblock, blk_number).inner[tar_face].volume = inp[iface].volume;
            }
        }
    }
    //=========================================================================================
    // 最后还需要添加切开的那一个面
    coordinate_inp inp_temp;
    coordinate_inp inp_temp_conj;
    inp_temp.is_parallel = 0;
    inp_temp_conj.is_parallel = 0;
    inp_temp.sub[direction - 1] = location;
    inp_temp.sup[direction - 1] = location;
    inp_temp.tar_sub[direction - 1] = 1;
    inp_temp.tar_sup[direction - 1] = 1;
    for (int32_t i = 0; i < 3; i++)
    {
        if (i == direction - 1)
            continue;
        inp_temp.sub[i] = 1;
        inp_temp.sup[i] = ptr->blk(blk_number).ijkmax[i];
        inp_temp.tar_sub[i] = 1;
        inp_temp.tar_sup[i] = ptr->blk(blk_number).ijkmax[i];
    }
    if (dimension == 3)
    {
        int temp_int = 0;
        if (direction - 1 == 0)
            temp_int = 1;
        inp_temp.sub[temp_int] = -1;
        inp_temp.sup[temp_int] = -ptr->blk(blk_number).ijkmax[temp_int];
        inp_temp.tar_sub[temp_int] = -1;
        inp_temp.tar_sup[temp_int] = -ptr->blk(blk_number).ijkmax[temp_int];
    }
    inp_temp.volume = 1;
    for (int32_t i = 0; i < 3; i++)
        inp_temp.volume = inp_temp.volume * (abs(abs(inp_temp.sub[i]) - abs(inp_temp.sup[i])) + 1);

    for (int32_t i = 0; i < 3; i++)
    {
        inp_temp_conj.sub[i] = inp_temp.tar_sub[i];
        inp_temp_conj.sup[i] = inp_temp.tar_sup[i];
        inp_temp_conj.tar_sub[i] = inp_temp.sub[i];
        inp_temp_conj.tar_sup[i] = inp_temp.sup[i];
    }
    inp_temp_conj.volume = inp_temp.volume;

    ptr->inp(blk_number, number_of_block).inner.push_back(inp_temp);
    ptr->inp(number_of_block, blk_number).inner.push_back(inp_temp_conj);
}

/*
 * @brief 根据sub sup的范围查找_inp->inner第几个面为所需要的面号
 */
int32_t find_index_inner_face(Inner *_inp, int32_t *_sub, int32_t *_sup)
{
    bool find_it = false;
    for (int32_t iface = 0; iface < _inp->inner.size(); iface++)
    {
        coordinate_inp &_temp = _inp->inner[iface];
        if (_temp.sub[0] == _sub[0] && _temp.sub[1] == _sub[1] && _temp.sub[2] == _sub[2] && _temp.sup[0] == _sup[0] && _temp.sup[1] == _sup[1] && _temp.sup[2] == _sup[2])
        {
            find_it = true;
            return iface;
        }
    }
    if (!find_it)
    {
        std::cout << "Error, cannot find the corresponding face during split the inp information ! ! !\n";
        exit(-1);
    }
    return -10086;
}
