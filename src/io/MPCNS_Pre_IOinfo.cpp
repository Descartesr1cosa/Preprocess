#include "io/MPCNS_Pre_IOinfo.h"
#include <iostream>
#include <iomanip>

// #######################################################################
//  DATE:2022/10/02
//  Remark:对于字符串不能直接输出，sizeof(string)结果不正确
// 正确的string输出二进制文件应该为：将长度和字符内容一同写入
//  string xuke="Hasaki";                    int32_t length=xuke.size()+1;
//  file.write((char*)&length,sizeof(int32_t));     file.write(xuke.c_str(),length);
// 正确的读取string方法为:先读入长度length，分配一个长度为length的字符数组，给string分配空间
// 读入字符数组，再转入string中，最后释放字符数组（注意这里字符数组到string是深拷贝，释放字符数组不会影响string）
//  string xuke;                          int32_t length;
//  file.read(&length,sizeof(int32_t));
//  char* temp=new char[length];          xuke.reserve(length);
//  file.read(temp,length);               xuke=temp;
//  delete[] temp;
//  没有那种new分配空间的类可以直接输出
//  class YASUO...
//  YASUO hasaki;
//  file.write((char*)&hasaki,sizeof(YASUO));
//  file.read((char*)&hasaki,sizeof(YASUO));
// #######################################################################

void output_split_group_info(Preprocess_Data_Structure *ptr, Preprocess_Group *group)
{
    std::ofstream file;
    // output_split_group_info_DEBUG_only(ptr, group); // DEBUG ONLY
    file.open("./OUTPUT/split_group_info.bin", std::ios_base::out | std::ios::binary);
    if (!file)
    {
        std::string command;
        command = "mkdir OUTPUT";
        std::system(command.c_str());
    }
    file.close();
    file.open("./OUTPUT/split_group_info.bin", std::ios_base::out | std::ios::binary);
    //=============================================================================================
    //  输出ptr中的信息，包含了剖分后的网格块自身与连接关系的全部信息
    //-------------------------------------------------------------------------
    // 输出blk信息
    int32_t length; // 用来计量string类型的长度
    int32_t nblock = ptr->blk.Getsize1();
    // 一共有多少个
    file.write((char *)&nblock, sizeof(int32_t));
    for (int32_t iiblock = 0; iiblock < nblock; iiblock++)
    {
        file.write((char *)&ptr->blk(iiblock).ijkmax[0], sizeof(int32_t));
        file.write((char *)&ptr->blk(iiblock).ijkmax[1], sizeof(int32_t));
        file.write((char *)&ptr->blk(iiblock).ijkmax[2], sizeof(int32_t));
        //---------------------------------------------------------------------
        // bound
        int32_t nbound = ptr->blk(iiblock).bound.size();
        file.write((char *)&nbound, sizeof(int32_t));
        for (int32_t iibound = 0; iibound < nbound; iibound++)
            file.write((char *)&ptr->blk(iiblock).bound[iibound], sizeof(coordinate_phy));
        //---------------------------------------------------------------------
        // hiera
        int32_t nhiera = ptr->blk(iiblock).hiera.size();
        file.write((char *)&nhiera, sizeof(int32_t));
        for (int32_t iihiera = 0; iihiera < nhiera; iihiera++)
        {
            for (int i = 0; i < 3; i++)
            {
                file.write((char *)&ptr->blk(iiblock).hiera[iihiera].sub[i], sizeof(int32_t));
                file.write((char *)&ptr->blk(iiblock).hiera[iihiera].sup[i], sizeof(int32_t));
                file.write((char *)&ptr->blk(iiblock).hiera[iihiera].oring_sub[i], sizeof(int32_t));
                file.write((char *)&ptr->blk(iiblock).hiera[iihiera].oring_sup[i], sizeof(int32_t));
                file.write((char *)&ptr->blk(iiblock).hiera[iihiera].Transform[i], sizeof(int32_t));
            }
            file.write((char *)&ptr->blk(iiblock).hiera[iihiera].oring_num, sizeof(int32_t));
            length = ptr->blk(iiblock).hiera[iihiera].blk_name.size() + 1;
            file.write((char *)&length, sizeof(int32_t));
            file.write(ptr->blk(iiblock).hiera[iihiera].blk_name.c_str(), length);
        }
    }
    //-------------------------------------------------------------------------
    // 输出inp信息
    for (int32_t iiblock = 0; iiblock < nblock; iiblock++)
    {
        for (int32_t jjblock = 0; jjblock < nblock; jjblock++)
        {
            file.write((char *)&ptr->inp(iiblock, jjblock).my_num, sizeof(int32_t));
            file.write((char *)&ptr->inp(iiblock, jjblock).tar_num, sizeof(int32_t));
            //-----------------------------------------------------------------
            // coordinate_inp
            int32_t ninner = ptr->inp(iiblock, jjblock).inner.size();
            file.write((char *)&ninner, sizeof(int32_t));
            for (int32_t iiinner = 0; iiinner < ninner; iiinner++)
                file.write((char *)&ptr->inp(iiblock, jjblock).inner[iiinner], sizeof(coordinate_inp));
        }
    }
    //-------------------------------------------------------------------------
    // 输出边界名称与标号的对应：phy_index,phy_name
    // 一共有多少个
    int32_t boundary_name_number = ptr->phy_index.size();
    file.write((char *)&boundary_name_number, sizeof(int32_t));
    // 依次输出map中的元素
    for (std::map<std::string, int32_t>::iterator i = ptr->phy_index.begin(); i != ptr->phy_index.end(); i++)
    {
        length = i->first.size() + 1;
        file.write((char *)&length, sizeof(int32_t));
        file.write(i->first.c_str(), length);

        file.write((char *)&(i->second), sizeof(int32_t));
    }
    // 输出标号与边界名称的对应：phy_name
    boundary_name_number = ptr->phy_name.size();
    file.write((char *)&boundary_name_number, sizeof(int32_t));
    // 依次输出phy_name中的元素
    for (int32_t i = 0; i < boundary_name_number; i++)
    {
        length = ptr->phy_name[i].size() + 1;
        file.write((char *)&length, sizeof(int32_t));
        file.write(ptr->phy_name[i].c_str(), length);
    }
    //=============================================================================================
    //  输出ptr中的信息，包含了组合分配信息
    //-------------------------------------------------------------------------
    // block_proc_index
    int32_t mysize = group->block_proc_index.Getsize1();
    file.write((char *)&mysize, sizeof(int32_t));
    for (int32_t i = 0; i < mysize; i++)
        file.write((char *)&group->block_proc_index(i), sizeof(int32_t));
    //-------------------------------------------------------------------------
    // block_proc_number
    mysize = group->block_proc_number.Getsize1();
    file.write((char *)&mysize, sizeof(int32_t));
    for (int32_t i = 0; i < mysize; i++)
        file.write((char *)&group->block_proc_number(i), sizeof(int32_t));
    //-------------------------------------------------------------------------
    // proc_num
    mysize = group->proc_num.Getsize1();
    file.write((char *)&mysize, sizeof(int32_t));
    for (int32_t i = 0; i < mysize; i++)
        file.write((char *)&group->proc_num(i), sizeof(int32_t));
    //-------------------------------------------------------------------------
    // pro_grid_number
    mysize = group->pro_grid_number.Getsize1();
    file.write((char *)&mysize, sizeof(int32_t));
    for (int32_t i = 0; i < mysize; i++)
        file.write((char *)&group->pro_grid_number(i), sizeof(int32_t));
    //-------------------------------------------------------------------------
    // pro_block_index
    mysize = group->pro_block_index.Getsize1();
    int32_t mysize1;
    file.write((char *)&mysize, sizeof(int32_t));
    for (int32_t i = 0; i < mysize; i++)
    {
        mysize1 = group->pro_block_index(i).size();
        file.write((char *)&mysize1, sizeof(int32_t));
        for (int32_t j = 0; j < mysize1; j++)
            file.write((char *)&group->pro_block_index(i)[j], sizeof(int32_t));
    }
    //=============================================================================================
    file.close();
}

void input_split_group_info(Preprocess_Data_Structure *ptr, Preprocess_Group *group)
{
    std::ifstream file;

    file.open("./OUTPUT/split_group_info.bin", std::ios_base::in | std::ios::binary);
    if (!file)
    {
        file.close();
        std::cout << "The file of 'split_group_info.bin' is not exist! ! ! Please resplit.\n";
        exit(-1);
    }
    //=============================================================================================
    //  读入ptr中的信息，包含了剖分后的网格块自身与连接关系的全部信息
    //-------------------------------------------------------------------------
    // 读入blk信息
    int32_t length; // 用来计量string类型的长度
    char *temp;     // 临时存储string
    int32_t nblock;
    // 一共有多少个
    file.read((char *)&nblock, sizeof(int32_t));
    ptr->blk.SetSize(nblock);
    for (int32_t iiblock = 0; iiblock < nblock; iiblock++)
    {
        file.read((char *)&ptr->blk(iiblock).ijkmax[0], sizeof(int32_t));
        file.read((char *)&ptr->blk(iiblock).ijkmax[1], sizeof(int32_t));
        file.read((char *)&ptr->blk(iiblock).ijkmax[2], sizeof(int32_t));
        //---------------------------------------------------------------------
        // bound
        int32_t nbound;
        file.read((char *)&nbound, sizeof(int32_t));
        for (int32_t iibound = 0; iibound < nbound; iibound++)
        {
            coordinate_phy temp_phy;
            file.read((char *)&temp_phy, sizeof(coordinate_phy));
            ptr->blk(iiblock).bound.push_back(temp_phy);
        }

        //---------------------------------------------------------------------
        // hiera
        int32_t nhiera;
        file.read((char *)&nhiera, sizeof(int32_t));
        for (int32_t iihiera = 0; iihiera < nhiera; iihiera++)
        {
            coordinate_blk temp_blk;
            for (int i = 0; i < 3; i++)
            {
                file.read((char *)&temp_blk.sub[i], sizeof(int32_t));
                file.read((char *)&temp_blk.sup[i], sizeof(int32_t));
                file.read((char *)&temp_blk.oring_sub[i], sizeof(int32_t));
                file.read((char *)&temp_blk.oring_sup[i], sizeof(int32_t));
                file.read((char *)&temp_blk.Transform[i], sizeof(int32_t));
            }
            file.read((char *)&temp_blk.oring_num, sizeof(int32_t));

            file.read((char *)&length, sizeof(int32_t));
            temp = new char[length];
            file.read(temp, length);
            temp_blk.blk_name.reserve(length);
            temp_blk.blk_name = temp;
            delete[] temp;
            ptr->blk(iiblock).hiera.push_back(temp_blk);
        }
    }
    //-------------------------------------------------------------------------
    // 读入inp信息
    ptr->inp.SetSize(nblock, nblock);
    for (int32_t iiblock = 0; iiblock < nblock; iiblock++)
    {
        for (int32_t jjblock = 0; jjblock < nblock; jjblock++)
        {
            file.read((char *)&ptr->inp(iiblock, jjblock).my_num, sizeof(int32_t));
            file.read((char *)&ptr->inp(iiblock, jjblock).tar_num, sizeof(int32_t));
            //-----------------------------------------------------------------
            // coordinate_inp
            int32_t ninner;
            file.read((char *)&ninner, sizeof(int32_t));
            for (int32_t iiinner = 0; iiinner < ninner; iiinner++)
            {
                coordinate_inp temp_inp;
                file.read((char *)&temp_inp, sizeof(coordinate_inp));
                ptr->inp(iiblock, jjblock).inner.push_back(temp_inp);
            }
        }
    }
    //-------------------------------------------------------------------------
    // 读入边界名称与标号的对应：phy_index,phy_name
    // 一共有多少个
    int32_t boundary_name_number;
    file.read((char *)&boundary_name_number, sizeof(int32_t));
    // 依次读入map中的元素
    ptr->phy_index.clear();
    for (int32_t i = 0; i < boundary_name_number; i++)
    {
        std::string name;
        int32_t index;

        file.read((char *)&length, sizeof(int32_t));
        temp = new char[length];
        name.reserve(length);
        file.read(temp, length);
        name = temp;
        delete[] temp;

        file.read((char *)&index, sizeof(int32_t));

        ptr->phy_index.insert(std::pair<std::string, int32_t>(name, index));
    }
    // 输出标号与边界名称的对应：phy_name
    file.read((char *)&boundary_name_number, sizeof(int32_t));
    // 依次输出phy_name中的元素
    for (int32_t i = 0; i < boundary_name_number; i++)
    {
        std::string name;

        file.read((char *)&length, sizeof(int32_t));
        temp = new char[length];
        name.reserve(length);
        file.read(temp, length);
        name = temp;
        delete[] temp;

        ptr->phy_name.push_back(name);
    }
    //=============================================================================================
    //  读入group中的信息，包含了组合分配信息
    //-------------------------------------------------------------------------
    // block_proc_index
    int32_t mysize;
    file.read((char *)&mysize, sizeof(int32_t));
    group->block_proc_index.SetSize(mysize);
    for (int32_t i = 0; i < mysize; i++)
        file.read((char *)&group->block_proc_index(i), sizeof(int32_t));
    //-------------------------------------------------------------------------
    // block_proc_number
    file.read((char *)&mysize, sizeof(int32_t));
    group->block_proc_number.SetSize(mysize);
    for (int32_t i = 0; i < mysize; i++)
        file.read((char *)&group->block_proc_number(i), sizeof(int32_t));
    //-------------------------------------------------------------------------
    // proc_num
    file.read((char *)&mysize, sizeof(int32_t));
    group->proc_num.SetSize(mysize);
    for (int32_t i = 0; i < mysize; i++)
        file.read((char *)&group->proc_num(i), sizeof(int32_t));
    //-------------------------------------------------------------------------
    // pro_grid_number
    file.read((char *)&mysize, sizeof(int32_t));
    group->pro_grid_number.SetSize(mysize);
    for (int32_t i = 0; i < mysize; i++)
        file.read((char *)&group->pro_grid_number(i), sizeof(int32_t));
    //-------------------------------------------------------------------------
    // pro_block_index
    file.read((char *)&mysize, sizeof(int32_t));
    group->pro_block_index.SetSize(mysize);
    int32_t mysize1;
    for (int32_t i = 0; i < mysize; i++)
    {
        file.read((char *)&mysize1, sizeof(int32_t));

        for (int32_t j = 0; j < mysize1; j++)
        {
            int32_t temp_int;
            file.read((char *)&temp_int, sizeof(int32_t));
            group->pro_block_index(i).push_back(temp_int);
        }
    }
    //=============================================================================================
    file.close();
    // output_split_group_info_DEBUG_only(ptr, group); // DEBUG ONLY
}

void output_split_group_info_DEBUG_only(Preprocess_Data_Structure *ptr, Preprocess_Group *group)
{
    //=============================================================================================
    //  输出ptr中的信息，包含了剖分后的网格块自身与连接关系的全部信息
    //-------------------------------------------------------------------------
    // 输出blk信息
    int32_t length; // 用来计量string类型的长度
    int32_t nblock = ptr->blk.Getsize1();
    // 一共有多少个
    std::cout << "nblock=" << nblock << std::endl;
    for (int32_t iiblock = 0; iiblock < nblock; iiblock++)
    {
        std::cout << "Block\t" << iiblock << std::endl;
        std::cout << ptr->blk(iiblock).ijkmax[0] << "\t" << ptr->blk(iiblock).ijkmax[1]
                  << "\t" << ptr->blk(iiblock).ijkmax[2] << std::endl;
        //---------------------------------------------------------------------
        // bound
        int32_t nbound = ptr->blk(iiblock).bound.size();
        std::cout << "Block-bound\t" << nbound << std::endl;

        for (int32_t iibound = 0; iibound < nbound; iibound++)
        {
            std::cout << "\tbound" << iibound << std::endl;
            std::cout << "\t" << ptr->blk(iiblock).bound[iibound].sub[0] << "\t" << ptr->blk(iiblock).bound[iibound].sub[1]
                      << "\t" << ptr->blk(iiblock).bound[iibound].sub[2] << std::endl;
            std::cout << "\t" << ptr->blk(iiblock).bound[iibound].sup[0] << "\t" << ptr->blk(iiblock).bound[iibound].sup[1]
                      << "\t" << ptr->blk(iiblock).bound[iibound].sup[2] << std::endl;
            std::cout << "\t" << ptr->blk(iiblock).bound[iibound].phy_kind << std::endl;
        }
        //---------------------------------------------------------------------
        // hiera
        int32_t nhiera = ptr->blk(iiblock).hiera.size();
        std::cout << "hiera\t" << nhiera << std::endl;
        for (int32_t iihiera = 0; iihiera < nhiera; iihiera++)
        {
            std::cout << "\t-hiera\t" << iihiera << std::endl;
            for (int i = 0; i < 3; i++)
            {
                std::cout << "\t" << ptr->blk(iiblock).hiera[iihiera].sub[i] << "\t" << ptr->blk(iiblock).hiera[iihiera].sup[i]
                          << "\t" << ptr->blk(iiblock).hiera[iihiera].oring_sub[i] << "\t" << ptr->blk(iiblock).hiera[iihiera].oring_sup[i]
                          << "\t" << ptr->blk(iiblock).hiera[iihiera].Transform[i] << std::endl;
            }
            std::cout << "\t\t->The blk_name is " << ptr->blk(iiblock).hiera[iihiera].blk_name << std::endl;
        }
    }
    std::cout << "--------------------------------------------------------------------------";
    //-------------------------------------------------------------------------
    // 输出inp信息
    for (int32_t iiblock = 0; iiblock < nblock; iiblock++)
    {
        for (int32_t jjblock = 0; jjblock < nblock; jjblock++)
        {
            std::cout << "\ti= " << iiblock << "\tj= " << jjblock << std::endl;
            int32_t ninner = ptr->inp(iiblock, jjblock).inner.size();
            std::cout << "\t\t ninner=\t" << ninner << std::endl;
            //-----------------------------------------------------------------
            // coordinate_inp
            for (int32_t iiinner = 0; iiinner < ninner; iiinner++)
            {
                std::cout << ptr->inp(iiblock, jjblock).inner[iiinner].sub[0] << "\t"
                          << ptr->inp(iiblock, jjblock).inner[iiinner].sub[1] << "\t"
                          << ptr->inp(iiblock, jjblock).inner[iiinner].sub[2] << "\t"
                          << ptr->inp(iiblock, jjblock).inner[iiinner].sup[0] << "\t"
                          << ptr->inp(iiblock, jjblock).inner[iiinner].sup[1] << "\t"
                          << ptr->inp(iiblock, jjblock).inner[iiinner].sup[2] << "\n"
                          << ptr->inp(iiblock, jjblock).inner[iiinner].tar_sub[0] << "\t"
                          << ptr->inp(iiblock, jjblock).inner[iiinner].tar_sub[1] << "\t"
                          << ptr->inp(iiblock, jjblock).inner[iiinner].tar_sub[2] << "\t"
                          << ptr->inp(iiblock, jjblock).inner[iiinner].tar_sup[0] << "\t"
                          << ptr->inp(iiblock, jjblock).inner[iiinner].tar_sup[1] << "\t"
                          << ptr->inp(iiblock, jjblock).inner[iiinner].tar_sup[2] << "\n"
                          << ptr->inp(iiblock, jjblock).inner[iiinner].is_parallel
                          << ptr->inp(iiblock, jjblock).inner[iiinner].volume << "\n\n";
            }
        }
    }
    //-------------------------------------------------------------------------
    // 输出边界名称与标号的对应：phy_index,phy_name
    // 一共有多少个
    std::cout << "--------------------------------------------------------------------------";
    int32_t boundary_name_number = ptr->phy_index.size();
    std::cout << "number of phy_index\t" << boundary_name_number << std::endl;
    // 依次输出map中的元素
    for (std::map<std::string, int32_t>::iterator i = ptr->phy_index.begin(); i != ptr->phy_index.end(); i++)
    {
        std::cout << "\t->The pair " << i->first << "\t" << i->second << std::endl;
    }
    // 输出标号与边界名称的对应：phy_name
    boundary_name_number = ptr->phy_name.size();
    std::cout << "number of phy_name\t" << boundary_name_number << std::endl;
    // 依次输出phy_name中的元素
    for (int32_t i = 0; i < boundary_name_number; i++)
    {
        std::cout << "\t->The vector " << i << "\t" << ptr->phy_name[i] << std::endl;
    }
    //=============================================================================================
    // //  输出ptr中的信息，包含了组合分配信息
    // //-------------------------------------------------------------------------
    // // block_proc_index
    // int32_t mysize = group->block_proc_index.Getsize1();
    // file.write((char *)&mysize, sizeof(int32_t));
    // for (int32_t i = 0; i < mysize; i++)
    //     file.write((char *)&group->block_proc_index(i), sizeof(int32_t));
    // //-------------------------------------------------------------------------
    // // block_proc_number
    // mysize = group->block_proc_number.Getsize1();
    // file.write((char *)&mysize, sizeof(int32_t));
    // for (int32_t i = 0; i < mysize; i++)
    //     file.write((char *)&group->block_proc_number(i), sizeof(int32_t));
    // //-------------------------------------------------------------------------
    // // proc_num
    // mysize = group->proc_num.Getsize1();
    // file.write((char *)&mysize, sizeof(int32_t));
    // for (int32_t i = 0; i < mysize; i++)
    //     file.write((char *)&group->proc_num(i), sizeof(int32_t));
    // //-------------------------------------------------------------------------
    // // pro_grid_number
    // mysize = group->pro_grid_number.Getsize1();
    // file.write((char *)&mysize, sizeof(int32_t));
    // for (int32_t i = 0; i < mysize; i++)
    //     file.write((char *)&group->pro_grid_number(i), sizeof(int32_t));
    // //-------------------------------------------------------------------------
    // // pro_block_index
    // mysize = group->pro_block_index.Getsize1();
    // int32_t mysize1;
    // file.write((char *)&mysize, sizeof(int32_t));
    // for (int32_t i = 0; i < mysize; i++)
    // {
    //     mysize1 = group->pro_block_index(i).size();
    //     file.write((char *)&mysize1, sizeof(int32_t));
    //     for (int32_t j = 0; j < mysize1; j++)
    //         file.write((char *)&group->pro_block_index(i)[j], sizeof(int32_t));
    // }
    // //=============================================================================================
}
