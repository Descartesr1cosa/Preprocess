#include "report/MPCNS_Pre_TopologyReport.h"

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>

namespace
{
enum FaceId
{
    XI_MINUS = 0,
    XI_PLUS = 1,
    ETA_MINUS = 2,
    ETA_PLUS = 3,
    ZETA_MINUS = 4,
    ZETA_PLUS = 5,
    FACE_UNKNOWN = 6
};

struct Patch
{
    int id;
    std::string kind;
    int source_block;
    int source_face;
    int source_box[3][2];
    int target_block;
    int target_face;
    int target_box[3][2];
    int boundary_kind;
    std::string boundary_name;
    int source_myid;
    int source_local;
    int target_myid;
    int target_local;
    int reverse_patch;
    int interface_id;
    bool periodic;
};

struct InterfaceInfo
{
    int id;
    int patch_a;
    int patch_b;
    int block_a;
    int block_b;
    int myid_a;
    int myid_b;
    bool inter_rank;
    bool periodic;
};

struct WarningInfo
{
    std::string kind;
    std::string message;
    std::string link;
};

struct FaceReport
{
    int block;
    int face;
    std::vector<int> patch_ids;
    std::vector<int> cuts_u;
    std::vector<int> cuts_v;
    int gaps;
    int overlaps;
    std::string markdown_link;
};

struct ReportData
{
    Preprocess_Data_Structure *ptr;
    Preprocess_Group *group;
    Param *param;
    int nblock;
    int nproc;
    int block_width;
    int rank_width;
    int patch_width;
    int interface_width;
    bool ascii_layout;
    std::vector<Patch> patches;
    std::vector<InterfaceInfo> interfaces;
    std::vector<WarningInfo> warnings;
    std::vector<std::vector<int> > block_patch_ids;
    std::vector<std::vector<int> > rank_patch_ids;
    std::vector<FaceReport> faces;
};

std::string ToString(int value)
{
    std::ostringstream stream;
    stream << value;
    return stream.str();
}

int WidthForMax(int max_value)
{
    int width = 1;
    while (max_value >= 10)
    {
        max_value /= 10;
        width++;
    }
    return width;
}

std::string Pad(int value, int width)
{
    std::ostringstream stream;
    stream << std::setw(width) << std::setfill('0') << value;
    return stream.str();
}

std::string MdEscape(const std::string &text)
{
    std::string out;
    for (std::size_t i = 0; i < text.size(); ++i)
    {
        if (text[i] == '|')
            out += "\\|";
        else
            out += text[i];
    }
    return out;
}

std::string JsonEscape(const std::string &text)
{
    std::string out;
    for (std::size_t i = 0; i < text.size(); ++i)
    {
        const char c = text[i];
        if (c == '"' || c == '\\')
            out += '\\';
        if (c == '\n')
            out += "\\n";
        else
            out += c;
    }
    return out;
}

void EnsureDirectory(const std::string &path)
{
    if (mkdir(path.c_str(), 0755) != 0 && errno != EEXIST)
    {
        std::cout << "#Fatal Error: Cannot create directory:\t" << path << std::endl;
        exit(-1);
    }
}

void EnsureReportDirectories()
{
    EnsureDirectory("OUTPUT");
    EnsureDirectory("OUTPUT/topology_report");
    EnsureDirectory("OUTPUT/topology_report/ranks");
    EnsureDirectory("OUTPUT/topology_report/blocks");
    EnsureDirectory("OUTPUT/topology_report/faces");
    EnsureDirectory("OUTPUT/topology_report/patches");
    EnsureDirectory("OUTPUT/topology_report/interfaces");
    EnsureDirectory("OUTPUT/topology_report/csv");
}

std::string FaceName(int face)
{
    switch (face)
    {
    case XI_MINUS:
        return "xi_minus";
    case XI_PLUS:
        return "xi_plus";
    case ETA_MINUS:
        return "eta_minus";
    case ETA_PLUS:
        return "eta_plus";
    case ZETA_MINUS:
        return "zeta_minus";
    case ZETA_PLUS:
        return "zeta_plus";
    default:
        return "unknown";
    }
}

std::string FaceLabel(int face)
{
    switch (face)
    {
    case XI_MINUS:
        return "xi-";
    case XI_PLUS:
        return "xi+";
    case ETA_MINUS:
        return "eta-";
    case ETA_PLUS:
        return "eta+";
    case ZETA_MINUS:
        return "zeta-";
    case ZETA_PLUS:
        return "zeta+";
    default:
        return "unknown";
    }
}

std::string FaceAxes(int face)
{
    switch (face)
    {
    case XI_MINUS:
    case XI_PLUS:
        return "eta-zeta";
    case ETA_MINUS:
    case ETA_PLUS:
        return "xi-zeta";
    case ZETA_MINUS:
    case ZETA_PLUS:
        return "xi-eta";
    default:
        return "unknown";
    }
}

int NormalAxis(int face)
{
    if (face == XI_MINUS || face == XI_PLUS)
        return 0;
    if (face == ETA_MINUS || face == ETA_PLUS)
        return 1;
    if (face == ZETA_MINUS || face == ZETA_PLUS)
        return 2;
    return -1;
}

void TangentialAxes(int face, int &u_axis, int &v_axis)
{
    if (face == XI_MINUS || face == XI_PLUS)
    {
        u_axis = 1;
        v_axis = 2;
    }
    else if (face == ETA_MINUS || face == ETA_PLUS)
    {
        u_axis = 0;
        v_axis = 2;
    }
    else
    {
        u_axis = 0;
        v_axis = 1;
    }
}

int DetectFace(const int box[3][2], const int size[3])
{
    for (int d = 0; d < 3; ++d)
    {
        const int lo = std::min(std::abs(box[d][0]), std::abs(box[d][1]));
        const int hi = std::max(std::abs(box[d][0]), std::abs(box[d][1]));
        if (lo == hi && lo == 1)
            return d == 0 ? XI_MINUS : (d == 1 ? ETA_MINUS : ZETA_MINUS);
        if (lo == hi && lo == size[d])
            return d == 0 ? XI_PLUS : (d == 1 ? ETA_PLUS : ZETA_PLUS);
    }
    return FACE_UNKNOWN;
}

std::string BoxText(const int box[3][2])
{
    std::ostringstream s;
    s << "[" << box[0][0] << ":" << box[0][1] << ", "
      << box[1][0] << ":" << box[1][1] << ", "
      << box[2][0] << ":" << box[2][1] << "]";
    return s.str();
}

int BoxExtent(const int box[3][2], int d)
{
    return std::abs(std::abs(box[d][1]) - std::abs(box[d][0])) + 1;
}

int PatchSize(const Patch &patch)
{
    return BoxExtent(patch.source_box, 0) * BoxExtent(patch.source_box, 1) * BoxExtent(patch.source_box, 2);
}

std::string BlockPath(const ReportData &data, int block, const std::string &prefix)
{
    return prefix + "blocks/block_" + Pad(block, data.block_width) + ".md";
}

std::string RankPath(const ReportData &data, int rank, const std::string &prefix)
{
    return prefix + "ranks/myid_" + Pad(rank, data.rank_width) + ".md";
}

std::string FacePath(const ReportData &data, int block, int face, const std::string &prefix)
{
    return prefix + "faces/block_" + Pad(block, data.block_width) + "_" + FaceName(face) + ".md";
}

std::string PatchPath(const ReportData &data, int patch, const std::string &prefix)
{
    return prefix + "patches/patch_" + Pad(patch, data.patch_width) + ".md";
}

std::string InterfacePath(const ReportData &data, int interface_id, const std::string &prefix)
{
    return prefix + "interfaces/interface_" + Pad(interface_id, data.interface_width) + ".md";
}

std::string OwnerLabel(const ReportData &data, const Patch &patch)
{
    if (patch.kind == "physical")
        return patch.boundary_name.empty() ? "BOUNDARY" : patch.boundary_name;
    return "B" + ToString(patch.target_block) + " " + FaceLabel(patch.target_face);
}

bool SameBox(const int a[3][2], const int b[3][2])
{
    for (int d = 0; d < 3; ++d)
        for (int e = 0; e < 2; ++e)
            if (a[d][e] != b[d][e])
                return false;
    return true;
}

void AddWarning(ReportData &data, const std::string &kind, const std::string &message, const std::string &link)
{
    WarningInfo warning;
    warning.kind = kind;
    warning.message = message;
    warning.link = link;
    data.warnings.push_back(warning);
}

void CopyBox(int target[3][2], const int sub[3], const int sup[3])
{
    for (int d = 0; d < 3; ++d)
    {
        target[d][0] = sub[d];
        target[d][1] = sup[d];
    }
}

std::string BoundaryName(Preprocess_Data_Structure *ptr, int kind)
{
    if (kind >= 0 && kind < static_cast<int>(ptr->phy_name.size()))
        return ptr->phy_name[kind];
    return "boundary_" + ToString(kind);
}

void CollectPatches(ReportData &data)
{
    data.block_patch_ids.assign(data.nblock, std::vector<int>());
    data.rank_patch_ids.assign(data.nproc, std::vector<int>());

    for (int iblock = 0; iblock < data.nblock; ++iblock)
    {
        block_info &block = data.ptr->blk(iblock);
        for (std::size_t ibound = 0; ibound < block.bound.size(); ++ibound)
        {
            Patch patch;
            patch.id = static_cast<int>(data.patches.size());
            patch.kind = "physical";
            patch.source_block = iblock;
            CopyBox(patch.source_box, block.bound[ibound].sub, block.bound[ibound].sup);
            patch.source_face = DetectFace(patch.source_box, block.ijkmax);
            patch.target_block = -1;
            patch.target_face = FACE_UNKNOWN;
            for (int d = 0; d < 3; ++d)
            {
                patch.target_box[d][0] = 0;
                patch.target_box[d][1] = 0;
            }
            patch.boundary_kind = block.bound[ibound].phy_kind;
            patch.boundary_name = BoundaryName(data.ptr, patch.boundary_kind);
            patch.source_myid = data.group->block_proc_number(iblock);
            patch.source_local = data.group->block_proc_index(iblock);
            patch.target_myid = -1;
            patch.target_local = -1;
            patch.reverse_patch = -1;
            patch.interface_id = -1;
            patch.periodic = false;
            data.block_patch_ids[iblock].push_back(patch.id);
            if (patch.source_myid >= 0 && patch.source_myid < data.nproc)
                data.rank_patch_ids[patch.source_myid].push_back(patch.id);
            data.patches.push_back(patch);
        }
    }

    for (int iblock = 0; iblock < data.nblock; ++iblock)
    {
        for (int jblock = 0; jblock < data.nblock; ++jblock)
        {
            Inner &inner = data.ptr->inp(iblock, jblock);
            for (std::size_t iface = 0; iface < inner.inner.size(); ++iface)
            {
                coordinate_inp &edge = inner.inner[iface];
                Patch patch;
                patch.id = static_cast<int>(data.patches.size());
                patch.kind = edge.is_parallel ? "periodic" : "inner";
                patch.source_block = iblock;
                CopyBox(patch.source_box, edge.sub, edge.sup);
                patch.source_face = DetectFace(patch.source_box, data.ptr->blk(iblock).ijkmax);
                patch.target_block = jblock;
                CopyBox(patch.target_box, edge.tar_sub, edge.tar_sup);
                patch.target_face = DetectFace(patch.target_box, data.ptr->blk(jblock).ijkmax);
                patch.boundary_kind = -1;
                patch.boundary_name = "";
                patch.source_myid = data.group->block_proc_number(iblock);
                patch.source_local = data.group->block_proc_index(iblock);
                patch.target_myid = data.group->block_proc_number(jblock);
                patch.target_local = data.group->block_proc_index(jblock);
                patch.reverse_patch = -1;
                patch.interface_id = -1;
                patch.periodic = edge.is_parallel != 0;
                data.block_patch_ids[iblock].push_back(patch.id);
                if (patch.source_myid >= 0 && patch.source_myid < data.nproc)
                    data.rank_patch_ids[patch.source_myid].push_back(patch.id);
                data.patches.push_back(patch);
            }
        }
    }
}

void DiagnosePatchRanges(ReportData &data)
{
    for (std::size_t i = 0; i < data.patches.size(); ++i)
    {
        Patch &patch = data.patches[i];
        if (patch.source_face == FACE_UNKNOWN)
            AddWarning(data, "face-location", "Patch " + ToString(patch.id) + " is not located on a logical source face.", PatchPath(data, patch.id, ""));
        for (int d = 0; d < 3; ++d)
        {
            const int lo = std::min(std::abs(patch.source_box[d][0]), std::abs(patch.source_box[d][1]));
            const int hi = std::max(std::abs(patch.source_box[d][0]), std::abs(patch.source_box[d][1]));
            if (lo < 1 || hi > data.ptr->blk(patch.source_block).ijkmax[d])
                AddWarning(data, "index-range", "Patch " + ToString(patch.id) + " source range is outside block " + ToString(patch.source_block) + ".", PatchPath(data, patch.id, ""));
        }
        if (patch.target_block >= 0)
        {
            if (patch.target_face == FACE_UNKNOWN)
                AddWarning(data, "face-location", "Patch " + ToString(patch.id) + " target is not located on a logical face.", PatchPath(data, patch.id, ""));
            for (int d = 0; d < 3; ++d)
            {
                const int lo = std::min(std::abs(patch.target_box[d][0]), std::abs(patch.target_box[d][1]));
                const int hi = std::max(std::abs(patch.target_box[d][0]), std::abs(patch.target_box[d][1]));
                if (lo < 1 || hi > data.ptr->blk(patch.target_block).ijkmax[d])
                    AddWarning(data, "index-range", "Patch " + ToString(patch.id) + " target range is outside block " + ToString(patch.target_block) + ".", PatchPath(data, patch.id, ""));
            }
            if (PatchSize(patch) != BoxExtent(patch.target_box, 0) * BoxExtent(patch.target_box, 1) * BoxExtent(patch.target_box, 2))
                AddWarning(data, "patch-size", "Patch " + ToString(patch.id) + " source and target patch sizes differ.", PatchPath(data, patch.id, ""));
        }
    }
}

void LinkReversePatchesAndInterfaces(ReportData &data)
{
    for (std::size_t i = 0; i < data.patches.size(); ++i)
    {
        Patch &patch = data.patches[i];
        if (patch.target_block < 0 || patch.reverse_patch >= 0)
            continue;

        int reverse_id = -1;
        for (std::size_t j = 0; j < data.patches.size(); ++j)
        {
            Patch &candidate = data.patches[j];
            if (candidate.target_block != patch.source_block || candidate.source_block != patch.target_block)
                continue;
            if (SameBox(candidate.source_box, patch.target_box) && SameBox(candidate.target_box, patch.source_box))
            {
                reverse_id = candidate.id;
                break;
            }
        }

        patch.reverse_patch = reverse_id;
        if (reverse_id >= 0)
            data.patches[reverse_id].reverse_patch = patch.id;

        if (reverse_id < 0 || patch.id <= reverse_id)
        {
            InterfaceInfo iface;
            iface.id = static_cast<int>(data.interfaces.size());
            iface.patch_a = patch.id;
            iface.patch_b = reverse_id;
            iface.block_a = patch.source_block;
            iface.block_b = patch.target_block;
            iface.myid_a = patch.source_myid;
            iface.myid_b = patch.target_myid;
            iface.inter_rank = iface.myid_a != iface.myid_b;
            iface.periodic = patch.periodic;
            patch.interface_id = iface.id;
            if (reverse_id >= 0)
                data.patches[reverse_id].interface_id = iface.id;
            data.interfaces.push_back(iface);
        }
    }
}

std::vector<int> FacePatchIds(const ReportData &data, int block, int face)
{
    std::vector<int> ids;
    const std::vector<int> &block_patches = data.block_patch_ids[block];
    for (std::size_t i = 0; i < block_patches.size(); ++i)
    {
        int patch_id = block_patches[i];
        if (data.patches[patch_id].source_face == face)
            ids.push_back(patch_id);
    }
    return ids;
}

void PrepareFaceReports(ReportData &data)
{
    for (int block = 0; block < data.nblock; ++block)
    {
        for (int face = 0; face < 6; ++face)
        {
            FaceReport report;
            report.block = block;
            report.face = face;
            report.patch_ids = FacePatchIds(data, block, face);
            report.gaps = 0;
            report.overlaps = 0;
            report.markdown_link = FacePath(data, block, face, "");
            data.faces.push_back(report);
        }
    }
}

std::string PatchRange2D(const Patch &patch, int face)
{
    int u_axis = 0;
    int v_axis = 1;
    TangentialAxes(face, u_axis, v_axis);
    const int ulo = std::min(std::abs(patch.source_box[u_axis][0]), std::abs(patch.source_box[u_axis][1]));
    const int uhi = std::max(std::abs(patch.source_box[u_axis][0]), std::abs(patch.source_box[u_axis][1]));
    const int vlo = std::min(std::abs(patch.source_box[v_axis][0]), std::abs(patch.source_box[v_axis][1]));
    const int vhi = std::max(std::abs(patch.source_box[v_axis][0]), std::abs(patch.source_box[v_axis][1]));
    std::ostringstream s;
    s << "[" << ulo << ":" << uhi << "] x [" << vlo << ":" << vhi << "]";
    return s.str();
}

void PatchCellInterval2D(const Patch &patch, int face, int &ulo, int &uhi, int &vlo, int &vhi)
{
    int u_axis = 0;
    int v_axis = 1;
    TangentialAxes(face, u_axis, v_axis);
    ulo = std::min(std::abs(patch.source_box[u_axis][0]), std::abs(patch.source_box[u_axis][1]));
    uhi = std::max(std::abs(patch.source_box[u_axis][0]), std::abs(patch.source_box[u_axis][1]));
    vlo = std::min(std::abs(patch.source_box[v_axis][0]), std::abs(patch.source_box[v_axis][1]));
    vhi = std::max(std::abs(patch.source_box[v_axis][0]), std::abs(patch.source_box[v_axis][1]));
}

bool PatchCoversCell(const Patch &patch, int face, int ulo, int uhi, int vlo, int vhi)
{
    int pu0 = 0;
    int pu1 = 0;
    int pv0 = 0;
    int pv1 = 0;
    PatchCellInterval2D(patch, face, pu0, pu1, pv0, pv1);
    if (pu0 >= pu1 || pv0 >= pv1)
        return false;
    return pu0 <= ulo && uhi <= pu1 && pv0 <= vlo && vhi <= pv1;
}

std::string ShortOwner(const std::string &owner)
{
    if (owner.size() <= 12)
        return owner;
    return owner.substr(0, 12);
}

void AnalyzeAndWriteFace(ReportData &data, FaceReport &face_report)
{
    const int block = face_report.block;
    const int face = face_report.face;
    int u_axis = 0;
    int v_axis = 1;
    TangentialAxes(face, u_axis, v_axis);
    const int umax = data.ptr->blk(block).ijkmax[u_axis];
    const int vmax = data.ptr->blk(block).ijkmax[v_axis];

    std::set<int> cut_u;
    std::set<int> cut_v;
    cut_u.insert(1);
    cut_u.insert(umax);
    cut_v.insert(1);
    cut_v.insert(vmax);
    for (std::size_t i = 0; i < face_report.patch_ids.size(); ++i)
    {
        const Patch &patch = data.patches[face_report.patch_ids[i]];
        int ulo = 0;
        int uhi = 0;
        int vlo = 0;
        int vhi = 0;
        PatchCellInterval2D(patch, face, ulo, uhi, vlo, vhi);
        if (ulo < uhi)
        {
            cut_u.insert(ulo);
            cut_u.insert(uhi);
        }
        if (vlo < vhi)
        {
            cut_v.insert(vlo);
            cut_v.insert(vhi);
        }
    }
    face_report.cuts_u.assign(cut_u.begin(), cut_u.end());
    face_report.cuts_v.assign(cut_v.begin(), cut_v.end());

    std::ofstream file(("OUTPUT/topology_report/" + face_report.markdown_link).c_str());
    file << "# Face block " << block << " " << FaceLabel(face) << "\n\n";
    file << "[Index](../index.md) | [Block](../" << BlockPath(data, block, "") << ")\n\n";
    file << "- 2D axes: `" << FaceAxes(face) << "`\n";
    file << "- Patch count: " << face_report.patch_ids.size() << "\n";
    file << "- Source block size: " << data.ptr->blk(block).ijkmax[0] << " x "
         << data.ptr->blk(block).ijkmax[1] << " x " << data.ptr->blk(block).ijkmax[2] << "\n\n";

    file << "## Patch Table\n\n";
    file << "| Patch | Kind | 2D range | Target / Boundary | Interface |\n";
    file << "|---|---|---|---|---|\n";
    for (std::size_t i = 0; i < face_report.patch_ids.size(); ++i)
    {
        const Patch &patch = data.patches[face_report.patch_ids[i]];
        file << "| [patch " << patch.id << "](../" << PatchPath(data, patch.id, "") << ") | "
             << patch.kind << " | `" << PatchRange2D(patch, face) << "` | ";
        if (patch.target_block >= 0)
            file << "[B" << patch.target_block << "](../" << BlockPath(data, patch.target_block, "") << ") " << FaceLabel(patch.target_face);
        else
            file << MdEscape(patch.boundary_name);
        file << " | ";
        if (patch.interface_id >= 0)
            file << "[interface " << patch.interface_id << "](../" << InterfacePath(data, patch.interface_id, "") << ")";
        file << " |\n";
    }

    file << "\n## Coverage Diagnostics\n\n";
    file << "Coverage is checked on face cells, not on face nodes. A patch node range `[lo:hi]` covers the half-open cell interval `[lo,hi)`, so shared nodes or shared edges are not reported as overlaps.\n\n";
    file << "- u cuts: ";
    for (std::size_t i = 0; i < face_report.cuts_u.size(); ++i)
        file << (i ? ", " : "") << face_report.cuts_u[i];
    file << "\n- v cuts: ";
    for (std::size_t i = 0; i < face_report.cuts_v.size(); ++i)
        file << (i ? ", " : "") << face_report.cuts_v[i];
    file << "\n\n";

    const int cells_u = static_cast<int>(face_report.cuts_u.size()) - 1;
    const int cells_v = static_cast<int>(face_report.cuts_v.size()) - 1;
    std::vector<std::vector<std::string> > cell_owner(cells_v, std::vector<std::string>(cells_u, ""));

    file << "| Cell | u half-open | v half-open | Coverage | Owner |\n";
    file << "|---|---|---|---|---|\n";
    for (int j = 0; j < cells_v; ++j)
    {
        for (int i = 0; i < cells_u; ++i)
        {
            const int ulo = face_report.cuts_u[i];
            const int uhi = face_report.cuts_u[i + 1];
            const int vlo = face_report.cuts_v[j];
            const int vhi = face_report.cuts_v[j + 1];
            std::vector<int> covering;
            for (std::size_t p = 0; p < face_report.patch_ids.size(); ++p)
                if (PatchCoversCell(data.patches[face_report.patch_ids[p]], face, ulo, uhi, vlo, vhi))
                    covering.push_back(face_report.patch_ids[p]);

            std::string owner;
            if (covering.empty())
            {
                face_report.gaps++;
                owner = "GAP";
            }
            else if (covering.size() > 1)
            {
                face_report.overlaps++;
                owner = "OVERLAP";
            }
            else
            {
                owner = OwnerLabel(data, data.patches[covering[0]]);
            }
            cell_owner[j][i] = ShortOwner(owner);
            file << "| (" << i << "," << j << ") | [" << ulo << "," << uhi << ") | ["
                 << vlo << "," << vhi << ") | " << covering.size() << " | " << MdEscape(owner) << " |\n";
        }
    }

    if (face_report.gaps > 0)
        AddWarning(data, "face-gap", "Face block " + ToString(block) + " " + FaceLabel(face) + " has " + ToString(face_report.gaps) + " uncovered cells.", face_report.markdown_link);
    if (face_report.overlaps > 0)
        AddWarning(data, "face-overlap", "Face block " + ToString(block) + " " + FaceLabel(face) + " has " + ToString(face_report.overlaps) + " overlapped cells.", face_report.markdown_link);

    file << "\n- Gap cells: " << face_report.gaps << "\n";
    file << "- Overlap cells: " << face_report.overlaps << "\n\n";

    file << "## ASCII Layout\n\n";
    const int max_ascii_cells = 400;
    if (!data.ascii_layout)
    {
        file << "ASCII layout disabled by `topology_report_ascii_layout`.\n";
    }
    else if (cells_u * cells_v > max_ascii_cells)
    {
        file << "ASCII layout skipped because the cut grid has " << cells_u * cells_v
             << " cells. The cell table above is still complete.\n";
    }
    else
    {
        file << "```text\n";
        for (int j = cells_v - 1; j >= 0; --j)
        {
            for (int i = 0; i < cells_u; ++i)
                file << std::setw(12) << std::left << cell_owner[j][i] << " ";
            file << "\n";
        }
        file << "```\n";
    }
}

void WriteCsv(ReportData &data)
{
    {
        std::ofstream file("OUTPUT/topology_report/csv/block_to_myid.csv");
        file << "global_block,myid,local_block,nx,ny,nz\n";
        for (int block = 0; block < data.nblock; ++block)
            file << block << "," << data.group->block_proc_number(block) << ","
                 << data.group->block_proc_index(block) << ","
                 << data.ptr->blk(block).ijkmax[0] << "," << data.ptr->blk(block).ijkmax[1] << ","
                 << data.ptr->blk(block).ijkmax[2] << "\n";
    }
    {
        std::ofstream file("OUTPUT/topology_report/csv/myid_to_blocks.csv");
        file << "myid,local_block,global_block,nx,ny,nz\n";
        for (int rank = 0; rank < data.nproc; ++rank)
            for (std::size_t i = 0; i < data.group->pro_block_index(rank).size(); ++i)
            {
                const int block = data.group->pro_block_index(rank)[i];
                file << rank << "," << i << "," << block << ","
                     << data.ptr->blk(block).ijkmax[0] << "," << data.ptr->blk(block).ijkmax[1] << ","
                     << data.ptr->blk(block).ijkmax[2] << "\n";
            }
    }
    {
        std::ofstream file("OUTPUT/topology_report/csv/inter_rank_interfaces.csv");
        file << "interface,patch_a,patch_b,block_a,block_b,myid_a,myid_b,periodic\n";
        for (std::size_t i = 0; i < data.interfaces.size(); ++i)
        {
            const InterfaceInfo &iface = data.interfaces[i];
            if (!iface.inter_rank)
                continue;
            file << iface.id << "," << iface.patch_a << "," << iface.patch_b << ","
                 << iface.block_a << "," << iface.block_b << "," << iface.myid_a << ","
                 << iface.myid_b << "," << (iface.periodic ? 1 : 0) << "\n";
        }
    }
}

void WriteIndexAndSummary(ReportData &data)
{
    int physical = 0;
    int inter_rank = 0;
    for (std::size_t i = 0; i < data.patches.size(); ++i)
        if (data.patches[i].kind == "physical")
            physical++;
    for (std::size_t i = 0; i < data.interfaces.size(); ++i)
        if (data.interfaces[i].inter_rank)
            inter_rank++;

    {
        std::ofstream file("OUTPUT/topology_report/index.md");
        file << "# Topology Report\n\n";
        file << "- [Summary](summary.md)\n";
        file << "- [Warnings](warnings.md)\n";
        file << "- [JSON](topology.json)\n";
        file << "- [CSV: block_to_myid](csv/block_to_myid.csv)\n";
        file << "- [CSV: myid_to_blocks](csv/myid_to_blocks.csv)\n";
        file << "- [CSV: inter_rank_interfaces](csv/inter_rank_interfaces.csv)\n\n";
        file << "## Ranks\n\n";
        for (int rank = 0; rank < data.nproc; ++rank)
            file << "- [myid " << rank << "](" << RankPath(data, rank, "") << ") blocks="
                 << data.group->pro_block_index(rank).size() << "\n";
    }

    {
        std::ofstream file("OUTPUT/topology_report/summary.md");
        file << "# Summary\n\n";
        file << "| Metric | Value |\n|---|---:|\n";
        file << "| Blocks | " << data.nblock << " |\n";
        file << "| Myids | " << data.nproc << " |\n";
        file << "| Patches | " << data.patches.size() << " |\n";
        file << "| Physical boundary patches | " << physical << " |\n";
        file << "| Interfaces | " << data.interfaces.size() << " |\n";
        file << "| Inter-rank interfaces | " << inter_rank << " |\n";
        file << "| Warnings | " << data.warnings.size() << " |\n";
    }
}

void WriteWarnings(ReportData &data)
{
    std::ofstream file("OUTPUT/topology_report/warnings.md");
    file << "# Warnings\n\n";
    if (data.warnings.empty())
    {
        file << "No topology diagnostics were reported.\n";
        return;
    }
    std::map<std::string, int> counts;
    for (std::size_t i = 0; i < data.warnings.size(); ++i)
        counts[data.warnings[i].kind]++;
    file << "## Counts\n\n";
    for (std::map<std::string, int>::iterator it = counts.begin(); it != counts.end(); ++it)
        file << "- " << it->first << ": " << it->second << "\n";
    file << "\n## Details\n\n";
    file << "| Kind | Message | Link |\n|---|---|---|\n";
    for (std::size_t i = 0; i < data.warnings.size(); ++i)
    {
        const WarningInfo &warning = data.warnings[i];
        file << "| " << warning.kind << " | " << MdEscape(warning.message) << " | ";
        if (!warning.link.empty())
            file << "[open](" << warning.link << ")";
        file << " |\n";
    }
}

void WriteRankPages(ReportData &data)
{
    for (int rank = 0; rank < data.nproc; ++rank)
    {
        std::ofstream file(("OUTPUT/topology_report/" + RankPath(data, rank, "")).c_str());
        file << "# myid " << rank << "\n\n";
        file << "[Index](../index.md)\n\n";
        file << "## Blocks\n\n";
        file << "| Global block | Local block | Size | Link |\n|---:|---:|---|---|\n";
        for (std::size_t i = 0; i < data.group->pro_block_index(rank).size(); ++i)
        {
            const int block = data.group->pro_block_index(rank)[i];
            file << "| " << block << " | " << i << " | "
                 << data.ptr->blk(block).ijkmax[0] << "x" << data.ptr->blk(block).ijkmax[1] << "x"
                 << data.ptr->blk(block).ijkmax[2] << " | [block](../" << BlockPath(data, block, "") << ") |\n";
        }

        file << "\n## Intra-Rank Interfaces\n\n";
        file << "| Patch | Global block | Local block | Face | Kind | Target / Boundary |\n|---|---:|---:|---|---|---|\n";
        for (std::size_t i = 0; i < data.rank_patch_ids[rank].size(); ++i)
        {
            const Patch &patch = data.patches[data.rank_patch_ids[rank][i]];
            if (patch.target_block < 0 || patch.target_myid != rank)
                continue;
            file << "| [patch " << patch.id << "](../" << PatchPath(data, patch.id, "") << ") | "
                 << patch.source_block << " | " << patch.source_local << " | " << FaceLabel(patch.source_face) << " | "
                 << patch.kind << " | ";
            file << "myid " << patch.target_myid << ", block " << patch.target_block << ", " << FaceLabel(patch.target_face);
            if (patch.interface_id >= 0)
                file << ", [interface " << patch.interface_id << "](../" << InterfacePath(data, patch.interface_id, "") << ")";
            file << " |\n";
        }

        file << "\n## Inter-Rank Interfaces\n\n";
        file << "| Patch | Global block | Local block | Face | Remote myid | Target block | Interface |\n|---|---:|---:|---|---:|---:|---|\n";
        for (std::size_t i = 0; i < data.rank_patch_ids[rank].size(); ++i)
        {
            const Patch &patch = data.patches[data.rank_patch_ids[rank][i]];
            if (patch.target_block < 0 || patch.target_myid == rank)
                continue;
            file << "| [patch " << patch.id << "](../" << PatchPath(data, patch.id, "") << ") | "
                 << patch.source_block << " | " << patch.source_local << " | " << FaceLabel(patch.source_face) << " | "
                 << patch.target_myid << " | " << patch.target_block << " | ";
            if (patch.interface_id >= 0)
                file << "[interface " << patch.interface_id << "](../" << InterfacePath(data, patch.interface_id, "") << ")";
            file << " |\n";
        }

        file << "\n## Physical Boundaries\n\n";
        file << "| Patch | Global block | Local block | Face | Boundary |\n|---|---:|---:|---|---|\n";
        for (std::size_t i = 0; i < data.rank_patch_ids[rank].size(); ++i)
        {
            const Patch &patch = data.patches[data.rank_patch_ids[rank][i]];
            if (patch.target_block >= 0)
                continue;
            file << "| [patch " << patch.id << "](../" << PatchPath(data, patch.id, "") << ") | "
                 << patch.source_block << " | " << patch.source_local << " | " << FaceLabel(patch.source_face) << " | "
                 << MdEscape(patch.boundary_name) << " |\n";
        }
    }
}

void WriteBlockPages(ReportData &data)
{
    for (int block = 0; block < data.nblock; ++block)
    {
        std::ofstream file(("OUTPUT/topology_report/" + BlockPath(data, block, "")).c_str());
        file << "# Block " << block << "\n\n";
        file << "[Index](../index.md) | [myid " << data.group->block_proc_number(block) << "](../"
             << RankPath(data, data.group->block_proc_number(block), "") << ")\n\n";
        file << "- myid: " << data.group->block_proc_number(block) << "\n";
        file << "- local block id: " << data.group->block_proc_index(block) << "\n";
        file << "- size: " << data.ptr->blk(block).ijkmax[0] << " x " << data.ptr->blk(block).ijkmax[1]
             << " x " << data.ptr->blk(block).ijkmax[2] << "\n\n";

        file << "## Hierarchy\n\n";
        file << "| Part | local sub/sup | original block | original sub/sup | Transform | name |\n|---:|---|---:|---|---|---|\n";
        for (std::size_t i = 0; i < data.ptr->blk(block).hiera.size(); ++i)
        {
            const coordinate_blk &h = data.ptr->blk(block).hiera[i];
            file << "| " << i << " | [" << h.sub[0] << ":" << h.sup[0] << ", "
                 << h.sub[1] << ":" << h.sup[1] << ", " << h.sub[2] << ":" << h.sup[2] << "] | "
                 << h.oring_num << " | [" << h.oring_sub[0] << ":" << h.oring_sup[0] << ", "
                 << h.oring_sub[1] << ":" << h.oring_sup[1] << ", " << h.oring_sub[2] << ":" << h.oring_sup[2]
                 << "] | [" << h.Transform[0] << "," << h.Transform[1] << "," << h.Transform[2] << "] | "
                 << MdEscape(h.blk_name) << " |\n";
        }

        file << "\n## Faces\n\n";
        file << "| Face | Patch count | Connections / boundaries | Link |\n|---|---:|---|---|\n";
        for (int face = 0; face < 6; ++face)
        {
            const std::vector<int> patch_ids = FacePatchIds(data, block, face);
            std::set<std::string> targets;
            for (std::size_t i = 0; i < patch_ids.size(); ++i)
                targets.insert(OwnerLabel(data, data.patches[patch_ids[i]]));
            file << "| " << FaceLabel(face) << " | " << patch_ids.size() << " | ";
            for (std::set<std::string>::iterator it = targets.begin(); it != targets.end(); ++it)
                file << (it == targets.begin() ? "" : ", ") << MdEscape(*it);
            file << " | [face](../" << FacePath(data, block, face, "") << ") |\n";
        }
    }
}

void WritePatchPages(ReportData &data)
{
    for (std::size_t i = 0; i < data.patches.size(); ++i)
    {
        const Patch &patch = data.patches[i];
        std::ofstream file(("OUTPUT/topology_report/" + PatchPath(data, patch.id, "")).c_str());
        file << "# Patch " << patch.id << "\n\n";
        file << "[Index](../index.md) | [Block](../" << BlockPath(data, patch.source_block, "") << ")";
        if (patch.source_face != FACE_UNKNOWN)
            file << " | [Face](../" << FacePath(data, patch.source_block, patch.source_face, "") << ")";
        file << "\n\n";
        file << "- kind: " << patch.kind << "\n";
        file << "- source: block " << patch.source_block << ", myid " << patch.source_myid
             << ", local " << patch.source_local << ", face " << FaceLabel(patch.source_face)
             << ", box `" << BoxText(patch.source_box) << "`\n";
        if (patch.target_block >= 0)
        {
            file << "- target: block " << patch.target_block << ", myid " << patch.target_myid
                 << ", local " << patch.target_local << ", face " << FaceLabel(patch.target_face)
                 << ", box `" << BoxText(patch.target_box) << "`\n";
            file << "- tangential mapping: source `" << FaceAxes(patch.source_face) << "` -> target `"
                 << FaceAxes(patch.target_face) << "`\n";
        }
        else
        {
            file << "- boundary: " << MdEscape(patch.boundary_name) << " (" << patch.boundary_kind << ")\n";
        }
        file << "- patch size: " << PatchSize(patch) << "\n";
        if (patch.reverse_patch >= 0)
            file << "- reverse patch: [patch " << patch.reverse_patch << "](../" << PatchPath(data, patch.reverse_patch, "") << ")\n";
        if (patch.interface_id >= 0)
            file << "- interface: [interface " << patch.interface_id << "](../" << InterfacePath(data, patch.interface_id, "") << ")\n";
    }
}

void WriteInterfacePages(ReportData &data)
{
    for (std::size_t i = 0; i < data.interfaces.size(); ++i)
    {
        const InterfaceInfo &iface = data.interfaces[i];
        std::ofstream file(("OUTPUT/topology_report/" + InterfacePath(data, iface.id, "")).c_str());
        file << "# Interface " << iface.id << "\n\n";
        file << "[Index](../index.md)\n\n";
        file << "- block A: [B" << iface.block_a << "](../" << BlockPath(data, iface.block_a, "") << "), myid " << iface.myid_a << "\n";
        file << "- block B: [B" << iface.block_b << "](../" << BlockPath(data, iface.block_b, "") << "), myid " << iface.myid_b << "\n";
        file << "- inter-rank: " << (iface.inter_rank ? "yes" : "no") << "\n";
        file << "- periodic: " << (iface.periodic ? "yes" : "no") << "\n";
        file << "- patch A: [patch " << iface.patch_a << "](../" << PatchPath(data, iface.patch_a, "") << ")\n";
        if (iface.patch_b >= 0)
            file << "- patch B: [patch " << iface.patch_b << "](../" << PatchPath(data, iface.patch_b, "") << ")\n";
        else
            file << "- patch B: not found\n";
    }
}

void WriteJson(ReportData &data)
{
    std::ofstream file("OUTPUT/topology_report/topology.json");
    file << "{\n";
    file << "  \"blocks\": [\n";
    for (int block = 0; block < data.nblock; ++block)
    {
        file << "    {\"id\": " << block << ", \"myid\": " << data.group->block_proc_number(block)
             << ", \"local_block\": " << data.group->block_proc_index(block)
             << ", \"size\": [" << data.ptr->blk(block).ijkmax[0] << ", " << data.ptr->blk(block).ijkmax[1]
             << ", " << data.ptr->blk(block).ijkmax[2] << "]}";
        file << (block + 1 == data.nblock ? "\n" : ",\n");
    }
    file << "  ],\n  \"patches\": [\n";
    for (std::size_t i = 0; i < data.patches.size(); ++i)
    {
        const Patch &p = data.patches[i];
        file << "    {\"id\": " << p.id << ", \"kind\": \"" << JsonEscape(p.kind) << "\", \"source_block\": "
             << p.source_block << ", \"source_face\": \"" << FaceLabel(p.source_face) << "\", \"target_block\": "
             << p.target_block << ", \"target_face\": \"" << FaceLabel(p.target_face) << "\", \"boundary\": \""
             << JsonEscape(p.boundary_name) << "\", \"reverse_patch\": " << p.reverse_patch
             << ", \"interface\": " << p.interface_id << "}";
        file << (i + 1 == data.patches.size() ? "\n" : ",\n");
    }
    file << "  ],\n  \"interfaces\": [\n";
    for (std::size_t i = 0; i < data.interfaces.size(); ++i)
    {
        const InterfaceInfo &iface = data.interfaces[i];
        file << "    {\"id\": " << iface.id << ", \"patch_a\": " << iface.patch_a << ", \"patch_b\": "
             << iface.patch_b << ", \"block_a\": " << iface.block_a << ", \"block_b\": " << iface.block_b
             << ", \"myid_a\": " << iface.myid_a << ", \"myid_b\": " << iface.myid_b
             << ", \"inter_rank\": " << (iface.inter_rank ? "true" : "false")
             << ", \"periodic\": " << (iface.periodic ? "true" : "false") << "}";
        file << (i + 1 == data.interfaces.size() ? "\n" : ",\n");
    }
    file << "  ],\n  \"warnings\": [\n";
    for (std::size_t i = 0; i < data.warnings.size(); ++i)
    {
        file << "    {\"kind\": \"" << JsonEscape(data.warnings[i].kind) << "\", \"message\": \""
             << JsonEscape(data.warnings[i].message) << "\", \"link\": \"" << JsonEscape(data.warnings[i].link) << "\"}";
        file << (i + 1 == data.warnings.size() ? "\n" : ",\n");
    }
    file << "  ]\n}\n";
}
} // namespace

void TopologyReport::Generate(Preprocess_Data_Structure *ptr, Preprocess_Group *group, Param *param)
{
    ReportData data;
    data.ptr = ptr;
    data.group = group;
    data.param = param;
    data.nblock = ptr->blk.Getsize1();
    data.nproc = group->proc_num.Getsize1();
    data.block_width = WidthForMax(std::max(0, data.nblock - 1));
    data.rank_width = WidthForMax(std::max(0, data.nproc - 1));
    data.patch_width = 1;
    data.interface_width = 1;
    data.ascii_layout = !param->HasBoo("topology_report_ascii_layout") || param->GetBoo("topology_report_ascii_layout");

    EnsureReportDirectories();
    CollectPatches(data);
    data.patch_width = WidthForMax(std::max(0, static_cast<int>(data.patches.size()) - 1));
    DiagnosePatchRanges(data);
    LinkReversePatchesAndInterfaces(data);
    data.interface_width = WidthForMax(std::max(0, static_cast<int>(data.interfaces.size()) - 1));
    PrepareFaceReports(data);

    for (std::size_t i = 0; i < data.faces.size(); ++i)
        AnalyzeAndWriteFace(data, data.faces[i]);

    WriteCsv(data);
    WriteRankPages(data);
    WriteBlockPages(data);
    WritePatchPages(data);
    WriteInterfacePages(data);
    WriteIndexAndSummary(data);
    WriteWarnings(data);
    WriteJson(data);

    std::cout << "\t-->Topology report has been generated in ./OUTPUT/topology_report\n";
}
