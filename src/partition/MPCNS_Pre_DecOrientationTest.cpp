#include "partition/MPCNS_Pre_DecOrientationTest.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
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
struct Transform
{
    int axis[3];
};

struct Box
{
    int sub[3];
    int sup[3];
};

struct DirectedPatch
{
    int source_block;
    int target_block;
    Box source_box;
    Box target_box;
};

struct CoverageResult
{
    std::set<std::string> keys;
    std::map<std::string, int> counts;
};

bool IsIdentity(const Transform &t)
{
    return t.axis[0] == 1 && t.axis[1] == 2 && t.axis[2] == 3;
}

std::string TransformText(const Transform &t, int dimension)
{
    std::ostringstream s;
    for (int i = 0; i < dimension; ++i)
    {
        if (i)
            s << " ";
        s << (t.axis[i] > 0 ? "+" : "") << t.axis[i];
    }
    if (dimension == 2)
        s << " +3";
    return s.str();
}

int PermutationSign(const std::vector<int> &perm)
{
    int inversions = 0;
    for (std::size_t i = 0; i < perm.size(); ++i)
        for (std::size_t j = i + 1; j < perm.size(); ++j)
            if (perm[i] > perm[j])
                inversions++;
    return (inversions % 2 == 0) ? 1 : -1;
}

std::vector<Transform> GenerateTransforms(int dimension)
{
    std::vector<Transform> transforms;
    if (dimension == 2)
    {
        const int values[4][3] = {
            {1, 2, 3},
            {2, -1, 3},
            {-1, -2, 3},
            {-2, 1, 3}};
        for (int i = 0; i < 4; ++i)
        {
            Transform t;
            for (int d = 0; d < 3; ++d)
                t.axis[d] = values[i][d];
            transforms.push_back(t);
        }
        return transforms;
    }

    std::vector<int> perm;
    perm.push_back(1);
    perm.push_back(2);
    perm.push_back(3);
    do
    {
        const int perm_sign = PermutationSign(perm);
        for (int sx = -1; sx <= 1; sx += 2)
            for (int sy = -1; sy <= 1; sy += 2)
                for (int sz = -1; sz <= 1; sz += 2)
                {
                    if (perm_sign * sx * sy * sz != 1)
                        continue;
                    Transform t;
                    t.axis[0] = sx * perm[0];
                    t.axis[1] = sy * perm[1];
                    t.axis[2] = sz * perm[2];
                    transforms.push_back(t);
                }
    } while (std::next_permutation(perm.begin(), perm.end()));
    return transforms;
}

std::string FaceLabel(int face)
{
    switch (face)
    {
    case 0:
        return "xi-";
    case 1:
        return "xi+";
    case 2:
        return "eta-";
    case 3:
        return "eta+";
    case 4:
        return "zeta-";
    case 5:
        return "zeta+";
    default:
        return "unknown";
    }
}

std::vector<std::string> TheoreticalKeys(int dimension)
{
    const int face_count = (dimension == 2) ? 4 : 6;
    std::vector<std::string> keys;
    for (int i = 0; i < face_count; ++i)
        for (int j = 0; j < face_count; ++j)
            keys.push_back(FaceLabel(i) + " -> " + FaceLabel(j));
    return keys;
}

int MapSignedIndex(int old_q, int old_size, int signed_axis)
{
    const int old_abs = std::abs(old_q);
    const int old_sign = old_q < 0 ? -1 : 1;
    if (signed_axis > 0)
        return old_sign * old_abs;
    return old_sign * (old_size + 1 - old_abs);
}

void MapBoxPreserveDirection(const int old_sub[3], const int old_sup[3], const int old_size[3],
                             const Transform &t, int new_sub[3], int new_sup[3])
{
    for (int new_axis = 0; new_axis < 3; ++new_axis)
    {
        const int old_axis = std::abs(t.axis[new_axis]) - 1;
        new_sub[new_axis] = MapSignedIndex(old_sub[old_axis], old_size[old_axis], t.axis[new_axis]);
        new_sup[new_axis] = MapSignedIndex(old_sup[old_axis], old_size[old_axis], t.axis[new_axis]);
    }
}

void MapBoxSortedPositive(const int old_sub[3], const int old_sup[3], const int old_size[3],
                          const Transform &t, int new_sub[3], int new_sup[3])
{
    int mapped_sub[3], mapped_sup[3];
    MapBoxPreserveDirection(old_sub, old_sup, old_size, t, mapped_sub, mapped_sup);
    for (int d = 0; d < 3; ++d)
    {
        new_sub[d] = std::min(std::abs(mapped_sub[d]), std::abs(mapped_sup[d]));
        new_sup[d] = std::max(std::abs(mapped_sub[d]), std::abs(mapped_sup[d]));
    }
}

int DetectFace(const int sub[3], const int sup[3], const int size[3], int dimension)
{
    for (int d = 0; d < dimension; ++d)
    {
        const int lo = std::min(std::abs(sub[d]), std::abs(sup[d]));
        const int hi = std::max(std::abs(sub[d]), std::abs(sup[d]));
        if (lo == hi && lo == 1)
            return 2 * d;
        if (lo == hi && lo == size[d])
            return 2 * d + 1;
    }
    return -1;
}

int TransformedFace(const Box &box, const int old_size[3], const Transform &t, int dimension)
{
    int new_sub[3], new_sup[3], new_size[3];
    MapBoxPreserveDirection(box.sub, box.sup, old_size, t, new_sub, new_sup);
    for (int d = 0; d < 3; ++d)
        new_size[d] = old_size[std::abs(t.axis[d]) - 1];
    return DetectFace(new_sub, new_sup, new_size, dimension);
}

std::vector<DirectedPatch> CollectDirectedPatches(Preprocess_Data_Structure *data)
{
    std::vector<DirectedPatch> patches;
    const int nblock = data->blk.Getsize1();
    for (int i = 0; i < nblock; ++i)
    {
        for (int j = 0; j < nblock; ++j)
        {
            Inner &inner = data->inp(i, j);
            for (std::size_t k = 0; k < inner.inner.size(); ++k)
            {
                DirectedPatch patch;
                patch.source_block = i;
                patch.target_block = j;
                for (int d = 0; d < 3; ++d)
                {
                    patch.source_box.sub[d] = inner.inner[k].sub[d];
                    patch.source_box.sup[d] = inner.inner[k].sup[d];
                    patch.target_box.sub[d] = inner.inner[k].tar_sub[d];
                    patch.target_box.sup[d] = inner.inner[k].tar_sup[d];
                }
                patches.push_back(patch);
            }
        }
    }
    return patches;
}

CoverageResult EvaluateCoverage(const std::vector<DirectedPatch> &patches,
                                const std::vector<Transform> &block_transforms,
                                const std::vector<std::vector<int> > &old_sizes,
                                int dimension)
{
    CoverageResult result;
    for (std::size_t i = 0; i < patches.size(); ++i)
    {
        const DirectedPatch &patch = patches[i];
        const int source_face = TransformedFace(patch.source_box, &old_sizes[patch.source_block][0],
                                                block_transforms[patch.source_block], dimension);
        const int target_face = TransformedFace(patch.target_box, &old_sizes[patch.target_block][0],
                                                block_transforms[patch.target_block], dimension);
        if (source_face < 0 || target_face < 0)
            continue;
        const std::string key = FaceLabel(source_face) + " -> " + FaceLabel(target_face);
        result.keys.insert(key);
        result.counts[key]++;
    }
    return result;
}

std::vector<int> MissingKeys(const std::set<std::string> &covered, int dimension)
{
    std::vector<int> missing_indices;
    const std::vector<std::string> all_keys = TheoreticalKeys(dimension);
    for (std::size_t i = 0; i < all_keys.size(); ++i)
        if (covered.count(all_keys[i]) == 0)
            missing_indices.push_back(static_cast<int>(i));
    return missing_indices;
}

std::vector<std::string> MissingKeyNames(const std::set<std::string> &covered, int dimension)
{
    std::vector<std::string> missing;
    const std::vector<std::string> all_keys = TheoreticalKeys(dimension);
    for (std::size_t i = 0; i < all_keys.size(); ++i)
        if (covered.count(all_keys[i]) == 0)
            missing.push_back(all_keys[i]);
    return missing;
}

std::string BoxText(const int sub[3], const int sup[3])
{
    std::ostringstream s;
    s << "[" << sub[0] << ":" << sup[0] << ", "
      << sub[1] << ":" << sup[1] << ", "
      << sub[2] << ":" << sup[2] << "]";
    return s.str();
}

int NegativeAxisCount(const int sub[3], const int sup[3])
{
    int count = 0;
    for (int d = 0; d < 3; ++d)
        if (sub[d] < 0 || sup[d] < 0)
            count++;
    return count;
}

void ValidateInnerBoxFormat(Preprocess_Data_Structure *data, int dimension)
{
    if (dimension == 2)
        return;

    const int nblock = data->blk.Getsize1();
    for (int iblock = 0; iblock < nblock; ++iblock)
        for (int jblock = 0; jblock < nblock; ++jblock)
        {
            Inner &inner = data->inp(iblock, jblock);
            for (std::size_t iface = 0; iface < inner.inner.size(); ++iface)
            {
                coordinate_inp &patch = inner.inner[iface];
                const int source_negative_axes = NegativeAxisCount(patch.sub, patch.sup);
                const int target_negative_axes = NegativeAxisCount(patch.tar_sub, patch.tar_sup);
                if (source_negative_axes != 1 || target_negative_axes != 1)
                {
                    std::cout << "#Fatal Error: DEC orientation test generated invalid 3D inner box sign format.\n";
                    std::cout << "\tEach inner output row must have exactly one negative ijk direction.\n";
                    std::cout << "\tsource_block=" << iblock << " target_block=" << jblock
                              << " patch_index=" << iface << "\n";
                    std::cout << "\tsource_box=" << BoxText(patch.sub, patch.sup)
                              << " negative_axis_count=" << source_negative_axes << "\n";
                    std::cout << "\ttarget_box=" << BoxText(patch.tar_sub, patch.tar_sup)
                              << " negative_axis_count=" << target_negative_axes << "\n";
                    exit(-1);
                }
            }
        }
}

void EnsureDirectory(const std::string &path)
{
    if (mkdir(path.c_str(), 0755) != 0 && errno != EEXIST)
    {
        std::cout << "#Fatal Error: Cannot create directory:\t" << path << std::endl;
        exit(-1);
    }
}

void WriteManifest(const std::vector<Transform> &selected,
                   const CoverageResult &coverage,
                   int theoretical_count,
                   int directed_patch_count,
                   int target_count,
                   int dimension,
                   bool success)
{
    EnsureDirectory("OUTPUT");
    std::ofstream file("OUTPUT/dec_orientation_manifest.txt");
    file << "if_dec_orientation_test=1\n";
    file << "final_block_count=" << selected.size() << "\n";
    file << "directed_interface_patch_count=" << directed_patch_count << "\n";
    file << "theoretical_key_count=" << theoretical_count << "\n";
    file << "target_count=" << target_count << "\n";
    file << "achieved_count=" << coverage.keys.size() << "\n";
    file << "status=" << (success ? "success" : "failed") << "\n\n";
    file << "[block_transforms]\n";
    for (std::size_t i = 0; i < selected.size(); ++i)
        file << "block " << i << " : " << TransformText(selected[i], dimension) << "\n";
    file << "\n[face_pair_coverage]\n";
    const std::vector<std::string> all_keys = TheoreticalKeys(dimension);
    for (std::size_t i = 0; i < all_keys.size(); ++i)
    {
        std::map<std::string, int>::const_iterator iter = coverage.counts.find(all_keys[i]);
        file << all_keys[i] << " " << (iter == coverage.counts.end() ? 0 : iter->second) << "\n";
    }
    file << "\n[missing_keys]\n";
    const std::vector<std::string> missing = MissingKeyNames(coverage.keys, dimension);
    for (std::size_t i = 0; i < missing.size(); ++i)
        file << missing[i] << "\n";

    EnsureDirectory("OUTPUT/topology_report");
    std::ofstream md("OUTPUT/topology_report/orientation_coverage.md");
    md << "# DEC Orientation Coverage\n\n";
    md << "- directed interface patches: " << directed_patch_count << "\n";
    md << "- theoretical key count: " << theoretical_count << "\n";
    md << "- target count: " << target_count << "\n";
    md << "- achieved count: " << coverage.keys.size() << "\n";
    md << "- status: " << (success ? "success" : "failed") << "\n\n";
    md << "## Block Transforms\n\n";
    md << "| Block | Transform |\n|---:|---|\n";
    for (std::size_t i = 0; i < selected.size(); ++i)
        md << "| " << i << " | `" << TransformText(selected[i], dimension) << "` |\n";
    md << "\n## Face-Pair Counts\n\n";
    md << "| Face pair | Count |\n|---|---:|\n";
    for (std::size_t i = 0; i < all_keys.size(); ++i)
    {
        std::map<std::string, int>::const_iterator iter = coverage.counts.find(all_keys[i]);
        md << "| `" << all_keys[i] << "` | " << (iter == coverage.counts.end() ? 0 : iter->second) << " |\n";
    }
    md << "\n## Missing Keys\n\n";
    for (std::size_t i = 0; i < missing.size(); ++i)
        md << "- `" << missing[i] << "`\n";
}

void ApplyTransforms(Preprocess_Data_Structure *data,
                     const std::vector<Transform> &selected,
                     const std::vector<std::vector<int> > &old_sizes,
                     int dimension)
{
    const int nblock = data->blk.Getsize1();
    for (int iblock = 0; iblock < nblock; ++iblock)
    {
        block_info &block = data->blk(iblock);
        const Transform &t = selected[iblock];
        for (std::size_t i = 0; i < block.bound.size(); ++i)
        {
            int old_sub[3], old_sup[3];
            for (int d = 0; d < 3; ++d)
            {
                old_sub[d] = block.bound[i].sub[d];
                old_sup[d] = block.bound[i].sup[d];
            }
            MapBoxPreserveDirection(old_sub, old_sup, &old_sizes[iblock][0], t,
                                    block.bound[i].sub, block.bound[i].sup);
        }

        for (std::size_t i = 0; i < block.hiera.size(); ++i)
        {
            coordinate_blk &h = block.hiera[i];
            int old_sub[3], old_sup[3], old_transform[3];
            for (int d = 0; d < 3; ++d)
            {
                old_sub[d] = h.sub[d];
                old_sup[d] = h.sup[d];
                old_transform[d] = h.Transform[d];
            }
            MapBoxSortedPositive(old_sub, old_sup, &old_sizes[iblock][0], t, h.sub, h.sup);
            for (int new_axis = 0; new_axis < 3; ++new_axis)
            {
                const int old_axis = std::abs(t.axis[new_axis]) - 1;
                h.Transform[new_axis] = (t.axis[new_axis] > 0 ? 1 : -1) * old_transform[old_axis];
            }
        }
    }

    for (int iblock = 0; iblock < nblock; ++iblock)
    {
        for (int d = 0; d < 3; ++d)
            data->blk(iblock).ijkmax[d] = old_sizes[iblock][std::abs(selected[iblock].axis[d]) - 1];
        if (dimension == 2)
            data->blk(iblock).ijkmax[2] = 1;
    }

    for (int iblock = 0; iblock < nblock; ++iblock)
        for (int jblock = 0; jblock < nblock; ++jblock)
        {
            Inner &inner = data->inp(iblock, jblock);
            for (std::size_t iface = 0; iface < inner.inner.size(); ++iface)
            {
                coordinate_inp &patch = inner.inner[iface];
                int old_sub[3], old_sup[3], old_tar_sub[3], old_tar_sup[3];
                for (int d = 0; d < 3; ++d)
                {
                    old_sub[d] = patch.sub[d];
                    old_sup[d] = patch.sup[d];
                    old_tar_sub[d] = patch.tar_sub[d];
                    old_tar_sup[d] = patch.tar_sup[d];
                }
                MapBoxPreserveDirection(old_sub, old_sup, &old_sizes[iblock][0], selected[iblock],
                                        patch.sub, patch.sup);
                MapBoxPreserveDirection(old_tar_sub, old_tar_sup, &old_sizes[jblock][0], selected[jblock],
                                        patch.tar_sub, patch.tar_sup);
                patch.volume = 1;
                for (int d = 0; d < 3; ++d)
                    patch.volume *= std::abs(std::abs(patch.sub[d]) - std::abs(patch.sup[d])) + 1;
            }
        }
}
} // namespace

void DecOrientationTest::ApplyIfEnabled(Preprocess_Data_Structure *data, Param *param)
{
    if (!param->HasBoo("if_dec_orientation_test") || !param->GetBoo("if_dec_orientation_test"))
        return;

    const int dimension = param->GetInt("dimension");
    const int nblock = data->blk.Getsize1();
    const std::vector<DirectedPatch> patches = CollectDirectedPatches(data);
    if (patches.empty())
    {
        std::cout << "#Fatal Error: if_dec_orientation_test=1 requires at least one directed inner interface patch after split.\n";
        exit(-1);
    }

    std::vector<std::vector<int> > old_sizes(nblock, std::vector<int>(3, 1));
    for (int i = 0; i < nblock; ++i)
        for (int d = 0; d < 3; ++d)
            old_sizes[i][d] = data->blk(i).ijkmax[d];

    const std::vector<Transform> candidates = GenerateTransforms(dimension);
    Transform identity;
    identity.axis[0] = 1;
    identity.axis[1] = 2;
    identity.axis[2] = 3;
    std::vector<Transform> selected(nblock, identity);

    const int theoretical_count = (dimension == 2) ? 16 : 36;
    const int target_count = std::min(theoretical_count, static_cast<int>(patches.size()));

    std::vector<int> incident(nblock, 0);
    for (std::size_t i = 0; i < patches.size(); ++i)
    {
        incident[patches[i].source_block]++;
        incident[patches[i].target_block]++;
    }
    std::vector<int> order;
    for (int i = 0; i < nblock; ++i)
        order.push_back(i);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        if (incident[a] != incident[b])
            return incident[a] > incident[b];
        return a < b;
    });

    CoverageResult current = EvaluateCoverage(patches, selected, old_sizes, dimension);
    for (std::size_t oi = 0; oi < order.size() && static_cast<int>(current.keys.size()) < target_count; ++oi)
    {
        const int block = order[oi];
        Transform best = selected[block];
        CoverageResult best_coverage = current;
        int best_new_keys = 0;

        for (std::size_t ci = 0; ci < candidates.size(); ++ci)
        {
            std::vector<Transform> trial = selected;
            trial[block] = candidates[ci];
            CoverageResult trial_coverage = EvaluateCoverage(patches, trial, old_sizes, dimension);
            int new_keys = 0;
            for (std::set<std::string>::const_iterator it = trial_coverage.keys.begin(); it != trial_coverage.keys.end(); ++it)
                if (current.keys.count(*it) == 0)
                    new_keys++;

            bool take = false;
            if (trial_coverage.keys.size() > best_coverage.keys.size())
                take = true;
            else if (trial_coverage.keys.size() == best_coverage.keys.size())
            {
                if (new_keys > best_new_keys)
                    take = true;
                else if (new_keys == best_new_keys)
                {
                    const bool candidate_non_identity = !IsIdentity(candidates[ci]);
                    const bool best_non_identity = !IsIdentity(best);
                    if (candidate_non_identity != best_non_identity)
                        take = candidate_non_identity;
                    else if (TransformText(candidates[ci], dimension) < TransformText(best, dimension))
                        take = true;
                }
            }

            if (take)
            {
                best = candidates[ci];
                best_coverage = trial_coverage;
                best_new_keys = new_keys;
            }
        }

        selected[block] = best;
        current = best_coverage;
    }

    const bool success = static_cast<int>(current.keys.size()) >= target_count;
    WriteManifest(selected, current, theoretical_count, static_cast<int>(patches.size()), target_count, dimension, success);
    if (!success)
    {
        std::cout << "#Fatal Error: DEC orientation test could not reach requested face-pair coverage.\n";
        std::cout << "\ttheoretical_key_count=" << theoretical_count << "\n";
        std::cout << "\tdirected_interface_patch_count=" << patches.size() << "\n";
        std::cout << "\ttarget_count=" << target_count << "\n";
        std::cout << "\tachieved_count=" << current.keys.size() << "\n";
        std::cout << "\tfinal_block_count=" << nblock << "\n";
        std::cout << "\tmissing face-pair keys:\n";
        const std::vector<std::string> missing = MissingKeyNames(current.keys, dimension);
        for (std::size_t i = 0; i < missing.size(); ++i)
            std::cout << "\t\t" << missing[i] << "\n";
        std::cout << "\tSuggestion: increase split density or use a grid with more block interfaces.\n";
        exit(-1);
    }

    ApplyTransforms(data, selected, old_sizes, dimension);
    ValidateInnerBoxFormat(data, dimension);
    std::cout << "\t-->DEC orientation test pass achieved " << current.keys.size()
              << "/" << target_count << " directed face-pair keys.\n";
}
