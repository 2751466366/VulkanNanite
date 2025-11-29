#pragma once
#include "meshoptimizer.h"
#include "metis.h"
#include "MyDataType.h"

/**
xadj 数组
    含义：xadj 是一个数组，用于存储每个顶点的邻接表的起始位置。
    格式：xadj[i] 表示顶点 i 的邻接表在 adjncy 数组中的起始索引。
    大小：xadj 的大小为 nvtxs + 1，其中 nvtxs 是图的顶点数。
        示例：假设图有 3 个顶点，邻接表如下：
        顶点 0 的邻接顶点：1, 2
        顶点 1 的邻接顶点：0, 2
        顶点 2 的邻接顶点：0, 1
        则 xadj 数组为：[0, 2, 4, 6]，表示：
        顶点 0 的邻接表从索引 0 开始，到索引 2 结束（不包括索引 2）。
        顶点 1 的邻接表从索引 2 开始，到索引 4 结束（不包括索引 4）。
        顶点 2 的邻接表从索引 4 开始，到索引 6 结束（不包括索引 6）。
adjncy 数组
    含义：adjncy 是一个数组，用于存储所有顶点的邻接顶点。
    格式：adjncy 数组中的元素是顶点的编号，表示与当前顶点相邻的顶点。
    大小：adjncy 的大小为图中所有顶点的度数之和。
        示例：对于上述图，adjncy 数组为：[1, 2, 0, 2, 0, 1]，表示：
        顶点 0 的邻接顶点：adjncy[0] = 1 和 adjncy[1] = 2。
        顶点 1 的邻接顶点：adjncy[2] = 0 和 adjncy[3] = 2。
        顶点 2 的邻接顶点：adjncy[4] = 0 和 adjncy[5] = 1。

示例
假设有一个图，顶点数为 4，边如下：
顶点 0 与顶点 1 和 2 相连
顶点 1 与顶点 0 和 3 相连
顶点 2 与顶点 0 和 3 相连
顶点 3 与顶点 1 和 2 相连
则：
xadj 数组为：[0, 2, 4, 6, 8]
adjncy 数组为：[1, 2, 0, 3, 0, 3, 1, 2]
*/

class NaniteMgr
{
public:
    const static size_t kClusterSize = 128;
    const static size_t kGroupSize = 8;
    const static int kMetisSlop = 2;

    static void RemapMesh(std::vector<Vertex>& vertices,
                          std::vector<uint32_t>& indices) {
        size_t vertex_size = sizeof(Vertex);
        size_t index_count = indices.size();
        size_t vertex_count = vertices.size();

        std::vector<unsigned int> remap(index_count);
        //remap.resize(vertex_count);
        size_t uniqueVertexCount = meshopt_generateVertexRemap(
            remap.data(),
            indices.data(),
            index_count,
            vertices.data(),
            vertex_count,
            sizeof(Vertex)
        );

        std::vector<Vertex> optimizedVertices(uniqueVertexCount);
        std::vector<uint32_t> optimizedIndices(index_count);

        meshopt_remapVertexBuffer(
            optimizedVertices.data(),
            vertices.data(),
            vertex_count,
            sizeof(Vertex),
            remap.data()
        );

        meshopt_remapIndexBuffer(
            optimizedIndices.data(),
            indices.data(),
            index_count,
            remap.data()
        );

        vertices = std::move(optimizedVertices);
        indices = std::move(optimizedIndices);
    }


    static std::vector<Cluster> ClusterizeMetis(
        const std::vector<uint32_t>& remap,
        const std::vector<unsigned int>& indices) {

        std::vector<std::vector<int> > trilist(remap.size());
        for (size_t i = 0; i < indices.size(); ++i)
            trilist[remap[i]].push_back(int(i / 3));

        std::vector<int> xadj(indices.size() / 3 + 1);
        std::vector<int> adjncy;
        std::vector<int> adjwgt;
        std::vector<int> part(indices.size() / 3);

        std::vector<int> scratch;

        // 找相邻的三角形
        for (size_t i = 0; i < indices.size() / 3; ++i)
        {
            unsigned int a = remap[i * 3 + 0], b = remap[i * 3 + 1], c = remap[i * 3 + 2];

            scratch.clear();
            scratch.insert(scratch.end(), trilist[a].begin(), trilist[a].end());
            scratch.insert(scratch.end(), trilist[b].begin(), trilist[b].end());
            scratch.insert(scratch.end(), trilist[c].begin(), trilist[c].end());
            std::sort(scratch.begin(), scratch.end());

            for (size_t j = 0; j < scratch.size(); ++j)
            {
                if (scratch[j] == int(i))
                    continue;

                if (j == 0 || scratch[j] != scratch[j - 1])
                {
                    adjncy.push_back(scratch[j]);
                    adjwgt.push_back(1);
                }
                else if (j != 0)
                {
                    assert(scratch[j] == scratch[j - 1]);
                    adjwgt.back()++;
                }
            }

            xadj[i + 1] = int(adjncy.size());
        }

        int options[METIS_NOPTIONS];
        METIS_SetDefaultOptions(options);
        options[METIS_OPTION_SEED] = 42;
        options[METIS_OPTION_UFACTOR] = 1; // minimize partition imbalance

        // since Metis can't enforce partition sizes, add a little slop to reduce the change we need to split results further
        int nvtxs = int(indices.size() / 3);
        int ncon = 1;
        int nparts = int(indices.size() / 3 + (kClusterSize - kMetisSlop) - 1) / (kClusterSize - kMetisSlop);
        int edgecut = 0;

        // not sure why this is a special case that we need to handle but okay metis
        if (nparts > 1)
        {
            int r = METIS_PartGraphRecursive(&nvtxs, &ncon, &xadj[0], &adjncy[0], NULL, NULL, &adjwgt[0], &nparts, NULL, NULL, options, &edgecut, &part[0]);
            assert(r == METIS_OK);
            (void)r;
        }

        std::vector<Cluster> result(nparts);

        for (size_t i = 0; i < indices.size() / 3; ++i)
        {
            result[part[i]].indices.push_back(indices[i * 3 + 0]);
            result[part[i]].indices.push_back(indices[i * 3 + 1]);
            result[part[i]].indices.push_back(indices[i * 3 + 2]);
        }

        for (int i = 0; i < nparts; ++i)
        {
            result[i].parent.error = FLT_MAX;

            // need to split the cluster further...
            // this could use meshopt but we're trying to get a complete baseline from metis
            if (result[i].indices.size() > kClusterSize * 3)
            {
                std::vector<Cluster> splits = ClusterizeMetis(remap, result[i].indices);
                assert(splits.size() > 1);

                result[i] = splits[0];
                for (size_t j = 1; j < splits.size(); ++j)
                    result.push_back(splits[j]);
            }
        }

        return result;
    }

    static std::vector<Cluster> Clusterize(const std::vector<Vertex>& vertices,
                                           const std::vector<uint32_t>& indices) {
        const size_t max_vertices = 192;
        const size_t max_triangles = kClusterSize;
        const size_t min_triangles = (kClusterSize / 3) & ~3;
        const float split_factor = 0.0f;
    
        size_t max_meshlets = meshopt_buildMeshletsBound(indices.size(), max_vertices, min_triangles);
    
        std::vector<meshopt_Meshlet> meshlets(max_meshlets);
        std::vector<unsigned int> meshlet_vertices(max_meshlets * max_vertices);
        std::vector<unsigned char> meshlet_triangles(max_meshlets * max_triangles * 3);
    
        // meshlets.resize(
        //     meshopt_buildMeshletsFlex(&meshlets[0], &meshlet_vertices[0], &meshlet_triangles[0],
        //         &indices[0], indices.size(), &vertices[0].position.x, vertices.size(), sizeof(Vertex),
        //         max_vertices, min_triangles, max_triangles, 0.f, split_factor
        //     )
        // );
        meshlets.resize(
            meshopt_buildMeshlets(&meshlets[0], &meshlet_vertices[0], &meshlet_triangles[0],
                &indices[0], indices.size(), &vertices[0].position.x, vertices.size(), sizeof(Vertex),
                max_vertices, max_triangles, 0.f)
        );
    
        std::vector<Cluster> clusters(meshlets.size());
    
        for (size_t i = 0; i < meshlets.size(); ++i) {
            const meshopt_Meshlet& meshlet = meshlets[i];
    
            meshopt_optimizeMeshlet(&meshlet_vertices[meshlet.vertex_offset], &meshlet_triangles[meshlet.triangle_offset], meshlet.triangle_count, meshlet.vertex_count);
    
            // note: for now we discard meshlet-local indices; they are valuable for shader code so in the future we should bring them back
            clusters[i].indices.resize(meshlet.triangle_count * 3);
            for (size_t j = 0; j < meshlet.triangle_count * 3; ++j)
                clusters[i].indices[j] = meshlet_vertices[meshlet.vertex_offset + meshlet_triangles[meshlet.triangle_offset + j]];
    
            clusters[i].parent.error = FLT_MAX;
        }
    
        return clusters;
    }

    static LODBounds Bounds(const std::vector<Vertex>& vertices, const std::vector<unsigned int>& indices, float error) {
        meshopt_Bounds bounds = meshopt_computeClusterBounds(&indices[0], indices.size(), &vertices[0].position.x, vertices.size(), sizeof(Vertex));

        LODBounds result;
        result.center[0] = bounds.center[0];
        result.center[1] = bounds.center[1];
        result.center[2] = bounds.center[2];
        result.radius = bounds.radius;
        result.error = error;
        return result;
    }

    static LODBounds BoundsMerge(const std::vector<Cluster>& clusters, const std::vector<int>& group) {
        std::vector<LODBounds> bounds(group.size());
        for (size_t j = 0; j < group.size(); ++j)
            bounds[j] = clusters[group[j]].self;

        meshopt_Bounds merged = meshopt_computeSphereBounds(&bounds[0].center[0], bounds.size(), sizeof(LODBounds), &bounds[0].radius, sizeof(LODBounds));

        LODBounds result = {};
        result.center[0] = merged.center[0];
        result.center[1] = merged.center[1];
        result.center[2] = merged.center[2];
        result.radius = merged.radius;

        // merged bounds error must be conservative wrt cluster errors
        result.error = 0.f;
        for (size_t j = 0; j < group.size(); ++j)
            result.error = std::max(result.error, clusters[group[j]].self.error);

        return result;
    }

    static void LockBoundary(std::vector<unsigned char>& locks,
        const std::vector<std::vector<int> >& groups,
        const std::vector<Cluster>& clusters) {
        // for each remapped vertex, keep track of index of the group it's in (or -2 if it's in multiple groups)
        std::vector<int> groupmap(locks.size(), -1);

        for (size_t i = 0; i < groups.size(); ++i)
            for (size_t j = 0; j < groups[i].size(); ++j)
            {
                const Cluster& cluster = clusters[groups[i][j]];

                for (size_t k = 0; k < cluster.indices.size(); ++k)
                {
                    unsigned int r = cluster.indices[k];

                    if (groupmap[r] == -1 || groupmap[r] == int(i))
                        groupmap[r] = int(i);
                    else
                        groupmap[r] = -2;
                }
            }

        // note: we need to consistently lock all vertices with the same position to avoid holes
        for (size_t i = 0; i < locks.size(); ++i)
        {
            locks[i] = (groupmap[i] == -2);
        }
    }

    static float boundsError(
        const LODBounds& bounds,
        const glm::mat4 worldPosMatrix,
        const glm::vec3& camera_pos,
        float camera_proj, float camera_znear) {
        glm::vec3 bounds_pos = glm::vec3(worldPosMatrix * glm::vec4(bounds.center[0], bounds.center[1], bounds.center[2], 1.0f));
        glm::vec3 bounds_dir = bounds_pos - camera_pos;
        float dis = glm::length(bounds_dir) - bounds.radius;
        return bounds.error / (dis > camera_znear ? dis : camera_znear) * (camera_proj * 0.5f);
        // float dx = bounds.center[0] - camera_pos.x, dy = bounds.center[1] - camera_pos.y, dz = bounds.center[2] - camera_pos.z;
        // float d = sqrtf(dx * dx + dy * dy + dz * dz) - bounds.radius;
        // return bounds.error / (d > camera_znear ? d : camera_znear) * (camera_proj * 0.5f);
    }

    static std::vector<std::vector<int> > PartitionMetis(
        const std::vector<Cluster>& clusters,
        const std::vector<int>& pending,
        const std::vector<Vertex>& in_vertices) {
        std::vector<std::vector<int> > result;
        std::vector<std::vector<int> > vertices(in_vertices.size());

        // find all cluster that vertices belong to
        for (size_t i = 0; i < pending.size(); ++i)
        {
            const Cluster& cluster = clusters[pending[i]];

            for (size_t j = 0; j < cluster.indices.size(); ++j)
            {
                int v = cluster.indices[j];

                std::vector<int>& list = vertices[v];
                if (list.empty() || list.back() != int(i))
                    list.push_back(int(i));
            }
        }

        std::map<std::pair<int, int>, int> adjacency;

        for (size_t v = 0; v < vertices.size(); ++v)
        {
            const std::vector<int>& list = vertices[v];

            for (size_t i = 0; i < list.size(); ++i)
                for (size_t j = i + 1; j < list.size(); ++j)
                    adjacency[std::make_pair(std::min(list[i], list[j]), std::max(list[i], list[j]))]++;
        }

        std::vector<std::vector<std::pair<int, int> > > neighbors(pending.size());

        // record adjacency weitht
        for (std::map<std::pair<int, int>, int>::iterator it = adjacency.begin(); it != adjacency.end(); ++it)
        {
            neighbors[it->first.first].push_back(std::make_pair(it->first.second, it->second));
            neighbors[it->first.second].push_back(std::make_pair(it->first.first, it->second));
        }

        std::vector<int> xadj(pending.size() + 1);
        xadj[0] = 0;
        std::vector<int> adjncy;
        std::vector<int> adjwgt;
        std::vector<int> part(pending.size());

        for (size_t i = 0; i < pending.size(); ++i)
        {
            for (size_t j = 0; j < neighbors[i].size(); ++j)
            {
                adjncy.push_back(neighbors[i][j].first);
                adjwgt.push_back(neighbors[i][j].second);
            }

            xadj[i + 1] = int(adjncy.size());
        }

        int options[METIS_NOPTIONS];
        METIS_SetDefaultOptions(options);
        options[METIS_OPTION_SEED] = 42;
        options[METIS_OPTION_UFACTOR] = 100;

        int nvtxs = int(pending.size());
        int ncon = 1;
        int nparts = int(pending.size() + kGroupSize - 1) / kGroupSize;
        // if (nparts == 1) {
        //     nparts = int(pending.size() + 4 - 1) / 4;
        //     if (nparts == 1) {
        //         nparts = int(pending.size() + 2 - 1) / 2;
        //     }
        // }
        int edgecut = 0;

        // not sure why this is a special case that we need to handle but okay metis
        if (nparts > 1)
        {
            int r = METIS_PartGraphRecursive(&nvtxs, (int *)&ncon, &xadj[0], &adjncy[0], NULL, NULL, &adjwgt[0], &nparts, NULL, NULL, options, &edgecut, &part[0]);
            assert(r == METIS_OK);
            (void)r;
        }

        result.resize(nparts);
        for (size_t i = 0; i < part.size(); ++i)
            result[part[i]].push_back(pending[i]);

        return result;
    }


    static std::vector<std::vector<int>> Partition(
        const std::vector<Cluster>& clusters,
        const std::vector<int>& pending,
        const std::vector<Vertex>& vertices) {

        std::vector<unsigned int> cluster_indices;
        std::vector<unsigned int> cluster_counts(pending.size());

        size_t total_index_count = 0;
        for (size_t i = 0; i < pending.size(); ++i)
            total_index_count += clusters[pending[i]].indices.size();

        cluster_indices.reserve(total_index_count);

        for (size_t i = 0; i < pending.size(); ++i)
        {
            const Cluster& cluster = clusters[pending[i]];

            cluster_counts[i] = unsigned(cluster.indices.size());

            for (size_t j = 0; j < cluster.indices.size(); ++j)
                cluster_indices.push_back(cluster.indices[j]); // 为了找cluster间公用的vertex
        }

        std::vector<unsigned int> cluster_part(pending.size());
        size_t partition_count = meshopt_partitionClusters(&cluster_part[0], &cluster_indices[0], cluster_indices.size(), &cluster_counts[0], cluster_counts.size(), vertices.size(), kGroupSize);

        std::vector<std::vector<int> > partitions(partition_count);
        for (size_t i = 0; i < partition_count; ++i)
            partitions[i].reserve(kGroupSize + 4);

        for (size_t i = 0; i < pending.size(); ++i)
            partitions[cluster_part[i]].push_back(pending[i]);

        return partitions;
    }

    static std::vector<unsigned int> Simplify(
        const std::vector<Vertex>& vertices,
        const std::vector<unsigned int>& indices,
        const std::vector<unsigned char>* locks,
        size_t target_count, float target_error, float* error = NULL) {
        if (target_count > indices.size()) {
            std::cout << "target_count > indices.size()" << std::endl;
            return indices;
        }

        std::vector<unsigned int> lod(indices.size());
        unsigned int options = 0;//meshopt_SimplifySparse/* | meshopt_SimplifyErrorAbsolute*/;
        //float normal_weights[3] = {0.5f, 0.5f, 0.5f};
        // if (kUseNormals)
        //     lod.resize(meshopt_simplifyWithAttributes(&lod[0], &indices[0], indices.size(), &vertices[0].px, vertices.size(), sizeof(Vertex), &vertices[0].nx, sizeof(Vertex), normal_weights, 3, locks ? &(*locks)[0] : NULL, target_count, FLT_MAX, options, error));
        lod.resize(meshopt_simplifySloppy(&lod[0], &indices[0], indices.size(), &vertices[0].position.x, vertices.size(), sizeof(Vertex), target_count, FLT_MAX, error));
        // if (locks != nullptr)
        //     lod.resize(meshopt_simplifyWithAttributes(&lod[0], &indices[0], indices.size(), &vertices[0].position.x, vertices.size(), sizeof(Vertex), NULL, 0, NULL, 0, &(*locks)[0], target_count, target_error, options, error));
        // else
        //     lod.resize(meshopt_simplify(&lod[0], &indices[0], indices.size(), &vertices[0].position.x, vertices.size(), sizeof(Vertex), target_count, target_error, options | meshopt_SimplifyLockBorder, error));
        return lod;
    }
};