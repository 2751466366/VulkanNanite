#pragma once
#include "assimp/Importer.hpp"
#include "assimp/scene.h"
#include "assimp/postprocess.h"
#include "MyDataType.h"
#include "NaniteMgr.hpp"
#include "MeshSimpler.hpp"

class MyModel {
public:
    //MyMesh myMesh;
    std::vector<Vertex> vertices;
    std::vector<uint32_t> indices;
    std::vector<Cluster> clusters;
    std::vector<std::vector<int>> groups;
    std::vector<std::vector<uint32_t>> simplifiedIndices;
    glm::mat4 worldPosMatrix = glm::mat4(1.f);
    MyModel() = default;
    MyModel(const std::string& path) {
        LoadModel(path);
    }
    void LoadModel(const std::string &path) {
        Assimp::Importer importer;

        // 导入STL模型
        const aiScene* scene = importer.ReadFile(path, aiProcess_Triangulate | aiProcess_GenNormals | aiProcess_FlipUVs);  // 翻转UV坐标（如果需要）

        // 检查是否导入成功
        if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) {
            std::cerr << "ERROR::ASSIMP::" << importer.GetErrorString() << std::endl;
            return;
        }

        for (unsigned int i = 0; i < scene->mNumMeshes; i++) {
            aiMesh* mesh = scene->mMeshes[i];

            int idx = 0;
            for (unsigned int j = 0; j < mesh->mNumVertices; j++) {
                Vertex vertex;
                vertex.position.x = mesh->mVertices[j].x;
                vertex.position.y = mesh->mVertices[j].y;
                vertex.position.z = mesh->mVertices[j].z;
                
                vertices.push_back(vertex);
            }

            for (unsigned int j = 0; j < mesh->mNumFaces; j++) {
                aiFace face = mesh->mFaces[j];
                for (unsigned int k = 0; k < face.mNumIndices; k++) {
                    indices.push_back(face.mIndices[k]);
                }
            }
            std::cout << "meshType: " << mesh->mPrimitiveTypes << std::endl;
            std::cout << "meshName: " << mesh->mName.C_Str() << std::endl;
        }
    }

    void Nanite() {
        NaniteMgr::RemapMesh(vertices, indices);

        // for (int i = 0; i < vertices.size(); i++) {
        //     SetColor(vertices[i]);
        // }

        std::cout << "vertices: " << vertices.size() << std::endl;
        std::cout << "indices: " << indices.size() << std::endl;
        clusters = NaniteMgr::Clusterize(vertices, indices);
        std::cout << "clusters: " << clusters.size() << std::endl;

        // for (const auto& cluster : clusters) {
        //     Vertex temp;
        //     SetColor(temp);
        //     for (const auto& index : cluster.indices) {
        //         vertices[index].color = temp.color;
        //     }
        // }

        for (size_t i = 0; i < clusters.size(); ++i) {
            clusters[i].self = NaniteMgr::Bounds(vertices, clusters[i].indices, 0.f);
        }

        std::vector<int> pending(clusters.size());
        for (size_t i = 0; i < clusters.size(); ++i)
            pending[i] = int(i);

        // 记录边界顶点，避免简化时被删除

        while (pending.size() > 1) {
            groups = NaniteMgr::PartitionMetis(clusters, pending, vertices);
            //groups = NaniteMgr::Partition(clusters, pending, vertices);
            if (groups.size() <= 1) {
                break;
            }
            std::vector<unsigned char> locks(vertices.size());
            NaniteMgr::LockBoundary(locks, groups, clusters);

            pending.clear();
    
            size_t triangles = 0;
            size_t stuck_triangles = 0;
    
            // every group needs to be simplified now
            for (size_t i = 0; i < groups.size(); ++i) {
                if (groups[i].empty())
                    continue; // metis shortcut
    
                std::vector<unsigned int> merged;
                for (size_t j = 0; j < groups[i].size(); ++j)
                    merged.insert(merged.end(), clusters[groups[i][j]].indices.begin(), clusters[groups[i][j]].indices.end());
    

                std::vector<uint32_t> outIndices;
                //float error = MeshSimpler::Simplify2(vertices, merged, locks, outIndices);
                std::vector<Vertex> outVertices;
                float error = MeshSimpler::Simplify(vertices, merged, locks, outVertices, outIndices);
                vertices = std::move(outVertices);
                simplifiedIndices.push_back(outIndices);
    

                LODBounds groupb = NaniteMgr::BoundsMerge(clusters, groups[i]);
                groupb.error += error;
    
                std::vector<Cluster> split = NaniteMgr::Clusterize(vertices, outIndices);
    
                for (size_t j = 0; j < groups[i].size(); ++j)
                {
                    assert(clusters[groups[i][j]].parent.error == FLT_MAX);
                    clusters[groups[i][j]].parent = groupb;
                }
    
                for (size_t j = 0; j < split.size(); ++j)
                {
                    split[j].self = groupb;

    
                    clusters.push_back(split[j]); // std::move
                    pending.push_back(int(clusters.size()) - 1);
                }
            }
        }
        for (int i = 0; i < clusters.size(); i++) {
            clusters[i].clusterId = i;
        }
        std::cout << "clusters after: " << clusters.size() << std::endl;
        std::cout << "vertices after: " << vertices.size() << std::endl;
    }

    void SetWorldPos(glm::vec3 pos = glm::vec3( 0.0f,  0.0f,  0.0f)) {
        worldPosMatrix = glm::translate(worldPosMatrix, pos);
    }
};