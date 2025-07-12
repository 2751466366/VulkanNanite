#pragma once
#include "MyDataType.h"
#include <unordered_map>
#include <queue>
#include <set>
#include <iomanip>
#include <iostream>

class MeshSimpler {
public:
    // 边结构体，用于优先级队列
    struct Edge {
        uint32_t v1, v2;        // 边的两个顶点
        float error;            // 折叠误差
        glm::vec3 newPos;
        uint32_t vnew;       // 折叠后的新顶点位置
        
        bool operator>(const Edge& other) const {
            return error > other.error;
        }
    };

    struct Triangle {
        uint32_t v1, v2, v3;    // 三角形的三个顶点
    };

    // 顶点到三角形的映射结构
    struct VertexTriangles {
        std::vector<Triangle> triangleIndices;  // 存储包含该顶点的三角形索引
    };

    struct Triangle2 {
        uint32_t v1, v2, v3;    // 三角形的三个顶点

        bool operator==(const Triangle2& other) const {
            return v1 == other.v1 && v2 == other.v2 && v3 == other.v3;
        }
        bool operator<(const Triangle2& other) const {
            return v1 < other.v1 || (v1 == other.v1 && (v2 < other.v2 || (v2 == other.v2 && v3 < other.v3)));
        }
        bool operator>(const Triangle2& other) const {
            return v1 > other.v1 || (v1 == other.v1 && (v2 > other.v2 || (v2 == other.v2 && v3 > other.v3)));
        }
        bool operator!=(const Triangle2& other) const {
            return !(*this == other);
        }
        Triangle2(uint32_t v1, uint32_t v2, uint32_t v3) {
            uint32_t min1 = std::min({v1, v2, v3});
            uint32_t max1 = std::max({v1, v2, v3});
            uint32_t mid1 = v1 + v2 + v3 - min1 - max1;
            this->v1 = min1;
            this->v2 = mid1;
            this->v3 = max1;
        }
    };
    struct VertexTriangles2 {
        std::set<Triangle2> triangleIndices;  // 存储包含该顶点的三角形索引
    };

    // 计算QEM矩阵
    static glm::mat4 computeQuadric(const glm::vec3& normal, const glm::vec3& point) {
        float d = -glm::dot(normal, point);
        return glm::mat4(
            normal.x * normal.x, normal.x * normal.y, normal.x * normal.z, normal.x * d,
            normal.x * normal.y, normal.y * normal.y, normal.y * normal.z, normal.y * d,
            normal.x * normal.z, normal.y * normal.z, normal.z * normal.z, normal.z * d,
            normal.x * d, normal.y * d, normal.z * d, d * d
        );
    }

    // 计算边折叠的误差和新顶点位置
    static std::pair<glm::vec3, float> computeEdgeCollapse(
        const glm::mat4& Q1, const glm::mat4& Q2,
        const glm::vec3& v1, const glm::vec3& v2) {
        
        glm::mat4 Q = Q1 + Q2;
        glm::mat4 Q4x4 = glm::mat4(
            Q[0][0], Q[0][1], Q[0][2], Q[0][3],
            Q[1][0], Q[1][1], Q[1][2], Q[1][3],
            Q[2][0], Q[2][1], Q[2][2], Q[2][3],
            0.0f, 0.0f, 0.0f, 1.0f
        );
        
        glm::vec3 newPos;
        float error;
        
        // if (glm::determinant(Q4x4) > 1e-6f) {
        //     glm::vec4 v = glm::inverse(Q4x4) * glm::vec4(0, 0, 0, 1);
        //     newPos = glm::vec3(v.x, v.y, v.z) / v.w;
        // } else {
            newPos = (v1 + v2) * 0.5f;
        // }
        
        error = glm::dot(glm::vec4(newPos, 1.0f), Q * glm::vec4(newPos, 1.0f));
        if (error < 0.0f) error = -error;
        
        return {newPos, error};
    }

    // 检查顶点是否是边界顶点
    static bool isBoundaryVertex(uint32_t vertexId, const std::vector<unsigned char>& boundary) {
        return vertexId < boundary.size() && boundary[vertexId] == 1;
    }

    // 检查边是否可以折叠
    static bool canCollapseEdge(const Edge& edge, 
                              const std::vector<bool>& vertexRemoved,
                              const std::vector<unsigned char>& boundary) {
        return !vertexRemoved[edge.v1] && !vertexRemoved[edge.v2] && 
               !isBoundaryVertex(edge.v1, boundary) && !isBoundaryVertex(edge.v2, boundary);
    }

    // 确保容器有足够空间
    static void ensureContainerCapacity(size_t newSize, 
                                      std::vector<glm::mat4>& vertexQuadrics,
                                      std::vector<bool>& vertexRemoved,
                                      std::vector<uint32_t>& vertexMap) {
        if (vertexQuadrics.size() <= newSize) {
            vertexQuadrics.resize(newSize + 1);
        }
        if (vertexRemoved.size() <= newSize) {
            vertexRemoved.resize(newSize + 1);
        }
        if (vertexMap.size() <= newSize) {
            vertexMap.resize(newSize + 1);
        }
    }

    static float Simplify(
        const std::vector<Vertex>& vertices, 
        const std::vector<uint32_t>& indices,
        const std::vector<unsigned char>& boundary,
        std::vector<Vertex>& outVertices,
        std::vector<uint32_t>& outIndices) {
        
        // 初始化输出顶点数组
        outVertices = vertices;
        
        // 1. 计算每个顶点的QEM矩阵并建立顶点到三角形的映射
        std::vector<glm::mat4> vertexQuadrics;
        std::vector<VertexTriangles> vertexTriangles(vertices.size() * 2);
        vertexQuadrics.reserve(vertices.size() * 2);
        vertexQuadrics.resize(vertices.size(), glm::mat4(0.0f));
        
        for (size_t i = 0; i < indices.size(); i += 3) {
            glm::vec3 v1 = outVertices[indices[i]].position;
            glm::vec3 v2 = outVertices[indices[i + 1]].position;
            glm::vec3 v3 = outVertices[indices[i + 2]].position;
            
            glm::vec3 normal = glm::normalize(glm::cross(v2 - v1, v3 - v1));
            glm::mat4 Kp = computeQuadric(normal, v1);
            
            // 更新QEM矩阵
            vertexQuadrics[indices[i]] += Kp;
            vertexQuadrics[indices[i + 1]] += Kp;
            vertexQuadrics[indices[i + 2]] += Kp;

            vertexTriangles[indices[i + 1]].
            triangleIndices.
            push_back(Triangle{indices[i], indices[i + 1], indices[i + 2]});
            vertexTriangles[indices[i + 2]].
            triangleIndices.
            push_back(Triangle{indices[i], indices[i + 1], indices[i + 2]});
            vertexTriangles[indices[i]].
            triangleIndices.
            push_back(Triangle{indices[i], indices[i + 1], indices[i + 2]});

        }

        // 2. 构建边折叠优先级队列
        std::priority_queue<Edge, std::vector<Edge>, std::greater<Edge>> edgeQueue;
        std::set<std::pair<uint32_t, uint32_t>> processedEdges;
        int errorCount = 0;
        
        // 遍历所有面片，构建边
        for (size_t i = 0; i < indices.size(); i += 3) {
            for (int j = 0; j < 3; j++) {
                uint32_t v1 = indices[i + j];
                uint32_t v2 = indices[i + (j + 1) % 3];
                
                if (v1 > v2) std::swap(v1, v2);
                if (processedEdges.find({v1, v2}) != processedEdges.end()) continue;
                processedEdges.insert({v1, v2});
                
                if (isBoundaryVertex(v1, boundary) || isBoundaryVertex(v2, boundary)) continue;
                
                auto [newPos, error] = computeEdgeCollapse(
                    vertexQuadrics[v1], vertexQuadrics[v2],
                    outVertices[v1].position, outVertices[v2].position
                );
                
                Edge edge{v1, v2, error, newPos};
                edgeQueue.push(edge);
                errorCount++;
            }
        }

        // 3. 执行边折叠
        std::vector<bool> vertexRemoved;
        std::vector<uint32_t> vertexMap;
        vertexRemoved.reserve(vertices.size() * 2);
        vertexMap.reserve(vertices.size() * 2);
        vertexRemoved.resize(vertices.size(), false);
        vertexMap.resize(vertices.size());
        for (size_t i = 0; i < vertices.size(); i++) {
            vertexMap[i] = i;
        }
        
        float totalError = 0.0f;
        int collapseCount = 0;
        int targetCollapseCount = errorCount / 2;
        
        while (!edgeQueue.empty() && collapseCount < targetCollapseCount) {
            Edge edge = edgeQueue.top();
            edgeQueue.pop();

            if (!canCollapseEdge(edge, vertexRemoved, boundary)) continue;
            
            // 更新顶点位置
            Vertex newPos;
            newPos.position = edge.newPos;
            SetColor(newPos);
            outVertices.push_back(newPos);
            uint32_t newVertexId = outVertices.size() - 1;
            
            // 确保容器有足够空间
            ensureContainerCapacity(newVertexId, vertexQuadrics, vertexRemoved, vertexMap);
            
            // 更新顶点映射
            vertexMap[edge.v1] = newVertexId;
            vertexMap[edge.v2] = newVertexId;
            vertexRemoved[edge.v1] = true;
            vertexRemoved[edge.v2] = true;
            
            // 更新QEM矩阵
            vertexQuadrics[newVertexId] = vertexQuadrics[edge.v1] + vertexQuadrics[edge.v2];
            
            // // 更新与新顶点相关的边
            updateConnectedEdgesOptimized(indices, outVertices, vertexQuadrics, boundary, 
                                       vertexRemoved, vertexMap, newVertexId, edgeQueue,
                                       vertexTriangles[edge.v1], vertexTriangles[edge.v2]);
            
            totalError += edge.error;
            collapseCount++;
        }

        //std::cout << "extra vertex size =" << collapseCount << std::endl;

        // 4. 重建索引
        for (size_t i = 0; i < indices.size(); i += 3) {
            uint32_t v1 = vertexMap[indices[i]];
            uint32_t v2 = vertexMap[indices[i + 1]];
            uint32_t v3 = vertexMap[indices[i + 2]];

            uint32_t newVertexId = outVertices.size() - 1;
            
            if (v1 == v2 || v2 == v3 || v1 == v3) {
                continue;
            }
            
            outIndices.push_back(v1);
            outIndices.push_back(v2);
            outIndices.push_back(v3);
        }

        return collapseCount > 0 ? totalError / collapseCount : 0.0f;
    }

    // 优化后的更新相关边函数
    static void updateConnectedEdgesOptimized(
        const std::vector<uint32_t>& indices,
        const std::vector<Vertex>& vertices,
        const std::vector<glm::mat4>& vertexQuadrics,
        const std::vector<unsigned char>& boundary,
        const std::vector<bool>& vertexRemoved,
        const std::vector<uint32_t>& vertexMap,
        uint32_t newVertexId,
        std::priority_queue<Edge, std::vector<Edge>, std::greater<Edge>>& edgeQueue,
        const VertexTriangles& v1Triangles,
        const VertexTriangles& v2Triangles) {
        
        std::set<std::pair<uint32_t, uint32_t>> newEdges;
        
        // 使用预存的三角形信息收集相关边
        auto processTriangles = [&](const VertexTriangles& triangles) {
            for (const Triangle& tri : triangles.triangleIndices) {
                uint32_t v1 = vertexMap[tri.v1];
                uint32_t v2 = vertexMap[tri.v2];
                uint32_t v3 = vertexMap[tri.v3];
                if (v1 == v2 || v2 == v3 || v1 == v3) {
                    continue;
                }
                if (v1 > v2) std::swap(v1, v2);
                if (v2 > v3) std::swap(v2, v3);
                if (v1 > v3) std::swap(v1, v3);
                newEdges.insert({v1, v2});
                newEdges.insert({v2, v3});
                newEdges.insert({v1, v3});
            }
        };
        
        processTriangles(v1Triangles);
        processTriangles(v2Triangles);
        
        // 重新计算新边的误差值
        for (const auto& [v1, v2] : newEdges) {
            if (vertexRemoved[v1] || vertexRemoved[v2]) continue;
            if (isBoundaryVertex(v1, boundary) || isBoundaryVertex(v2, boundary)) continue;
            
            auto [newPos, error] = computeEdgeCollapse(
                vertexQuadrics[v1], vertexQuadrics[v2],
                vertices[v1].position, vertices[v2].position
            );
            
            Edge newEdge{v1, v2, error, newPos};
            edgeQueue.push(newEdge);
        }
    }

    // 该方法不产生新顶点，而是将边坍缩到误差最小的一个顶点上
    // 该方法在循环处理坍缩边时不会更新新边
    static float Simplify2(
        const std::vector<Vertex>& vertices, 
        const std::vector<uint32_t>& indices,
        const std::vector<unsigned char>& boundary,
        std::vector<uint32_t>& outIndices) {
        
        // 1. 计算每个顶点的QEM矩阵并建立顶点到三角形的映射
        std::vector<glm::mat4> vertexQuadrics;
        vertexQuadrics.resize(vertices.size(), glm::mat4(0.0f));
        std::vector<VertexTriangles> vertexTriangles(vertices.size());
        
        for (size_t i = 0; i < indices.size(); i += 3) {
            glm::vec3 v1 = vertices[indices[i]].position;
            glm::vec3 v2 = vertices[indices[i + 1]].position;
            glm::vec3 v3 = vertices[indices[i + 2]].position;
            
            glm::vec3 normal = glm::normalize(glm::cross(v2 - v1, v3 - v1));
            glm::mat4 Kp = computeQuadric(normal, v1);
            
            vertexQuadrics[indices[i]] += Kp;
            vertexQuadrics[indices[i + 1]] += Kp;
            vertexQuadrics[indices[i + 2]] += Kp;

            vertexTriangles[indices[i]].triangleIndices.push_back(Triangle{indices[i], indices[i + 1], indices[i + 2]});
            vertexTriangles[indices[i + 1]].triangleIndices.push_back(Triangle{indices[i], indices[i + 1], indices[i + 2]});
            vertexTriangles[indices[i + 2]].triangleIndices.push_back(Triangle{indices[i], indices[i + 1], indices[i + 2]});
        }

        // 2. 构建边折叠优先级队列
        std::priority_queue<Edge, std::vector<Edge>, std::greater<Edge>> edgeQueue;
        int errorCount = 0;
        
        for (size_t i = 0; i < indices.size(); i += 3) {
            for (int j = 0; j < 3; j++) {
                uint32_t v1 = indices[i + j];
                uint32_t v2 = indices[i + (j + 1) % 3];
                
                if (v1 > v2) std::swap(v1, v2);
                
                if (isBoundaryVertex(v1, boundary) || isBoundaryVertex(v2, boundary)) {
                    continue;
                }
                
                // 计算两个端点的误差
                float error1 = glm::dot(glm::vec4(vertices[v1].position, 1.0f), 
                                      vertexQuadrics[v2] * glm::vec4(vertices[v1].position, 1.0f));
                float error2 = glm::dot(glm::vec4(vertices[v2].position, 1.0f), 
                                      vertexQuadrics[v1] * glm::vec4(vertices[v2].position, 1.0f));

                error1 = error1 < 0.0f ? -error1 : error1;
                error2 = error2 < 0.0f ? -error2 : error2;
                
                // 选择误差较小的点作为保留点
                glm::vec3 newPos = error1 <= error2 ? vertices[v1].position : vertices[v2].position;
                uint32_t vnew = error1 <= error2 ? v1 : v2;
                float error = std::min(error1, error2);
                
                Edge edge{v1, v2, error, newPos, vnew};
                edgeQueue.push(edge);
                errorCount++;
            }
        }

        // 3. 执行边折叠
        std::vector<bool> vertexRemoved(vertices.size(), false);
        std::vector<uint32_t> vertexMap(vertices.size());
        for (size_t i = 0; i < vertices.size(); i++) {
            vertexMap[i] = i;
        }

        auto isExpired = [&](uint32_t v1, uint32_t v2) {
            return vertexRemoved[v1] || vertexRemoved[v2] || 
                   vertexMap[v1] == vertexMap[v2];
        };
        
        float totalError = 0.0f;
        int collapseCount = 0;
        int targetCollapseCount = errorCount / 3;
        uint32_t removeVertex = -1;
        uint32_t keepVertex = -1;
        while (!edgeQueue.empty() && collapseCount < targetCollapseCount) {
            Edge edge = edgeQueue.top();
            edgeQueue.pop();

            if (isExpired(edge.v1, edge.v2)) continue;

            removeVertex = edge.vnew == edge.v1 ? edge.v2 : edge.v1;
            keepVertex = edge.vnew == edge.v1 ? edge.v1 : edge.v2;
            
            // 更新顶点映射
            vertexMap[removeVertex] = keepVertex;
            vertexRemoved[removeVertex] = true;

            std::set<std::pair<uint32_t, uint32_t>> newEdges;

            for (const Triangle& tri : vertexTriangles[removeVertex].triangleIndices) {
                uint32_t v1 = vertexMap[tri.v1];
                uint32_t v2 = vertexMap[tri.v2];
                uint32_t v3 = vertexMap[tri.v3];
                if (v1 == v2 || v2 == v3 || v1 == v3) {
                    continue;
                }
                if (v1 > v2) std::swap(v1, v2);
                if (v2 > v3) std::swap(v2, v3);
                if (v1 > v3) std::swap(v1, v3);
                newEdges.insert({v1, v2});
                newEdges.insert({v2, v3});
                newEdges.insert({v1, v3});
            }
            
            // 重新计算新边的误差值
            for (const auto& [v1, v2] : newEdges) {
                if (isExpired(edge.v1, edge.v2)) continue;
                if (isBoundaryVertex(v1, boundary) || isBoundaryVertex(v2, boundary)) continue;
                
                 // 计算两个端点的误差
                float error1 = glm::dot(glm::vec4(vertices[v1].position, 1.0f), 
                                      vertexQuadrics[v2] * glm::vec4(vertices[v1].position, 1.0f));
                float error2 = glm::dot(glm::vec4(vertices[v2].position, 1.0f), 
                                      vertexQuadrics[v1] * glm::vec4(vertices[v2].position, 1.0f));

                error1 = error1 < 0.0f ? -error1 : error1;
                error2 = error2 < 0.0f ? -error2 : error2;
                
                // 选择误差较小的点作为保留点
                glm::vec3 newPos = error1 <= error2 ? vertices[v1].position : vertices[v2].position;
                uint32_t vnew = error1 <= error2 ? v1 : v2;
                float error = std::min(error1, error2);
                
                Edge newEdge{v1, v2, error, newPos};
                edgeQueue.push(newEdge);
            }

            
            totalError += edge.error;
            collapseCount++;
        }

        // 4. 重建索引
        for (size_t i = 0; i < indices.size(); i += 3) {
            uint32_t v1 = vertexMap[indices[i]];
            uint32_t v2 = vertexMap[indices[i + 1]];
            uint32_t v3 = vertexMap[indices[i + 2]];
            
            if (v1 == v2 || v2 == v3 || v1 == v3) {
                continue;
            }
            
            outIndices.push_back(v1);
            outIndices.push_back(v2);
            outIndices.push_back(v3);
        }

        return collapseCount > 0 ? totalError / collapseCount : 0.0f;
        //return totalError;
    }
};

