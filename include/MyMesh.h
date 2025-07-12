#pragma once
#include "MyDataType.h"

class MyMesh {
public:
    std::vector<Vertex> vertices;
    std::vector<uint32_t> indices;
    std::vector<Cluster> clusters;
    MyMesh() = default;
    void LoadCube() {
        vertices = {
            {{-0.5f, -0.5f,  0.5f}, {1.0f, 0.0f, 0.0f}},
            {{ 0.5f, -0.5f,  0.5f}, {0.0f, 1.0f, 0.0f}},
            {{ 0.5f,  0.5f,  0.5f}, {0.0f, 0.0f, 1.0f}},
            {{-0.5f,  0.5f,  0.5f}, {1.0f, 1.0f, 0.0f}},
            {{-0.5f, -0.5f, -0.5f}, {1.0f, 0.0f, 1.0f}},
            {{ 0.5f, -0.5f, -0.5f}, {0.0f, 1.0f, 1.0f}},
            {{ 0.5f,  0.5f, -0.5f}, {0.5f, 0.5f, 0.5f}},
            {{-0.5f,  0.5f, -0.5f}, {0.5f, 0.0f, 0.5f}}
        };
        indices = {
            // 正面
            0, 1, 2, 2, 3, 0,
            // 背面
            4, 7, 6, 6, 5, 4,
            // 左侧
            0, 3, 7, 7, 4, 0,
            // 右侧
            1, 5, 6, 6, 2, 1,
            // 顶部
            3, 2, 6, 6, 7, 3,
            // 底部
            0, 4, 5, 5, 1, 0
        };
    }
};