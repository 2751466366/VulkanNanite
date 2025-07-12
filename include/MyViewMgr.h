#pragma once
#include "common.h"

class MyView {
public:
    glm::mat4 projection;
    glm::mat4 view;
    glm::vec3 cameraPos   = glm::vec3(0.0f, 0.0f, 3.0f);
    glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
    glm::vec3 cameraUp    = glm::vec3(0.0f, 1.0f, 0.0f);
    void SetProj(int width, int height) {
        projection = glm::perspective(45.f, float(width) / float(height), 0.1f, 100.f);
        //projection[1][1] *= -1; // Invert Y axis

        view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);
    }
    void updateView() {
        view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);
    }
    void setCameraPos(glm::vec3 pos) {
        cameraPos = pos;
    }
};