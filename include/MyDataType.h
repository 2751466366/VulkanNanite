#pragma once
#include "common.h"

struct Vertex {
    glm::vec3 position;
    glm::vec3 color;
};

struct LODBounds
{
	float center[3];
	float radius;
	float error;
};

struct Cluster
{
	std::vector<unsigned int> indices;

	LODBounds self;
	LODBounds parent;
	uint32_t clusterId;
};

void SetColor(Vertex& vertex) {
	static int color = 0;
	vertex.color.x = (float)(color == 0 ? 1.0f : 0.0f);
	vertex.color.y = (float)(color == 1 ? 1.0f : 0.0f);
	vertex.color.z = (float)(color == 2 ? 1.0f : 0.0f);
	if (color >= 2) {
		color = 0;
	} else {
		color++;
	}
}

void SetLodColor(Vertex& vertex, int lod, int maxlod) {
	static int color = 0;
	float k = lod *  1.0f / maxlod;
	vertex.color.x = (float)(lod % 3 == 0 ? 1.0f * k : 1.0f * (1 - k)) * (color == 0 ? 1.0f : 0.5f);
	vertex.color.y = (float)(lod % 3 == 1 ? 1.0f * (1 - k) : 1.0f * k) * (color == 0 ? 1.0f : 0.5f);
	vertex.color.z = (float)(lod % 3 == 2 ? 1.0f * k: 1.0f * (1 - k)) * (color == 0 ? 1.0f : 0.5f);
	if (color >= 1) {
		color = 0;
	} else {
		color++;
	}
}

uint32_t MurmurMix(uint32_t Hash){
	Hash^=Hash>>16;
	Hash*=0x45d9f3b;
	Hash^=Hash>>13;
	Hash*=0xc2b2ae35;
	Hash^=Hash>>12;
	return Hash;
}