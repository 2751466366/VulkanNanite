// vertex_shader.glsl
#version 450

layout(location = 0) in vec3 inPosition;
layout(location = 1) in vec3 inColor;

layout(location = 0) out vec3 fragColor;

layout(set = 0, binding = 0) uniform UniformBufferObject {
    mat4 model;
    mat4 view;
    mat4 proj;
} ubo;

layout(set = 1, binding = 0) uniform clusterInfo {
	uint clusterId;
} cluster;

uint MurmurMix(uint Hash){
	Hash^=Hash>>16;
	Hash*=0x45d9f3b;
	Hash^=Hash>>13;
	Hash*=0xc2b2ae35;
	Hash^=Hash>>12;
	return Hash;
}

vec3 to_color(uint idx)
{
	uint Hash = MurmurMix(idx + 1);

	vec3 color=vec3(
		(Hash>>0)&255,
		(Hash>>8)&255,
		(Hash>>16)&255
	);

	return color*(1.0f/255.0f);
}

void main() {
    gl_Position = ubo.proj * ubo.view * ubo.model * vec4(inPosition, 1.0);
    fragColor = to_color(cluster.clusterId);
}