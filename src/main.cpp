#define STB_IMAGE_IMPLEMENTATION
#include "GlfwGeneral.hpp"
#include "EasyVulkan.hpp"
#include "MyModel.hpp"
#include "MyView.h"
#include <ctime>
using namespace vulkan;

#define WIDHT 1920
#define HEIGHT 1080

pipelineLayout pipelineLayout_ubo;
pipeline pipeline_cube;
descriptorSetLayout descriptorSetLayout_view;
descriptorSetLayout descriptorSetLayout_cluster;
bool showGroup = false;
const auto& RenderPassAndFramebuffers() {
	static const auto& rpwf = easyVulkan::CreateRpwf_ScreenWithDS();
	return rpwf;
}
void CreateLayout() {
	VkDescriptorSetLayoutBinding descriptorSetLayoutBindingViewport = {
        .binding = 0,                                       //描述符被绑定到0号binding
        .descriptorType = VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER,
        .descriptorCount = 1,                               //个数是1个
        .stageFlags = VK_SHADER_STAGE_VERTEX_BIT            //在顶点着色器阶段读取uniform缓冲区
    };
	VkDescriptorSetLayoutCreateInfo descriptorSetLayoutCreateInfoView = {
        .bindingCount = 1,
        .pBindings = &descriptorSetLayoutBindingViewport
    };
	VkDescriptorSetLayoutBinding descriptorSetLayoutBindingCluster = {
        .binding = 0,                                       //描述符被绑定到0号binding
        .descriptorType = VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER_DYNAMIC,//类型为动态uniform缓冲区
        .descriptorCount = 1,                               //个数是1个
        .stageFlags = VK_SHADER_STAGE_VERTEX_BIT            //在顶点着色器阶段读取uniform缓冲区
    };
	VkDescriptorSetLayoutCreateInfo descriptorSetLayoutCreateInfoCluster = {
        .bindingCount = 1,
        .pBindings = &descriptorSetLayoutBindingCluster
    };
    descriptorSetLayout_view.Create(descriptorSetLayoutCreateInfoView);
	descriptorSetLayout_cluster.Create(descriptorSetLayoutCreateInfoCluster);
	VkDescriptorSetLayout descriptorSetLayouts[2] = {
		*descriptorSetLayout_view.Address(), *descriptorSetLayout_cluster.Address()};
    VkPipelineLayoutCreateInfo pipelineLayoutCreateInfo = {
        .setLayoutCount = 2,
        .pSetLayouts = descriptorSetLayouts
    };
	pipelineLayout_ubo.Create(pipelineLayoutCreateInfo);
}
void CreatePipeline() {
	static shaderModule vert("shaders/virtualMesh.vert.spv");
	static shaderModule frag("shaders/virtualMesh.frag.spv");
	static VkPipelineShaderStageCreateInfo shaderStageCreateInfos_triangle[2] = {
		vert.StageCreateInfo(VK_SHADER_STAGE_VERTEX_BIT),
		frag.StageCreateInfo(VK_SHADER_STAGE_FRAGMENT_BIT)
	};
	auto Create = [] {
		graphicsPipelineCreateInfoPack pipelineCiPack;
		pipelineCiPack.createInfo.layout = pipelineLayout_ubo;
		pipelineCiPack.createInfo.renderPass = RenderPassAndFramebuffers().renderPass;

		pipelineCiPack.vertexInputBindings.emplace_back(0, sizeof(Vertex), VK_VERTEX_INPUT_RATE_VERTEX);
		pipelineCiPack.vertexInputAttributes.emplace_back(0, 0, VK_FORMAT_R32G32B32_SFLOAT, offsetof(Vertex, position));
		pipelineCiPack.vertexInputAttributes.emplace_back(1, 0, VK_FORMAT_R32G32B32_SFLOAT, offsetof(Vertex, color));

		pipelineCiPack.inputAssemblyStateCi.topology = VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST;
		pipelineCiPack.viewports.emplace_back(0.f, 0.f, float(windowSize.width), float(windowSize.height), 0.f, 1.f);
		pipelineCiPack.scissors.emplace_back(VkOffset2D{}, windowSize);
		pipelineCiPack.rasterizationStateCi.cullMode = VK_CULL_MODE_BACK_BIT;
		pipelineCiPack.rasterizationStateCi.frontFace = VK_FRONT_FACE_COUNTER_CLOCKWISE;
		pipelineCiPack.rasterizationStateCi.lineWidth = 1.2;
		//pipelineCiPack.rasterizationStateCi.polygonMode = VK_POLYGON_MODE_LINE;
		pipelineCiPack.depthStencilStateCi.depthTestEnable = VK_TRUE;
		pipelineCiPack.depthStencilStateCi.depthWriteEnable = VK_TRUE;
		pipelineCiPack.depthStencilStateCi.depthCompareOp = VK_COMPARE_OP_LESS;
		pipelineCiPack.multisampleStateCi.rasterizationSamples = VK_SAMPLE_COUNT_1_BIT;
		pipelineCiPack.colorBlendAttachmentStates.push_back({ .colorWriteMask = 0b1111 });
		pipelineCiPack.UpdateAllArrays();
		pipelineCiPack.createInfo.stageCount = 2;
		pipelineCiPack.createInfo.pStages = shaderStageCreateInfos_triangle;
		pipeline_cube.Create(pipelineCiPack);
	};
	auto Destroy = [] {
		pipeline_cube.~pipeline();
	};
	graphicsBase::Base().AddCallback_CreateSwapchain(Create);
	graphicsBase::Base().AddCallback_DestroySwapchain(Destroy);
	Create();
}


MyView myView;
float deltaTime = 0.0f;	// time between current frame and last frame
float lastFrame = 0.0f;
bool firstMouse = true;
double lastX = WIDHT / 2.0f;
double lastY = HEIGHT / 2.0f;
float yaw = -90.0f; // yaw is initialized to -90.0 degrees, to rotate the camera to the right
float pitch = 0.0f; // pitch is initialized to 0.0 degrees, to keep the camera level
void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
    if (firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top
    lastX = xpos;
    lastY = ypos;

    float sensitivity = 0.05f; // change this value to your liking
    xoffset *= sensitivity;
    yoffset *= sensitivity;

    yaw += xoffset;
    pitch += yoffset;

    // make sure that when pitch is out of bounds, screen doesn't get flipped
    if (pitch > 89.0f)
        pitch = 89.0f;
    if (pitch < -89.0f)
        pitch = -89.0f;

    glm::vec3 front;
    front.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
    front.y = sin(glm::radians(pitch));
    front.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
    myView.cameraFront = glm::normalize(front);
}

void processInput(GLFWwindow *window, time_t time)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

	float currentFrame = glfwGetTime();
	deltaTime = currentFrame - lastFrame;
	lastFrame = currentFrame;
    float cameraSpeed = 5 * deltaTime;
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
		myView.cameraPos += cameraSpeed * myView.cameraFront;
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
		myView.cameraPos -= cameraSpeed * myView.cameraFront;
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
		myView.cameraPos -= glm::normalize(glm::cross(myView.cameraFront, myView.cameraUp)) * cameraSpeed;
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        myView.cameraPos += glm::normalize(glm::cross(myView.cameraFront, myView.cameraUp)) * cameraSpeed;
	if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS && std::time(nullptr) - time > 1.0f) {
		showGroup = !showGroup;
	}
}

int main() {
	if (!InitializeWindow({ WIDHT, HEIGHT }))
		return -1;
	glfwSetCursorPosCallback(pWindow, mouse_callback);
	glfwSetInputMode(pWindow, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

	const auto& [renderPass, framebuffers] = RenderPassAndFramebuffers();
	CreateLayout();
	CreatePipeline();

	myView.SetProj(windowSize.width, windowSize.height);
	MyModel myModel;
	myModel.LoadModel("asset/bunny_10k.obj");
	myModel.Nanite();


	uniformBuffer uniformBufferMat(sizeof(glm::mat4) * 3, VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);
	uniformBufferMat.TransferData(&myModel.worldPosMatrix, sizeof(glm::mat4), 0);
	uniformBufferMat.TransferData(&myView.view, sizeof(glm::mat4), sizeof(glm::mat4));
	uniformBufferMat.TransferData(&myView.projection, sizeof(glm::mat4), sizeof(glm::mat4) * 2);

	VkDeviceSize minUboAlignment = graphicsBase::Base().PhysicalDeviceProperties().limits.minUniformBufferOffsetAlignment;
	VkDeviceSize dynamicAlignment = (sizeof(uint32_t) + minUboAlignment - 1) & ~(minUboAlignment - 1);
	VkDeviceSize bufferSize = dynamicAlignment * myModel.clusters.size();
	uniformBuffer uniformBufferCluster(bufferSize, VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);
	// 更新所有cluster的uniform数据
	for (int i = 0; i < myModel.clusters.size(); i++) {
		VkDeviceSize offset = i * dynamicAlignment;
		uniformBufferCluster.TransferData(&myModel.clusters[i].clusterId, sizeof(uint32_t), offset);
	}

	vertexBuffer vertexBufferCube(myModel.vertices.size() * sizeof(Vertex));
	vertexBufferCube.TransferData(myModel.vertices.data(), myModel.vertices.size() * sizeof(Vertex));
	
	int size = 0;
	for (int i = 0; i < myModel.clusters.size(); i++) {
		size += myModel.clusters[i].indices.size();
	}
	indexBuffer indexBuffer(size * sizeof(uint32_t));
	int lastSize = 0;


	fence fence;
	semaphore semaphore_imageIsAvailable;
	semaphore semaphore_renderingIsOver;

	commandBuffer commandBuffer;
	commandPool commandPool(graphicsBase::Base().QueueFamilyIndex_Graphics(), VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT);
	commandPool.AllocateBuffers(commandBuffer);

	VkDescriptorPoolSize descriptorPoolSizes[] = {
		{ VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1 },
		{ VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER_DYNAMIC, 1 }
	};
	descriptorPool descriptorPool(2, descriptorPoolSizes);
	descriptorSet descriptorSet_view;
	descriptorSet descriptorSet_cluster;
	descriptorPool.AllocateSets(descriptorSet_view, descriptorSetLayout_view);
	descriptorPool.AllocateSets(descriptorSet_cluster, descriptorSetLayout_cluster);
	VkDescriptorBufferInfo bufferInfoMat = {
		.buffer = uniformBufferMat,
		.offset = 0,
		.range = VK_WHOLE_SIZE
	};
	VkDescriptorBufferInfo bufferInfoCluster = {
		.buffer = uniformBufferCluster,
		.offset = 0,
		.range = dynamicAlignment
	};
	descriptorSet_view.Write(bufferInfoMat, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 0);
	descriptorSet_cluster.Write(bufferInfoCluster, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER_DYNAMIC, 0);

	VkClearValue clearValues[2] = {
		{ .color = { 0.12f, 0.15f, 0.1f, 1.f } },
		{ .depthStencil = { 1.f, 0 } }
	};

	time_t time = std::time(nullptr);
	int g = 0;
	int idxSize = size;

	float threshold = 2e-4f; // 2 pixels at 1080p
	int drawCount = 0;

	while (!glfwWindowShouldClose(pWindow)) {
		while (glfwGetWindowAttrib(pWindow, GLFW_ICONIFIED))
			glfwWaitEvents();
		processInput(pWindow, time);

		drawCount = 0;
		std::vector<uint32_t> clusterIds;

		for (int i = 0; i < myModel.clusters.size(); i++) {
			float selferr = NaniteMgr::boundsError(
				myModel.clusters[i].self,
				myModel.worldPosMatrix,
				myView.cameraPos,
				-myView.projection[1][1],
				myView.znear);
			float parenterr = NaniteMgr::boundsError(
				myModel.clusters[i].parent,
				myModel.worldPosMatrix,
				myView.cameraPos,
				-myView.projection[1][1],
				myView.znear);
			if ((selferr < threshold && parenterr > threshold) ||
				(parenterr > threshold <= threshold && myModel.clusters[i].parent.error == FLT_MAX)) {
			// if (myModel.clusters[i].clusterId == 112) {
				indexBuffer.TransferData(
					myModel.clusters[i].indices.data(),
					myModel.clusters[i].indices.size() * sizeof(uint32_t),
					drawCount * sizeof(uint32_t));
				drawCount += myModel.clusters[i].indices.size();
				clusterIds.push_back(i);
			}
		}

		myView.updateView();
		// 更新view矩阵
		uniformBufferMat.TransferData(&myView.view, sizeof(glm::mat4), sizeof(glm::mat4));

		// 更新cluster的uniform数据
		for (int i = 0; i < myModel.clusters.size(); i++) {
			VkDeviceSize offset = i * dynamicAlignment;
			uniformBufferCluster.TransferData(&myModel.clusters[i].clusterId, sizeof(uint32_t), offset);
		}

		graphicsBase::Base().SwapImage(semaphore_imageIsAvailable);
		auto i = graphicsBase::Base().CurrentImageIndex();

		commandBuffer.Begin(VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT);
		renderPass.CmdBegin(commandBuffer, framebuffers[i], { {}, windowSize }, clearValues);

		VkDeviceSize offset = 0;
		vkCmdBindVertexBuffers(commandBuffer, 0, 1, vertexBufferCube.Address(), &offset);
		vkCmdBindIndexBuffer(commandBuffer, indexBuffer, 0, VK_INDEX_TYPE_UINT32);

		vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, pipeline_cube);

		vkCmdBindDescriptorSets(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS,
			pipelineLayout_ubo, 0, 1, descriptorSet_view.Address(), 0, nullptr);

		uint32_t indexOffset = 0;
		for (int c = 0; c < clusterIds.size(); c++) {
			auto& cluster = myModel.clusters[clusterIds[c]];
			uint32_t dynamicOffset = clusterIds[c] * dynamicAlignment;
			
			vkCmdBindDescriptorSets(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS,
				pipelineLayout_ubo, 1, 1, descriptorSet_cluster.Address(), 1, &dynamicOffset);

			vkCmdDrawIndexed(commandBuffer, cluster.indices.size(), 1, indexOffset, 0, 0);
			indexOffset += cluster.indices.size();
		}

		renderPass.CmdEnd(commandBuffer);
		commandBuffer.End();

		graphicsBase::Base().SubmitCommandBuffer_Graphics(commandBuffer, semaphore_imageIsAvailable, semaphore_renderingIsOver, fence);
		graphicsBase::Base().PresentImage(semaphore_renderingIsOver);

		glfwPollEvents();
		TitleFps();

		fence.WaitAndReset();
	}
	TerminateWindow();
	return 0;
}