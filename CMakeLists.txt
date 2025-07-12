cmake_minimum_required(VERSION 3.10)
project(VulkanTriangle)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 设置本地Vulkan和GLFW路径
set(VULKAN_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/opensource/vulkan")
#set(VULKAN_LIBRARIES "${CMAKE_CURRENT_SOURCE_DIR}/lib/vulkan-1.lib")

set(GLFW_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/opensource/glfw")
#set(GLFW_LIBRARIES "${CMAKE_CURRENT_SOURCE_DIR}/lib/glfw3.lib")

file(GLOB LIBRARIES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/lib/*.lib)
message(STATUS "lib文件列表: ${LIBRARIES}")

# 包含目录
include_directories(
    ${VULKAN_INCLUDE_DIRS}
    ${GLFW_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/opensource
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

file(GLOB sources CONFIGURE_DEPENDS src/*.cpp)
#message(STATUS "源文件列表: ${sources}")
# 添加源文件
add_executable(VulkanTriangle ${sources})

# 链接库
#target_link_libraries(VulkanTriangle ${VULKAN_LIBRARIES} ${GLFW_LIBRARIES})
target_link_libraries(VulkanTriangle ${LIBRARIES})