# Description
This project implements a Vulkan-based nanite effect.

# Running
<img src="./asset/screenShot.gif"/>

## Features

- **Mesh Simplification**: Simplifies 3D meshes by collapsing edges based on QEM.
- **Boundary Preservation**: Ensures that boundary vertices are not collapsed to maintain the integrity of the mesh.

## Key Directories

- **`asset/`**: Contains 3D models and other assets (e.g., `bunny_10k.obj`).
- **`include/`**: Header files, including the main mesh simplification logic in `MeshSimpler.hpp`.
- **`shaders/`**: Vulkan shader programs.
- **`src/`**: Source code for the Vulkan application.

## Dependencies

The project relies on the following libraries:

- [Vulkan SDK](https://vulkan.lunarg.com/sdk/home)
- [GLFW](https://www.glfw.org/)
- [Assimp](https://github.com/assimp/assimp)

## Build Instructions

### Prerequisites

- **CMake**: Ensure CMake is installed and available in your system's PATH.
- **Vulkan SDK**: Install the Vulkan SDK and set up the environment variables.
- **Compiler**: A C++ compiler that supports C++17 or later.

### Steps

1, run build.bat

2, run run.bat