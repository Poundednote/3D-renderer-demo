# 3D Renderer Demo
Custom C++ Software Rasterizer & Physics Engine

A high-performance, scratch-built 3D game engine written in C++ without any graphics APIs (OpenGL, DirectX, or Vulkan). This project implements a Scanline-based Software Rasterizer to render 3D geometry directly to a memory buffer, coupled with a custom physics engine simulating orbital mechanics and spring dynamics.

The engine prioritizes Data-Oriented Design and SIMD optimizations to achieve real-time performance on the CPU.
Core Architecture
1. Scanline Rasterization Pipeline

Unlike modern GPU pipelines that process triangles in parallel, this engine implements a classic scanline algorithm to rasterize geometry.

    Triangle Decomposition: Triangles are sorted and split into Flat-Top and Flat-Bottom segments to simplify traversal.

    Edge Walking: The rasterizer walks down the left and right edges of the triangle, calculating gradients for attributes (depth, normals, color) per scanline.

    Perspective Correct Interpolation: Attributes are interpolated across the scanline to ensure texture/color accuracy at depth.

    Z-Buffering: Implements a depth buffer for hidden surface removal, with an Early Depth Test optimization to discard occluded pixels before expensive shading calculations.

2. Hand-Optimized SIMD Math

The core math library (particle_math.h) is built on SSE Intrinsics (__m128) to vectorize heavy arithmetic operations.

    4-Wide Vertex Processing: Transformations (World -> View -> Projection) are applied to 4 vertices simultaneously using SIMD registers.

    Quaternion Rotation: Rotations are handled via custom quaternion math, avoiding Gimbal lock and reducing matrix overhead.

3. Data-Oriented Physics System

The particle system uses a Structure of Arrays (SoA) layout instead of Array of Structures (AoS) to maximize CPU cache locality.

    Contiguous Streams: Position, velocity, and mass are stored in separate arrays (V3 pos[MAX], float mass[MAX]), allowing the physics solver to stream through memory linearly without cache pollution from unrelated data (like render flags).

    Simulation Features:

        N-Body Gravity: Simulates planetary orbits where bodies exert gravitational forces on each other.

        Spring Dynamics: Implements Hooke's Law with damping for cloth-like or elastic simulations.

4. Memory Management

    Zero-Allocation Runtime: The engine allocates two large blocks of memory at startup (Permanent and Transient).

    Linear Allocation: Per-frame data uses a bump-pointer allocator in the Transient block, which is reset instantly at the end of every frame, completely eliminating memory fragmentation and malloc/free overhead.

# Things I've learned since then
* GPUs don't do scanline rendering because its hard to parallelise. Instead they use edge detection algorithms to check if pixels are inside triangles. You can do these edge detection equations in parallel across a block of pixels.
* Quaternions are cool but harder to parallelise with simd vs matrices and column major matrices make vector matrix operations much faster in SIMD.
* You can bin triangles into a screen sector so that each CPU thread can work on a specific part of the screen when rasterizing pixels. 

# Learning Resources I have found useful so far
* [Handmade Hero](https://handmadehero.org/) - This is probably the best resource for understanding how games are made
### Physics Simulation
* [Physically Based Modelling](http://www.cs.cmu.edu/~baraff/pbm/pbm.html) - Understanding physics simulation
* [Differential Equations](https://youtube.com/playlist?list=PLZHQObOWTQDNPOjrT6KVlfJuKtYTftqH6) - Most of the maths behind physics simulation is differential equations
* [Multivarible Calculus](https://youtube.com/playlist?list=PLSQl0a2vh4HC5feHa6Rc5c0wbRTx56nF7) - Also helpful to understand some of the Physically Based Modelling documents

### 3D Rendering
* [3D Graphics Pipeline](https://youtu.be/7qUuzRY5YwI) - Understanding the steps the renderer has to go through to go to get an image on the screen
* [3D Programming Fundamentals](https://youtube.com/playlist?list=PLqCJpWy5Fohe8ucwhksiv9hTF5sfid8lA) - More in-depth on implementation of 3D rendering
* [Quaternions and Rotation](https://graphics.stanford.edu/courses/cs348a-17-winter/Papers/quaternion.pdf) - Quaternions provide an alternative way to represent rotations in 3d space that is more efficient that using matrices

### Optimisation
* [Single Instruction Multiple Data (SIMD)](https://en.wikipedia.org/wiki/Single_instruction,_multiple_data) - Performance and optimisation
* [Modern x64 Architecture and the Cache](https://youtu.be/tk5P7mt2fAw) - It seems that in most cases ineffecient use of the cache is the cause of bad performance
* [Data Oriented Design](https://youtu.be/rX0ItVEVjHc) - Good practises for making good use of CPU cache
* [CPU Caches and Why You Care](https://youtu.be/WDIkqP4JbkE) - This talk is a bit like the previous outlines on the importance of CPU cache

I have also used Wikipedia extensively and have found almost all the relevant Wikipedia articles to be helpful
