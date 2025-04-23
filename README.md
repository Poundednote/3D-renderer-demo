# 3D Renderer Demo
I originally had planned this to be a full 3d game engine however I want to put this project on hold here while I focus on other things. What I have managed to accomplish though, is a full 3d software rasteriser that has rudimentary mesh support. I primarily did this to better understand how computers handle 3d graphics. I also have implemented a procedurally generated world, with a physics engine. The world models planets and suns where the planets orbit the suns. The suns use Area lighting with a radial linear falloff and have a bloom effect applied to them.



## Learning Resources I have found useful so far
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
