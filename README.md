# SomeGame
This project is very much in early development and has been embarked on primarily as a learning excersise
I am trying to implement a full 3d game engine in c++ from scratch, inspired partly by the Handmade Hero series on Youtube and Tiwtch.

## What I am trying to acomplish
* A working physics engine complete with 3d collision
* My own (performant) software renderer
* creating a very large, procedurally generated, world

## Learning Rescources I have found useful so far
* [Handmade Hero](https://handmadehero.org/) - This is probably the best resource for understanding how games are made
### Physics Simulation
* [Physically Based Modelling](http://www.cs.cmu.edu/~baraff/pbm/pbm.html) - Understanding physics simulation
* [Differential Equations](https://youtube.com/playlist?list=PLZHQObOWTQDNPOjrT6KVlfJuKtYTftqH6) - Most of the maths behind physics simulation is differential equations
* [Multivarible Calculus](https://youtube.com/playlist?list=PLSQl0a2vh4HC5feHa6Rc5c0wbRTx56nF7) - Also helpful to understand some of the Physically Based Modelling documents

### 3D Rendering
* [3D Graphics Pipeline](https://youtu.be/7qUuzRY5YwI) - Understanding the steps renderer has to go through to go to get an image on screen
* [3D Programming Fundamentals](https://youtube.com/playlist?list=PLqCJpWy5Fohe8ucwhksiv9hTF5sfid8lA) - More in depth on implemntation of 3D rendering
* [Quaternions and Rotation](https://graphics.stanford.edu/courses/cs348a-17-winter/Papers/quaternion.pdf) - Quaternions provide an alternative way to represent roatations in 3d space that is more effecient that using matrices

### Optimisation
* [Single Instruction Multiple Data (SIMD)](https://en.wikipedia.org/wiki/Single_instruction,_multiple_data) - Performance and optimisation
* [Modern x64 Archietures and the Cache](https://youtu.be/tk5P7mt2fAw) - It seems that in most cases ineffecient use of the cache is the cause of bad performance
* [Data Oriented Design](https://youtu.be/rX0ItVEVjHc) - Good practises for making good use of CPU cache
* [CPU Caches and Why You Care](https://youtu.be/WDIkqP4JbkE) - This talk a bit like the previous outlines the importance of CPU cache

I have also used wikipedia extensively and have found almost all the relevant wikipedia articles to be helpful
