# RustSPH

## Description
This project is a Smooth Particle Hydrodynamics (SPH) fluid simulation implementation written in Rust and Bevy.
## Demo

Here is a demo of the final implementation.

[![Demo](http://img.youtube.com/vi/8ooJT1wxFFo/0.jpg)](https://www.youtube.com/watch?v=8ooJT1wxFFo "Video Title")

This demo was run on a Surface Laptop 4.

## Overview

As described above, Rust was chosen as the programming langauge for this project. I have been attempting to learn Rust for several years, but always felt that I wasn't quite getting it. This project was an attempt to work through that. Bevy was chosen with Rust, as it has a simple way of rendering spheres, and contains an ECS system. 

This project was written for CSC473 at the Univeristy of Victoria, taught by Brandon Haworth in spring 2023.

### Paper

[CSC473_SPHFinalPaper.pdf](https://github.com/Syyreign/RustSPH/files/11589816/CSC473_SPHFinalPaper.pdf)

### Particles

The particles of this simulation are simple entities spheres controlled by a fluid system. Since this simulation runs primarily on the CPU, the amount of particles that can be rendered is rather limited.

### Collision

The collision is a custon point-plane collision as it was the simplest to implement. Unfortunatenly this means that the collision cannot have concave geometry, and in our case was a simple cube.

### UI

The UI is using a EGUI crate for Bevy. This crate allows for nice looking UI quickly. There was a concern that this UI would eat into the performance of the simulation, as EGUI is a immediate mode GUI. Fortunately, there was not a too noticable drop from tests performed.
https://github.com/mvlabat/bevy_egui

![image](https://github.com/Syyreign/RustSPH/assets/7028156/d4929214-5f66-4e8c-96d8-a81404aea512)

### Quirks
>* The collision can failed sometimes when the simulation is first started. This can be fixed by restarting the simulation. This is most likely a system ordering problem.
>* The particle quantity is limited, as this runs on the CPU. For my simulation computer n<700 particles was found to be stable.
>* As this simulation uses synplectic euler integration, the particles will begin jittering when the FPS gets too low. This could be fixed by using a better intergration methods, but wheres the fun in that.

## Future Work
As it stand this project is being shelved indefinitely. It was a fantastic learning tool, but with the amount of rewrite needed to improve this project, a new one would be more beneficial. In the future, it would be interesting to work on a PIC/FLIP fluid simulation and compare the two.
