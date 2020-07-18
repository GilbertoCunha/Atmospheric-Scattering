# The Physics of the color of the sky: Atmospheric Light Scattering
This repository is a project that was made for a Computational Physics class at University of Minho. It is an algorithm that calculates the color of the sky using light scattering physics.

In this repo you will be able to learn the physics behind it and how to implement it in Python code, producing some beautiful results.

You can calculate how the sunset looks:
![Sky](Images/sky.png "Sky color calculation")

Or maybe you want to look at the sky at around midday:
![Sky](Images/NewSky.png "Sky color calculation")

Or any other direction you want. You could even play around with some values to test how the sky would look like in other planets.

# How to navigate through the code
I subdivided the code in 3 files:

1. Vector.py: A file that contain a Vector class with some useful properties to make the development of the main algorithm simpler
2. Algorithm.py: Where the main algorithm is. It contains a function that calculates the intensity of a certain color for one pixel and some auxiliary functions to help out in this process.
3. MakeImage.py: This is where you will be able to generate the images. It calculates the color of each pixel for each RGB channel producing an image at the end. It also parallelizes the code of Algorithm.py to make calculations faster, even though they will still take some time.

I am also working on a Jupyter Notebook that thoroughly explains the whole project.

# How to generate images
First go to the MakeImage.py file. In there you will see many variables. I will not explain each and every one of them, just the main ones to help you tinker around and produce some nice results.

- 
