# Atmospheric Scattering

This is an algorithm for calculating the color of the sky in python. It takes in account both rayleigh and mie scattering. 
The code is subdivided in three chunks:
- Classes.py: Where I create a vector class with vector properties that will be useful throughout the code
- Algorithm.py: Where I create the algorithm that we will be using for calculating the color of one pixel
- Results.py: Where I generate a matrix of pixels (the image grid) and calculate an image of the sky according to the observer's viewpoint and the direction of the sun.
