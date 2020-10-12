# The Physics of the color of the sky: Atmospheric Light Scattering
This repository is the final project that I developed for a 
Computational Physics class at University of Minho. It is a Python 
algorithm that calculates the color of the sky using light scattering 
physics, more specifically both Rayleigh and Mie scattering.

It can calculate how the sunset looks:
![Sky](Images/sky.png "Sky color calculation")

Or maybe you want to look at the sky at around midday:
![Sky](Images/NewSky.png "Sky color calculation")

Or any other direction you want. You can even change the values of the 
earth's atmosphere to see how it would look like in a particular scenario.

# How to navigate through the code
I subdivided the code in 4 files under the folder Code:

1. Vector.py: A file that contains a Vector class with some useful 
properties to make the development of the main algorithm simpler.
2. Algorithm.py: Where the main algorithm is. It contains a function 
that calculates the intensity of a certain color for one pixel and some
auxiliary functions to help out in this process.
3. Aux.py: This file contains an auxiliary function to calculate the 
shooting direction Vectors for each pixel in the image plane.
4. MakeImage.py: This is where you will be able to generate the images. 
It calculates the color of each pixel for each RGB channel producing 
an image at the end. Is is also a parallelization of the code of 
Algorithm.py to make calculations faster, even though they will still
take some time.

# The Vector class
The vector class is just that, a class that represents vectors. 
I also use it to represent directions as unit vectors. It supports 
both cartesian and spherical coordinates, since these two coordinate 
systems are very useful for the development of this project. To create 
a vector object, you simply do the following:

vector = Vector(first_coordinate,
                second_coordinate,
                third_coordinate,
                coordinate_type)
                
Coordinate type can be either 'cart' or 'sph' representing the two
supported coordinate systems.

For the cartesian system the three coordinates are, respectively, x, y 
and z. For the spherical coordinates they are the distance to the origin,
the azimuthal and the polar angle, respectively.

This class also includes some useful methods for vector operations.

# How to generate images
First go to the MakeImage.py file under the Code folder. In there you
will see many variables. These are the main ones that will help you 
play around and produce some nice results:

- direction: The viewing direction of the observer
- sun_direction: The direction of the sun
- fov: The field of view of the observer
- width and length: number of rows and columns of pixels, respectively.
Since for now the program only supports square images, you should give
both of them the same value

You can then produce images with two simples steps:

1. Modify desired variables in MakeImage.py
2. Run MakeImage.py

Since calculating large images is time
consuming, you usually should preview it 
first with, say, 100x100 pixels to see what
to expect and then produce the final image.

Note: If you want to play around with the position of the observer, 
you might have to change the exposure of the image to be able to 
see it clearly.
