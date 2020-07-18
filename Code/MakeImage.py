# Import internal libraries
from Code.Algorithm import *
from Code.Vector import Vector
from Code.Aux import angle_matrix

# Import external libraries
import matplotlib.pyplot as plt
from itertools import product
import concurrent.futures
from tqdm import tqdm
import time

# Position and direction variables
position = Vector(earth_radius + 1.8, 0, 0, 'sph')
direction, sun_direction = Vector(1, 0, pi / 2, 'sph'), Vector(1, 0, pi / 2 - 0.05, 'sph')

# Light wavelengths
wavelengths = [(6.8e-7, 0.685), (5.35e-7, 0.81), (4.6e-7, 0.775)]

# pixel width and height, number of light beams shot per pixel, field of view
width, length, bpp, fov = 50, 50, 1, pi / 12

# Tone mapping exposure and gamma correction
exposure, gamma_corr = 1.8e4, 2.2

# Generating the shooting angles for each pixel
start = time.time()
shoot_angles = angle_matrix(width, length, fov, bpp, direction)
image = []


def do_pixel(arg):
    """
    This function calculates one pixel of the final image given its row and column
    It also applies tone mapping to give the image a more dynamic color range and to
    brighten the image, since these physical calculations tend to return small
    intensity values

    args -> (int, int):
        arg[0]: pixel row
        arg[1]: pixel column

    return -> [float, float, float]: RGB color of the pixel
    """

    row, column = int(arg[0]), int(arg[1])
    pixel = []
    for (w, sun_emit) in wavelengths:
        intensity = 0
        for k in range(bpp):
            intensity += main_cycle(position, shoot_angles[row][column][k], sun_direction, w)
        tone_map = (1 - np.exp(-sun_emit * intensity * exposure)) ** (1 / gamma_corr)
        pixel += [tone_map]

    return pixel


# Parallelize the calculation for each pixel
if __name__ == '__main__':
    with concurrent.futures.ProcessPoolExecutor() as executor:
        wid_list, len_list = np.linspace(0, width - 1, width), np.linspace(0, length - 1, length)

        args = tqdm(product(wid_list, len_list), total=width * length, desc="Pre-processing")

        results = list(tqdm(executor.map(do_pixel, args), total=width * length, desc="Calculating"))

        for result in results:
            image.append(result)


# Reshaping the result as a np array
image = np.asarray(image).reshape((width, length, 3))
stop = time.time()
print(f'This took {stop-start:.2f} seconds | {(stop-start)/60:.2f} minutes | {(stop-start)/3600:.2f} hours')

# Plotting the final image
fig, ax = plt.subplots(dpi=400)
ax.imshow(image)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_yticks([])
ax.set_xticks([])
plt.savefig('sky1T.png', dpi=800)
plt.show()
