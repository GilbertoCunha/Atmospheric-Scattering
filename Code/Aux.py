import numpy as np
from tqdm import tqdm
from Code.Vector import Vector, Transformation

# Numpy constants
cos, sin, tan, pi = np.cos, np.sin, np.tan, np.pi


def angle_matrix(width, length, fov, bpp, direction):
    """
    This function generates a matrix of shooting angles for the image plane

    width -> int: number of pixel rows
    length -> int: number of pixel columns
    fov -> float: field of vision
    bpp -> int: light beams shot per pixel
    direction -> Vector: observer's sight direction unit vector

    return -> [[Vector]]: Matrix containing shooting angles or each pixel
    """

    # Polar angle rotation matrix
    r1x = Vector(cos(direction.z - pi / 2), 0, sin(direction.z - pi / 2), 'cart')
    r1y = Vector(0, 1, 0, 'cart')
    r1z = Vector(-sin(direction.z - pi / 2), 0, cos(direction.z - pi / 2), 'cart')
    rotation1 = Transformation(r1x, r1y, r1z)

    # Azimuthal angle rotation matrix
    r2x = Vector(cos(direction.y), sin(direction.y), 0, 'cart')
    r2y = Vector(-sin(direction.y), cos(direction.y), 0, 'cart')
    r2z = Vector(0, 0, 1, 'cart')
    rotation2 = Transformation(r2x, r2y, r2z)

    # Determine image plain vertices for an image plane center in the x axis
    ul_point = Vector(1, tan(fov), tan(fov), 'cart')
    dr_point = Vector(1, -tan(fov), -tan(fov), 'cart')
    dl_point = Vector(1, tan(fov), -tan(fov), 'cart')

    # Apply transformations to get vertices in image plane according to observer viewing direction
    ul_point = ul_point.transform(rotation1).transform(rotation2).cart2sph()
    dr_point = dr_point.transform(rotation1).transform(rotation2).cart2sph()
    dl_point = dl_point.transform(rotation1).transform(rotation2).cart2sph()

    # Image plane unit vectors
    x_unit_vector = dr_point.add(-1 * dl_point).normalize()
    y_unit_vector = ul_point.add(-1 * dl_point).normalize()

    matrix = []
    for ii in tqdm(range(width), desc="Generating angles"):
        line = []
        for jj in range(length):
            angles = []

            # Calculate pixel direction vector
            x_displacement = 2 * (jj / width) * tan(fov) * x_unit_vector
            y_displacement = 2 * (ii / length) * tan(fov) * y_unit_vector
            pixel_dir = dl_point.add(x_displacement).add(y_displacement).normalize()

            for _ in range(bpp):
                angles += [pixel_dir]

            line += [angles]
        matrix += [line]

    return matrix[::-1]
