from Algorithm import *
import matplotlib.pyplot as plt

# Numpy constants
pi, cos, sin, sqrt, tan, arctan = np.pi, np.cos, np.sin, np.sqrt, np.tan, np.arctan2


def angle_matrix(width, length, fov, dist_to_plane, ppp, direc):
    """
    This function creates a matrix of angles for the image grid

    width {int}: number of pixel lines
    length {int}: number of pixel columns
    fov {float}: field of view
    dist_to_plane {float}: distance from observer to image grid
    ppp {int}: light rays shot per pixel
    direc {Vector}: direction unitary vector of the observer's viewpoint

    return {[[Vector]]}: Matrix with all the firing angles for each light ray in each pixel
    """

    matrix, angles = [], []
    pixel_side = (2 * dist_to_plane / width) * tan(fov)
    for ii in range(width):
        line = []
        for jj in range(length):
            angles = []
            # Number of pixels in x and y axis'
            pp_x = jj - 0.5 * (length - 1)
            pp_y = ii - 0.5 * (width - 1)

            # Defining the boundary angles of each pixel
            po_up = arctan(-(pp_x + 0.4) * pixel_side, dist_to_plane)
            po_down = arctan(-(pp_x - 0.4) * pixel_side, dist_to_plane)
            az_up = arctan((pp_y + 0.4) * pixel_side, dist_to_plane)
            az_down = arctan((pp_y - 0.4) * pixel_side, dist_to_plane)

            for _ in range(ppp):
                # Generating random angles inside each pixel
                polar = (po_up - po_down) * np.random.random() + po_down + direc.y
                azimuth = (az_up - az_down) * np.random.random() + az_down + direc.z
                angles += [Vector(1, polar, azimuth, 'sph')]

            line += [angles]
        matrix += [line]

    return matrix


# Position and direction variables
posicao = Vector(raio_terra + 1.8, 0, 0, 'sph')
direcao, sun_direcao = Vector(1, 0, pi/2, 'sph'), Vector(1, 0, pi/2 - 0.05, 'sph')

# Image plane and tone mapping parameters
DIST_TO_PLANE, FOV = 1, pi/12
WIDTH, LENGTH, PPP = 20, 20, 1
EXPOSURE, STRETCH = 18000, 2.2
BOOST = 50000

# Matrix of angles
ANGLE_MATRIX = angle_matrix(WIDTH, LENGTH, FOV, DIST_TO_PLANE, PPP, direcao)

# Cicle to obtain results matrix, the final image
wavelengths = [(6.8e-7, 0.685), (5.5e-7, 0.81), (4.4e-7, 0.775)]
LUMINANCE, RESULTS = [], []
for i in range(WIDTH):
    linha, lum_line = [], []
    for j in range(LENGTH):
        pixel = []
        for (w, sun_emit) in wavelengths:
            intensity = 0
            for k in range(PPP):
                intensity += main_cicle(posicao, ANGLE_MATRIX[i][j][k], sun_direcao, w)
            tone_map = (1 - np.exp(-sun_emit * intensity * EXPOSURE))**(1/STRETCH)
            pixel += [tone_map]
        linha += [pixel]
    RESULTS += [linha]

# Show the final image
plt.figure(dpi=200)
plt.imshow(RESULTS, cmap='hot', interpolation='bilinear', extent=[-1, 1, -1, 1])
plt.savefig('goodsky.png')
plt.show()
