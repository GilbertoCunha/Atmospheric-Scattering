# Import internal libraries
from Algorithm import *
from Vector import Vector
from Aux import angle_matrix

# Import external libraries
import matplotlib.pyplot as plt
from itertools import product
import concurrent.futures
from tqdm import tqdm
import time

# Position and direction variables
position = Vector(earth_radius + 1.8, 0, 0, 'sph')
direction, sun_direction = Vector(1, 0, pi / 2, 'sph'), Vector(1, 0, pi / 2 - 0.05, 'sph')

# Parâmetros de ajuste da imagem
wavelengths = [(6.8e-7, 0.685), (5.35e-7, 0.81), (4.6e-7, 0.775)]
DIST_TO_PLANE, FOV = 1, pi/12
WIDTH, LENGTH, PPP = 100, 100, 1
EXPOSURE, STRETCH = 1.8e4, 2.2


def do_pixel(arg):
    """
    Esta função calcula uma linha da imagem total

    args {tuple}:
        mat_pos: coluna e linha do píxel na matriz
        wavs: comprimentos de onda a considerar
        ppp: número de raios de luz disparados por píxel
        pos: posição do observador
        sun_dir: direção do sol

    return {[(float, float, float)]}: Cálculo de cor para a correspondente linha da imagem
    """

    width, length = int(arg[0]), int(arg[1])
    pixel = []
    for (w, sun_emit) in wavelengths:
        intensity = 0
        for k in range(PPP):
            intensity += main_cycle(position, ANGLE_MATRIX[width][length][k], sun_direction, w)
        tone_map = (1 - np.exp(-sun_emit * intensity * EXPOSURE)) ** (1 / STRETCH)
        pixel += [tone_map]

    return pixel


# Criar a matriz de ângulos e resultados
start = time.time()
ANGLE_MATRIX = angle_matrix(WIDTH, LENGTH, FOV, PPP, direction)
RESULTS = []

# Paralelizar o cálculo de cada píxel
if __name__ == '__main__':
    with concurrent.futures.ProcessPoolExecutor() as executor:
        wid_list, len_list = np.linspace(0, WIDTH-1, WIDTH), np.linspace(0, LENGTH-1, LENGTH)

        args = tqdm(product(wid_list, len_list), total=WIDTH*LENGTH, desc="Pre-processing")

        results = list(tqdm(executor.map(do_pixel, args), total=WIDTH*LENGTH, desc="Calculating"))

        for result in results:
            RESULTS.append(result)


RESULTS = np.asarray(RESULTS).reshape((WIDTH, LENGTH, 3))

stop = time.time()
print(f'This took {stop-start:.2f} seconds | {(stop-start)/60:.2f} minutes | {(stop-start)/3600:.2f} hours')

# Mostrar a imagem gerada
plt.figure(dpi=800)
plt.imshow(RESULTS, cmap='hot', interpolation='bilinear', extent=[-1, 1, -1, 1])
plt.savefig('sky1T.png')
plt.show()
