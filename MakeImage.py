# Import my libraries
from Algorithm import *
from Vector import *

# Import external libraries
import matplotlib.pyplot as plt
from itertools import product
import concurrent.futures
from tqdm import tqdm
import time

# Constantes numpy
pi, cos, sin, sqrt, tan, arctan = np.pi, np.cos, np.sin, np.sqrt, np.tan, np.arctan2

# Variáveis de posição e direção
posicao = Vector(earth_radius + 8000, 0, 0, 'sph')
direcao, sun_direcao = Vector(1, 0, pi/2, 'sph'), Vector(1, 0, pi/2 - 0.05, 'sph')
# direcao, sun_direcao = Vector(1, 0, pi/4, 'sph'), Vector(1, 0, pi/4, 'sph')

# Parâmetros de ajuste da imagem
wavelengths = [(6.8e-7, 0.685), (5.35e-7, 0.81), (4.6e-7, 0.775)]
DIST_TO_PLANE, FOV = 1, pi/12
WIDTH, LENGTH, PPP = 100, 100, 1
EXPOSURE, STRETCH = 1.8e4, 2.2


def angle_matrix(width, length, fov, ppp, direc):
    """
    Esta função cria a matriz de ângulos para o plano de imagem

    width {int}: número de linhas de píxeis
    length {int}: número de colunas de píxeis
    fov {float}: campo de visão
    dist_to_plane {float}: distância do observador ao plano de imagem
    ppp {int}: fotões a disparar por píxel
    direc {Vector}: vetor direção (unitário) do olhar do observador

    return {[[Vector]]}: Matriz com todos os ângulos de disparo para cada fotão
    """

    # Campos de visão em ambas as direções
    fovx = fov
    fovy = fov * (width / length)

    # Matriz de rotação do ângulo polar
    r1x = Vector(cos(direc.z-pi/2), 0, sin(direc.z-pi/2), 'cart')
    r1y = Vector(0, 1, 0, 'cart')
    r1z = Vector(-sin(direc.z-pi/2), 0, cos(direc.z-pi/2), 'cart')
    rotation1 = Transformation(r1x, r1y, r1z)

    # Matriz de rotação do ângulo azimutal
    r2x = Vector(cos(direc.y), sin(direc.y), 0, 'cart')
    r2y = Vector(-sin(direc.y), cos(direc.y), 0, 'cart')
    r2z = Vector(0, 0, 1, 'cart')
    rotation2 = Transformation(r2x, r2y, r2z)

    # Determinar vértices do plano imagem caso este esteja no eixo x
    ul_point = Vector(1, tan(fovx), tan(fovy), 'cart')
    dr_point = Vector(1, -tan(fovx), -tan(fovy), 'cart')
    dl_point = Vector(1, tan(fovx), -tan(fovy), 'cart')

    # Aplicar transformações aos pontos
    ul_point = ul_point.transform(rotation1).transform(rotation2).cart2sph()
    dr_point = dr_point.transform(rotation1).transform(rotation2).cart2sph()
    dl_point = dl_point.transform(rotation1).transform(rotation2).cart2sph()

    # Versores do plano imagem
    x_direc = dr_point.add(-1 * dl_point).normalize()
    y_direc = ul_point.add(-1 * dl_point).normalize()

    matrix = []
    for ii in tqdm(range(width), desc="Generating angles"):
        line = []
        for jj in range(length):
            angles = []

            # Calcular vetor direção do píxel
            x_disp = 2 * (jj / width) * tan(fov) * x_direc
            y_disp = 2 * (ii / length) * tan(fov) * y_direc
            pixel_dir = dl_point.add(x_disp).add(y_disp).normalize()

            for _ in range(ppp):
                angles += [pixel_dir]

            line += [angles]
        matrix += [line]

    return matrix[::-1]


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
            intensity += main_cycle(posicao, ANGLE_MATRIX[width][length][k], sun_direcao, w)
        tone_map = (1 - np.exp(-sun_emit * intensity * EXPOSURE)) ** (1 / STRETCH)
        pixel += [tone_map]

    return pixel


# Criar a matriz de ângulos e resultados
start = time.time()
ANGLE_MATRIX = angle_matrix(WIDTH, LENGTH, FOV, PPP, direcao)
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
