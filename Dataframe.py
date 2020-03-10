from Algorithm import *
import matplotlib.pyplot as plt

# Constantes numpy
pi, cos, sin, sqrt, tan, arctan = np.pi, np.cos, np.sin, np.sqrt, np.tan, np.arctan2


def angle_matrix(width, length, fov, dist_to_plane, ppp, direc):
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

    matrix, angles = [], []
    pixel_side = (2 * dist_to_plane / width) * tan(fov)
    for ii in range(width):
        line = []
        for jj in range(length):
            angles = []
            # Número dos píxeis nos eixos x e y
            pp_x = jj - 0.5 * (length - 1)
            pp_y = ii - 0.5 * (width - 1)

            # Definir os ângulos que delimitam cada píxel
            po_up = arctan(-(pp_x + 0.4) * pixel_side, dist_to_plane)
            po_down = arctan(-(pp_x - 0.4) * pixel_side, dist_to_plane)
            az_up = arctan((pp_y + 0.4) * pixel_side, dist_to_plane)
            az_down = arctan((pp_y - 0.4) * pixel_side, dist_to_plane)

            for _ in range(ppp):
                # Gerar ângulos aleatórios dentro da área de cada píxel
                polar = (po_up - po_down) * np.random.random() + po_down + direc.y
                azimuth = (az_up - az_down) * np.random.random() + az_down + direc.z
                angles += [Vector(1, polar, azimuth, 'sph')]

            line += [angles]
        matrix += [line]

    return matrix


# Variáveis de posição e direção
posicao = Vector(raio_terra + 1.8, 0, 0, 'sph')
direcao, sun_direcao = Vector(1, 0, pi/2, 'sph'), Vector(1, 0, pi/2 - 0.05, 'sph')

# Parâmetros de ajuste do plano de imagem e tone mapping
DIST_TO_PLANE, FOV = 1, pi/12
WIDTH, LENGTH, PPP = 20, 20, 1
EXPOSURE, STRETCH = 18000, 2.2
BOOST = 50000

# Criar a matriz de ângulos
ANGLE_MATRIX = angle_matrix(WIDTH, LENGTH, FOV, DIST_TO_PLANE, PPP, direcao)

# Ciclo da matriz de resultados
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

# Mostrar a imagem gerada
plt.figure(dpi=200)
plt.imshow(RESULTS, cmap='hot', interpolation='bilinear', extent=[-1, 1, -1, 1])
plt.savefig('goodsky.png')
plt.show()
