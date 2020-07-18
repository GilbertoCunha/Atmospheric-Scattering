import numpy as np
from tqdm import tqdm
from Vector import Vector, Transformation

cos, sin, tan, pi = np.cos, np.sin, np.tan, np.pi


def angle_matrix(width, length, fov, ppp, direction):
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
    r1x = Vector(cos(direction.z - pi / 2), 0, sin(direction.z - pi / 2), 'cart')
    r1y = Vector(0, 1, 0, 'cart')
    r1z = Vector(-sin(direction.z - pi / 2), 0, cos(direction.z - pi / 2), 'cart')
    rotation1 = Transformation(r1x, r1y, r1z)

    # Matriz de rotação do ângulo azimutal
    r2x = Vector(cos(direction.y), sin(direction.y), 0, 'cart')
    r2y = Vector(-sin(direction.y), cos(direction.y), 0, 'cart')
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
