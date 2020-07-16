import numpy as np
from Vector import Vector

# Constantes numpy
pi, cos, sin, sqrt = np.pi, np.cos, np.sin, np.sqrt

# Constantes de ajuste
g, steps, prob_scatter = 0.996, 16, 0.5

# Constantes globais
n = 1.00029
# T = 20  # Temperatura ao nível do mar em graus Celsius
# P = 101  # Pressão ao nível do mar em kPa
raio_terra, altura_atm = 6.3781e6, 8e4
beta_mie, beta_rayleigh = 5.76e-7, (8 / (3 * 2.504e25)) * (pi ** 3) * (n ** 2 - 1) ** 2
# beta_mie, beta_rayleigh = 0, (8 / (3 * 2.504e25)) * (pi ** 3) * (n ** 2 - 1) ** 2
H_M, H_R = 1200, 8500

# Construção das camadas esféricas da atmosfera
shell_height = np.geomspace(10, altura_atm - 800, steps)
log_shell_height = [altura_atm - altura for altura in shell_height]
log_shell_height = log_shell_height[::-1]
shell_height += [altura_atm]


# log_shell_height = np.linspace(10, altura_atm - 800, steps)


def dist_min_log_shell(r, d, num):
    """
    Esta função retorna a distância entre um fotão e uma das camadas esféricas da atmosfera
    de acordo com a sua direção

    r {Vector}: vetor posição em coordenadas esféricas
    dir {Vector}: versor direção em coordenadas esféricas
    num {int}: índice de shell_height, camada da atmosfera à qual se deve calcular a distância

    return: distância à camada esférica
    """
    angle = pi - r.angle(d)
    r2 = raio_terra + log_shell_height[num]
    discriminant = r.x ** 2 * cos(angle) ** 2 + r2 ** 2 - r.x ** 2
    dist_atm = r.x * cos(angle) + sqrt(discriminant)
    return dist_atm


def dist_min_shell(r, d, num):
    """
    Esta função retorna a distância entre um fotão e uma das camadas esféricas da atmosfera
    de acordo com a sua direção

    r {Vector}: vetor posição em coordenadas esféricas
    dir {Vector}: versor direção em coordenadas esféricas
    num {int}: índice de shell_height, camada da atmosfera à qual se deve calcular a distância

    return: distância à camada esférica
    """
    angle = pi - r.angle(d)
    r2 = raio_terra + shell_height[num]
    discriminant = r.x ** 2 * cos(angle) ** 2 + r2 ** 2 - r.x ** 2
    dist_atm = r.x * cos(angle) + sqrt(discriminant)
    return dist_atm


def scatter_transmission(h, angle, w):
    """
    Esta função calcula os índices de transmissão de um fotão num espalhamento

    h {float}: altura do ponto
    angle {float}: ângulo de espalhamento
    w {float}: comprimento de onda

    return {[Float]}:
    return[0]: índice de transmissão de Rayleigh
    return[1]: índice de transmissão de Mie
    """

    # Cálculo das densidades de mie e de rayleigh
    rho_mie = np.exp(-h / H_M)
    rho_rayleigh = np.exp(-h / H_R)

    # Cálculo da fase de rayleigh e de mie
    phase_rayleigh = (3 / 16*pi) * (1 + (cos(angle)) ** 2)
    phase_mie = 1.5 * ((1 - g ** 2) / (2 + g ** 2)) * (1 + (cos(angle)) ** 2) / (
                (1 + g ** 2 - 2 * g * cos(angle)) ** 1.5)

    # Soma das contribuições dos espalhamentos de mie e de rayleigh
    beta_ray = beta_rayleigh / (w ** 4)
    ratio_rayleigh = beta_ray * rho_rayleigh * phase_rayleigh
    ratio_mie = beta_mie * rho_mie * phase_mie

    return ratio_rayleigh, ratio_mie


def sum_to_opt_dep(r, dis, w, one_step=False):
    """
    Esta função calcula a optical depth de um fotão num determinado trajeto através de interseção de camadas esféricas

    r {Vector}: vetor posição em coordenadas esféricas
    dis {Vector}: vetor deslocamento do fotão em coordenadas esféricas
    w {float}: comprimento de onda do fotão

    return {float}: optical depth no percurso descrito pelos vetores posição e deslocamento
    """

    beta_ray = beta_rayleigh / (w ** 4)
    # Caso estejamos a atualizar por um único step
    if one_step:
        # Cálculo das densidades para cada passo
        rho_mie = np.exp(- (r.x - raio_terra) / H_M)
        rho_rayleigh = np.exp(- (r.x - raio_terra) / H_R)

        # Optical depth in this step
        opt_dep = 4 * pi * (beta_mie * rho_mie + beta_ray * rho_rayleigh) * dis.x
        return opt_dep

    # Calcular altura mínima e máxima
    r_f = r.add(dis)
    max_height = max(r_f.x, r.x) - raio_terra
    min_height = r.x * sin(pi - r.angle(dis)) - raio_terra

    # Determinar as camadas esféricas a percorrer:
    n_min, n_max = 0, steps - 1
    for step in range(steps):
        if shell_height[step] > max_height - raio_terra:
            n_max = step - 1
            break
        elif shell_height[step] < min_height - raio_terra:
            n_min = step + 1

    # Ciclo principal do cálculo da optical depth
    opt_dep, num = 0, n_min
    while num <= n_max:
        # Determinar a camada esférica a perfurar
        angle = pi - r.angle(dis)
        min_angle = np.arcsin(r.x / shell_height[num - 1])
        if angle < min_angle:
            num -= 1

        # Calcular o vetor deslocamento e atualizar o vetor posição
        step_displace = Vector(dist_min_shell(r, dis.normalize(), num), dis.y, dis.z, 'sph')
        r = r.add(step_displace)

        # Cálculo das densidades para cada passo
        rho_mie = np.exp(- (r.x - raio_terra) / H_M)
        rho_rayleigh = np.exp(- (r.x - raio_terra) / H_R)

        # Incremento da opt_depth
        opt_dep += 4 * pi * (beta_mie * rho_mie + beta_ray * rho_rayleigh) * step_displace.x
        # opt_dep += 4 * pi * beta_ray * rho_rayleigh * step_displace.x

        # Diminuir a camada para intersetar novamente a camada esférica inferior
        if angle < min_angle:
            num -= 1

    return opt_dep


def scatter_per_step(r, direction, sun_dir, w, opt_t):
    """
    Esta função é responsável por calcular as transmissões dos fotões espalhados em cada
    posição onde um fotão seja espalhado. Utiliza também um método comum no raytracing de
    enviar um fotão direto para a fonte por cada fotão espalhado, de modo a aumentar a
    intensidade de luz que chega ao olho

    r {Vector}: Vector posição do fotão em coordenadas esféricas
    direction {Vector}: versor direção do fotão em coordenadas esféricas
    sun_dir {Vector}: versor direção do sol em coordenadas esféricas
    w {float}: comprimento de onda do fotão

    return [[Float]]:
        return[0]: Coeficientes de transmissão dos fotões espalhados:
            return[0][0]: Coeficiente de transmissão de Rayleigh
            return[0][1]: Coeficiente de transmissão de Mie
        return[1]: Coeficientes de transmissão dos fotões diretos:
            return[1][0]: Coeficiente de transmissão de Rayleigh
            return[1][1]: Coeficiente de transmissão de Mie
    """

    # Cálculo das transmissões de cada fotão
    scatter_direction = Vector(1, 2 * pi * np.random.random(), pi * np.random.random(), 'sph')
    scatter_angle = direction.angle(scatter_direction)
    direct_angle = direction.angle(sun_dir)
    scatter_tr, scatter_tm = scatter_transmission(r.x - raio_terra, scatter_angle, w)
    direct_tr, direct_tm = scatter_transmission(r.x - raio_terra, direct_angle, w)

    # Incrementar a optical depth neste percurso
    distance = min(min(r.x - raio_terra, dist_min_shell(r, scatter_direction, -1)), 0.5 * altura_atm)
    distance *= np.random.random()
    scatter_displacement = Vector(distance, scatter_direction.y, scatter_direction.z, 'sph')
    direct_displacement = Vector(dist_min_shell(r, sun_dir, -1), sun_dir.y, sun_dir.z, 'sph')
    scatter_opt_t = opt_t + sum_to_opt_dep(r, scatter_displacement, w)
    direct_opt_t = opt_t + sum_to_opt_dep(r, direct_displacement, w)

    # Cálculo das transmitâncias de rayleigh e mie para o raio direto ao sol
    direct_tr *= np.exp(-direct_opt_t)
    direct_tm *= np.exp(-direct_opt_t)

    # Caso o fotão espalhado se espalhe novamente
    prob = np.random.random()
    r = r.add(scatter_displacement)
    if prob < prob_scatter:
        # Novos coeficientes de transmissão
        scatter_vector = scatter_per_step(r, scatter_direction, sun_dir, w, scatter_opt_t)

        # Somar as contribuições dos outros fotões diretos ao sol
        direct_tr += direct_tr * scatter_vector[1][0]
        direct_tm += direct_tm * scatter_vector[1][1]

        # Atualizar o valor da transmitância do fotão espalhado
        scatter_tr *= scatter_vector[0][0]
        scatter_tm *= scatter_vector[0][1]
        return [scatter_tr, scatter_tm], [direct_tr, direct_tm]

    # Caso o fotão não se espalhe, enviá-lo para o sol
    else:
        # Atualizar a transmitância do fotão espalhado agora para o sol
        final_angle = scatter_direction.angle(sun_dir)
        final_t = scatter_transmission(r.x - raio_terra, final_angle, w)
        scatter_tr *= final_t[0]
        scatter_tm *= final_t[1]

        # Atualizar a optical depth devido ao novo percurso e calcular as transmitâncias finais
        scatter_displacement = Vector(dist_min_shell(r, sun_dir, -1), sun_dir.y, sun_dir.z, 'sph')
        scatter_opt_t += sum_to_opt_dep(r, scatter_displacement, w)
        scatter_vector = [scatter_tr * np.exp(-scatter_opt_t), scatter_tm * np.exp(-scatter_opt_t)]
        return scatter_vector, [direct_tr, direct_tm]


def main_cicle(r, direction, sun_dir, w):
    """
    Esta função é o ciclo principal que calcula a transmitância total dos fotões ao longo da linha de visão até
    perfurar a atmosfera, tendo em conta eventos de scattering em cada superfície esférica divisória da atmosfera

    r {Vector}: vetor posição do fotão em coordenadas esféricas
    direction {Vector}: versor direção do fotão em coordenadas esféricas
    sun_dir {Vector}: versor direção do sol em coordenadas esféricas
    w {float}: comprimento de onda do fotão

    return {float}: Transmitância total dos fotões que atingem o olho
    """

    # Se o raio for para baixo da terra
    min_angle = np.arcsin(raio_terra / r.x)
    if (pi - r.angle(direction)) < min_angle:
        return 0

    # Calcular as superfícies esféricas a serem atravessadas
    n_min = 0
    for num in range(steps):
        if log_shell_height[num] > (r.x - raio_terra):
            break
        else:
            n_min = num + 1
    nums = [num for num in range(n_min, steps - 1)]

    # Percorrer cada uma das superfícies esféricas logaritmicas
    opt_dep, t = 0, 0
    for num in nums:
        # Calcular o vetor deslocamento e atualizar o vetor posição
        displacement = Vector(dist_min_log_shell(r, direction, num), direction.y, direction.z, 'sph')
        r = r.add(displacement)

        # Calcular a optical depth neste deslocamento
        opt_dep += sum_to_opt_dep(r, displacement, w, one_step=True)

        # Efetuar o espalhamento do raio
        ts = scatter_per_step(r, direction, sun_dir, w, opt_dep)
        t += sum(ts[0]) + sum(ts[1])

    # Percorrer a última camada da atmosfera
    displacement = Vector(dist_min_log_shell(r, direction, -1), direction.y, direction.z, 'sph')
    r = r.add(displacement)
    opt_dep += sum_to_opt_dep(r, displacement, w, one_step=True)

    # Espalhamento final para o sol
    scatter_angle = direction.angle(sun_dir)
    scatter_tr, scatter_tm = scatter_transmission(r.x, scatter_angle, w)
    t += scatter_tr * np.exp(-opt_dep) + scatter_tm * np.exp(-opt_dep)

    return t
