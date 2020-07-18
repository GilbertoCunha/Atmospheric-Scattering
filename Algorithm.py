import numpy as np
from Vector import Vector

# Numpy constants
pi, cos, sin, sqt, tan = np.pi, np.cos, np.sin, np.sqrt, np.tan

# Algorithmic constants
g, steps, prob_scatter = 0.996, 16, 0.5

"""Global Constants - you can play with these to pretend you are on different planets
beta_mie and beta_rayleigh are a pre-calculation of the extinction coefficients
without accounting for the wavelength"""
n = 1.00029
earth_radius, atm_height = 6.3781e6, 8e4
beta_mie, beta_rayleigh = 5.76e-7, (8 / (3 * 2.504e25)) * (pi ** 3) * (n ** 2 - 1) ** 2
H_M, H_R = 1200, 8500

# Construction of atmosphere spherical shells
shell_height = np.geomspace(10, atm_height - 800, steps)
log_shell_height = [atm_height - height for height in shell_height]
log_shell_height = log_shell_height[::-1]
shell_height += [atm_height]


def dist_min_shell(r, d, num, shell):
    """
    This function returns the distance between the traveling light and one of the atmospheric
    shells according to the light direction

    r -> Vector: position vector in spherical coordinates
    dir -> Vector: direction unit vector in spherical coordinates
    num -> int: shell_height index, atmospheric shell to calculate distance with respect to
    shell -> String: the type of shell (two types were created), 'log' or 'exp' to calculate
    distance to

    return -> float: distance to atmospheric shell
    """

    angle = pi - r.angle(d)

    r2 = 0
    if shell == 'exp':
        r2 = earth_radius + shell_height[num]
    elif shell == 'log':
        r2 = earth_radius + log_shell_height[num]

    discriminant = r.x ** 2 * cos(angle) ** 2 + r2 ** 2 - r.x ** 2
    dist_atm = r.x * cos(angle) + sqt(discriminant)
    return dist_atm


def scatter_transmission(h, angle, w):
    """
    This function calculates the transmission coefficients of light in a scattering

    h -> float: scattering height
    angle -> float: scattering angle
    w -> float: wavelength of light

    return -> [Float]:
        return[0]: Rayleigh transmission coefficient
        return[1]: Mie transmission coefficient
    """

    # Mie and Rayleigh density
    rho_mie = np.exp(-h / H_M)
    rho_rayleigh = np.exp(-h / H_R)

    # Mie and Rayleigh phase
    phase_rayleigh = (3 / 16*pi) * (1 + (cos(angle)) ** 2)
    phase_mie = 1.5 * ((1 - g ** 2) / (2 + g ** 2)) * (1 + (cos(angle)) ** 2) / (
                (1 + g ** 2 - 2 * g * cos(angle)) ** 1.5)

    # Joining the functions to calculate the coefficients
    beta_ray = beta_rayleigh / (w ** 4)
    ratio_rayleigh = beta_ray * rho_rayleigh * phase_rayleigh
    ratio_mie = beta_mie * rho_mie * phase_mie

    return ratio_rayleigh, ratio_mie


def sum_to_opt_dep(r, dis, w, one_step=False):
    """
    This function calculates the optical depth of light in a certain trajectory, updating
    the value at each intersection with the exponential atmospheric shells

    r -> Vector: position vector in spherical coordinates
    dis -> Vector: light's displacement vector in spherical coordinates
    w -> float: wavelength of light

    return -> float: optical depth in the path given by the position and displacement vectors
    """

    beta_ray = beta_rayleigh / (w ** 4)

    # If we update only by one step
    if one_step:
        # Calculating density
        rho_mie = np.exp(- (r.x - earth_radius) / H_M)
        rho_rayleigh = np.exp(- (r.x - earth_radius) / H_R)

        # Optical depth in this step
        opt_dep = 4 * pi * (beta_mie * rho_mie + beta_ray * rho_rayleigh) * dis.x
        return opt_dep

    # Maximum and minimum height in the trajectory
    r_f = r.add(dis)
    max_height = max(r_f.x, r.x) - earth_radius
    min_height = r.x * sin(pi - r.angle(dis)) - earth_radius

    # Determine spherical shells to go through
    n_min, n_max = 0, steps - 1
    for step in range(steps):
        if shell_height[step] > max_height - earth_radius:
            n_max = step - 1
            break
        elif shell_height[step] < min_height - earth_radius:
            n_min = step + 1

    # Main calculation cycle
    opt_dep, num = 0, n_min
    while num <= n_max:
        # next shell to go through
        angle = pi - r.angle(dis)
        min_angle = np.arcsin(r.x / shell_height[num - 1])
        if angle < min_angle:
            num -= 1

        # Calculate displacement to that shell and update position vector
        step_displace = Vector(dist_min_shell(r, dis.normalize(), num, 'exp'), dis.y, dis.z, 'sph')
        r = r.add(step_displace)

        # Calculate density
        rho_mie = np.exp(- (r.x - earth_radius) / H_M)
        rho_rayleigh = np.exp(- (r.x - earth_radius) / H_R)

        # calculate optical depth increment
        opt_dep += 4 * pi * (beta_mie * rho_mie + beta_ray * rho_rayleigh) * step_displace.x

        if angle < min_angle:
            # Select the previous shell to intersect the same shell again
            # after increment num
            num -= 1

    return opt_dep


def scatter_per_step(r, direction, sun_dir, w, opt_t):
    """
    This function calculates the transmission coefficients of scattered light at a certain
    position taking into account optical depth and a ray-tracing method of sending light
    directly to the sun for each scattering event, so as to increase the amount of light
    that reaches the observer, brightening the whole final image

    r -> Vector: light's position vector in spherical coordinates
    direction -> Vector: light's direction unit vector in spherical coordinates
    sun_dir -> Vector: sun's direction unit vector in spherical coordinates
    w -> float: light's wavelength
    opt_t -> float: optical depth of the light before reaching this position

    return -> [[Float]]:
        return[0]: Transmission coefficients of scattered light:
            return[0][0]: Rayleigh transmission coefficient
            return[0][1]: Mie transmission coefficient
        return[1]: Transmission coefficients of light sent to the sun:
            return[1][0]: Rayleigh transmission coefficient
            return[1][1]: Mie transmission coefficient
    """

    # Calculating scattering transmission of each light beam
    scatter_direction = Vector(1, 2 * pi * np.random.random(), pi * np.random.random(), 'sph')
    scatter_angle = direction.angle(scatter_direction)
    direct_angle = direction.angle(sun_dir)
    scatter_tr, scatter_tm = scatter_transmission(r.x - earth_radius, scatter_angle, w)
    direct_tr, direct_tm = scatter_transmission(r.x - earth_radius, direct_angle, w)

    # Increment optical depth in this path
    distance = min(min(r.x - earth_radius, dist_min_shell(r, scatter_direction, -1, 'exp')), 0.5 * atm_height)
    distance *= np.random.random()
    scatter_displacement = Vector(distance, scatter_direction.y, scatter_direction.z, 'sph')
    direct_displacement = Vector(dist_min_shell(r, sun_dir, -1, 'exp'), sun_dir.y, sun_dir.z, 'sph')
    scatter_opt_t = opt_t + sum_to_opt_dep(r, scatter_displacement, w)
    direct_opt_t = opt_t + sum_to_opt_dep(r, direct_displacement, w)

    # Rayleigh and Mie transmission coefficients for direct light beam
    direct_tr *= np.exp(-direct_opt_t)
    direct_tm *= np.exp(-direct_opt_t)

    # Taking into account multiple scattering
    prob = np.random.random()
    r = r.add(scatter_displacement)
    if prob < prob_scatter:
        # New transmission coefficients
        scatter_vector = scatter_per_step(r, scatter_direction, sun_dir, w, scatter_opt_t)

        # Summing the contributions of the new direct light beams
        direct_tr += direct_tr * scatter_vector[1][0]
        direct_tm += direct_tm * scatter_vector[1][1]

        # Updating the transmission coefficient of the scattered light
        scatter_tr *= scatter_vector[0][0]
        scatter_tm *= scatter_vector[0][1]
        return [scatter_tr, scatter_tm], [direct_tr, direct_tm]

    # In case light doesn't scatter, send it to the sun
    else:
        # Update the transmission of the scattered light
        final_angle = scatter_direction.angle(sun_dir)
        final_t = scatter_transmission(r.x - earth_radius, final_angle, w)
        scatter_tr *= final_t[0]
        scatter_tm *= final_t[1]

        # Update optical depth in new path and calculate final transmission coefficients
        scatter_displacement = Vector(dist_min_shell(r, sun_dir, -1, 'exp'), sun_dir.y, sun_dir.z, 'sph')
        scatter_opt_t += sum_to_opt_dep(r, scatter_displacement, w)
        scatter_vector = [scatter_tr * np.exp(-scatter_opt_t), scatter_tm * np.exp(-scatter_opt_t)]
        return scatter_vector, [direct_tr, direct_tm]


def main_cycle(r, direction, sun_dir, w):
    """
    This function calculates the final transmission coefficient of light with a specific
    wavelength starting in the observer and scattering throughout the whole atmosphere until
    reaching the sun.
    Scattering events occur in the intersection of the light beam with each of the logarithmic
    shells defined above.

    r -> Vector: light's position vector in spherical coordinates
    direction -> Vector: light's direction unit vector in spherical coordinates
    sun_dir -> Vector: sun's direction unit vector in spherical coordinates
    w -> float: light's wavelength

    return -> float: Total transmission coefficient of light that reaches the observer
    """

    # If the light beam goes below the earth
    min_angle = np.arcsin(earth_radius / r.x)
    if (pi - r.angle(direction)) < min_angle:
        return 0

    # Determine all intersected spherical logarithmic shells
    n_min = 0
    for num in range(steps):
        if log_shell_height[num] > (r.x - earth_radius):
            break
        else:
            n_min = num + 1
    nums = [num for num in range(n_min, steps - 1)]

    # Go through each of the shells
    opt_dep, t = 0, 0
    for num in nums:
        # Calculate displacement and update position vector
        displacement = Vector(dist_min_shell(r, direction, num, 'log'), direction.y, direction.z, 'sph')
        r = r.add(displacement)

        # Calculate optical depth in this displacement
        opt_dep += sum_to_opt_dep(r, displacement, w, one_step=True)

        # Scatter the light at this position
        ts = scatter_per_step(r, direction, sun_dir, w, opt_dep)
        t += sum(ts[0]) + sum(ts[1])

    # Go through the last shell
    displacement = Vector(dist_min_shell(r, direction, -1, 'log'), direction.y, direction.z, 'sph')
    r = r.add(displacement)
    opt_dep += sum_to_opt_dep(r, displacement, w, one_step=True)

    # Send the light beam to the sun
    scatter_angle = direction.angle(sun_dir)
    scatter_tr, scatter_tm = scatter_transmission(r.x, scatter_angle, w)
    t += scatter_tr * np.exp(-opt_dep) + scatter_tm * np.exp(-opt_dep)

    return t
