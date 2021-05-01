from itertools import product

from astropy import units as u
from jinja2 import Template
from numpy.linalg import norm
from poliastro.bodies import Jupiter as Jupiter_poliastro
from poliastro.bodies import Mars as Mars_poliastro
from poliastro.bodies import Mercury as Mercury_poliastro
from poliastro.bodies import Moon as Moon_poliastro
from poliastro.bodies import Neptune as Neptune_poliastro
from poliastro.bodies import Saturn as Saturn_poliastro
from poliastro.bodies import Sun as Sun_poliastro
from poliastro.bodies import Uranus as Uranus_poliastro
from poliastro.bodies import Venus as Venus_poliastro

# All interesting 3D directions
R_SET, V_SET = [list(product([-1, 0, 1], repeat=3)) for _ in range(2)]

# Name of the bodies
BODIES_NAMES = [
    "Sun",
    "Mercury",
    "Venus",
    "Moon",
    "Mars",
    "Jupiter",
    "Saturn",
    "Uranus",
    "Neptune",
]

POLIASTRO_BODIES = [
    Sun_poliastro,
    Mercury_poliastro,
    Venus_poliastro,
    Moon_poliastro,
    Mars_poliastro,
    Jupiter_poliastro,
    Saturn_poliastro,
    Uranus_poliastro,
    Neptune_poliastro,
]

BODY_FROM_NAME = dict(zip(BODIES_NAMES, POLIASTRO_BODIES))

if __name__ == "__main__":

    # Open the template
    with open("template_frames.script") as template_file:
        template_raw = Template(template_file.read())

        # Generate a particular configuration of body + r_vec + v_vec
        for body_name in BODIES_NAMES:
            for r_vec in R_SET:
                for v_vec in V_SET:

                    # Generate a filename hosting simulation information for the script
                    (rx, ry, rz), (vx, vy, vz) = r_vec, v_vec
                    filename = f"gmat_validate_{body_name}_frames_R{rx}{ry}{rz}_V{vx}{vy}{vz}.script"

                    # Compute Body radius in kilometers
                    r_norm, v_norm = [norm(vec) for vec in [r_vec, v_vec]]
                    R = BODY_FROM_NAME[body_name].R.to(u.km).value

                    # Correction factor to normalize position and velocity vectors
                    k_r = R / r_norm if r_norm != 0 else 1.00
                    k_v = 1 / v_norm if v_norm != 0 else 1.00

                    # Make a position vector who's norm is equal to the body's radius. Make a
                    # unitary velocity vector. Units are in [m] and [m / s].
                    rx, ry, rz = [float(k_r * r_i) for r_i in r_vec]
                    vx, vy, vz = [float(k_v * v_i) for v_i in v_vec]

                    # Fill the template
                    filled_template = template_raw.render(
                        body=body_name,
                        rx=rx,
                        ry=ry,
                        rz=rz,
                        vx=vx,
                        vy=vy,
                        vz=vz,
                    )

                    with open(f"frames_tests/{filename}", "w") as filled_template_file:
                        print(filled_template, file=filled_template_file)
