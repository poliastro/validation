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
R_SET = [
    vector for vector in product([-1, 0, 1], repeat=3) if list(vector) != [0, 0, 0]
]
V_SET = R_SET.copy()

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


def main():

    # Open the template
    with open("template_frames.script.tpl") as template_file:
        template_raw = Template(template_file.read())

    # Generate a particular configuration of body + r_vec + v_vec
    for body_name in BODIES_NAMES:
        for r_vec in R_SET:
            for v_vec in V_SET:

                # Generate a filename hosting simulation information for the script
                (rx, ry, rz), (vx, vy, vz) = r_vec, v_vec
                filename = f"gmat_validate_{body_name}_frames_R{rx}{ry}{rz}_V{vx}{vy}{vz}.script"

                # Compute vector norms and body radius in kilometers
                r_norm, v_norm = [norm(vec) for vec in [r_vec, v_vec]]
                R = BODY_FROM_NAME[body_name].R.to(u.km).value

                # Make a position vector who's norm is equal to the body's radius. Make a
                # unitary velocity vector. Units are in [km] and [km / s].
                rx, ry, rz = [r_i * R / r_norm for r_i in r_vec]
                vx, vy, vz = [v_i / v_norm for v_i in v_vec]

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
                    # print(filled_template, file=filled_template_file)
                    filled_template_file.write(filled_template)

                # Close the filled template_file
                filled_template_file.close()

    # Close the template file
    template_file.close()


if __name__ == "__main__":
    main()
