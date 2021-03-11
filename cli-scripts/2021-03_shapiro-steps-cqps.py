import logging
from typing import List

import numpy as np
import matplotlib.pyplot as plt

plt.style.use("../support-files/qubit.mplstyle")
from matplotlib import cm
from IPython.display import display
import ipywidgets
from scipy import integrate
import scipy.special as spl
from pyprind import ProgBar
import plotly.graph_objects as go
import plotly.figure_factory as ff
from scipy.spatial import Delaunay

eCharge = 1.6 * 10 ** (-19)
TWO_PI = 2 * np.pi
Q_0 = 2 * eCharge
uV = 10 ** (-6)
GHz = 10 ** 9
ns = 10 ** (-9)
uH = 10 ** (-6)
fF = 10 ** (-15)

###############################################################################
#                                    Kernel                                   #
###############################################################################
STATE_INITIAL = [0, 0]
CURRENT_IDX = 0


def slimmed_cqps_kernel(
    t: float,
    state_vector: List[float],
    Vb_over_L: float,
    Vmw_over_L: float,
    omega_mw: float,
    Vs_over_L: float,
    R_over_L,
) -> List[float]:
    """Kernel evalautes the derivitives for each components of the state_vector to make the next step:
            x_99 = x_98 + dx_98/dt*Δt
    - Evaluate the various fractions, Vb/L in order to speed up performace.
    - state_vector must have the order: (i, ncqps)
    """

    (i, ncqps) = state_vector

    di_dt = (
        Vb_over_L
        + (Vmw_over_L * np.cos(omega_mw * t))
        - R_over_L * i
        - Vs_over_L * np.sin(TWO_PI * ncqps)
    )
    dncqps_dt = i / Q_0

    return [di_dt, dncqps_dt]


Vb = 12 * uV
Vmw = 10 * uV
fmw = 1 * GHz
Vs = 6 * uV
R = 20_000
L = 10 * uH

t_end = 8 * ns
t_points = 201

# Initial simulation and plot
# t_list = np.linspace(0, t_end, t_points)
# simulation = integrate.odeint(
#     func=slimmed_cqps_kernel,
#     y0=STATE_INITIAL,
#     t=t_list,
#     args=(Vb / L, Vmw / L, fmw * TWO_PI, Vs / L, R / L),
#     tfirst=True,
# )

# fig, ax = plt.subplots(1, 1, figsize=(4, 3))
# (simulation_graph,) = ax.plot(t_list / ns, simulation[:, CURRENT_IDX] / Q_0 / fmw)
# ax.set_xlabel("Time, t (ns)", fontsize=12)
# ax.set_ylabel("Current, $I/2ef_{mw}$", fontsize=12)
# plt.tight_layout()
# plt.show()

###############################################################################
#                                   I-V plot                                  #
###############################################################################
def average_in_range(
    x_input: List[float], y_input: List[float], x_start: float, x_end: float
) -> float:
    """
    Takes average of y_input ([1,2,3,4,5....]) by specifying:
    - x_start
    - x_end
    use values from y_input that correspond to the selected x_input scope
    """

    x_start_idx = np.searchsorted(x_input, x_start)
    x_end_idx = np.searchsorted(x_input, x_end)
    return np.mean(y_input[x_start_idx:x_end_idx])


def average_current_slimmed_cqps_kernel(
    t_list: List[float],
    t_start: float,
    t_end: float,
    Vb: float,
    Vmw: float,
    fmw: float,
    Vs: float,
    R: float,
    L: float,
):
    """
    Kernel evaluates the average current in the system for times t_list
    """

    normalised_current = (
        integrate.odeint(
            func=slimmed_cqps_kernel,
            y0=STATE_INITIAL,
            t=t_list,
            args=(Vb / L, Vmw / L, fmw * TWO_PI, Vs / L, R / L),
            tfirst=True,
        )[:, CURRENT_IDX]
        / Q_0
        / fmw
    )

    return average_in_range(t_list, normalised_current, t_start, t_end)


# print("Generating IV-curve")
# Vmw = 10 * uV
# fmw = 1 * GHz
# Vs = 6 * uV
# R = 20_000
# L = 10 * uH

# t_start = 0
# t_end = 8 * ns
# t_points = 201

# vb_start = 0
# vb_end = 20 * uV
# vb_points = 10

# # Deriving parameters and running simulation
# # Simulations must always start from t=0, otherwise you will simply be shifting the response when changing t_start
# t_list = np.linspace(0, t_end, t_points)
# vb_list = np.linspace(vb_start, vb_end, vb_points)
# i_list = [
#     average_current_slimmed_cqps_kernel(t_list, t_start, t_end, vb, Vmw, fmw, Vs, R, L)
#     for vb in vb_list
# ]

# fig, ax = plt.subplots(1, 1, figsize=(4, 3))
# (simulation_graph,) = ax.plot(vb_list / uV, i_list)
# ax.set_xlabel("$V_{bias}$ ($\mu{V}$)", fontsize=12)
# ax.set_ylabel("Average Current, $I/2ef_{mw}$", fontsize=12)
# ax.set_title(
#     f"""
#     $V_{{mw}}={Vmw/uV:.3f}\mu{{V}}$, $f_{{mw}}={fmw/GHz:.3f}GHz$,
#     $V_s={Vs/uV:.2f}\mu{{V}}$
#     $R={R}\Omega$, $L={L/uH:.3f}\mu{{H}}$
#     """,
#     fontsize=12,
# )
# plt.tight_layout()
# plt.show()


###############################################################################
#                                   3D-Plot                                   #
###############################################################################
fmw = 1 * GHz
Vs = 5 * uV
R = 10_000
L = 2 * uH

t_start = 0
t_end = 8 * ns
t_points = 201

vb_start = 0
vb_end = 20 * uV
vb_points = 41

vmw_start = 0
vmw_end = 100 * uV
vmw_points = 101

# Deriving parameters and running simulation
# Simulations must always start from t=0, otherwise you will simply be shifting the response when changing t_start
t_list = np.linspace(0, t_end, t_points)
vb_list = np.linspace(vb_start, vb_end, vb_points)
vmw_list = np.linspace(vmw_start, vmw_end, vmw_points)
i_average_mesh = np.empty((vb_points, vmw_points))

progress_bar = ProgBar(vb_points * vmw_points, bar_char="█")
for Vb_idx, Vb in enumerate(vb_list):
    for Vmw_idx, Vmw in enumerate(vmw_list):
        i_average_mesh[Vb_idx][Vmw_idx] = average_current_slimmed_cqps_kernel(
            t_list, t_start, t_end, Vb, Vmw, fmw, Vs, R, L
        )
        progress_bar.update()

# Arrays for differnet plotting
vb_mesh, vmw_mesh = np.meshgrid(vb_list, vmw_list)
vb_mesh_flat = vb_mesh.flatten()
vmw_mesh_flat = vmw_mesh.flatten()
i_average_mesh_flat = i_average_mesh.T.flatten()

# Triangular plot #############################################################
# simplices = Delaunay(np.vstack([vb_mesh_flat, vmw_mesh_flat]).T).simplices

# fig = ff.create_trisurf(
#     x=vb_mesh_flat / uV,
#     y=vmw_mesh_flat / uV,
#     z=i_average_mesh_flat,
#     simplices=simplices,
#     width=1000,
#     height=1000,
# )
# fig.update_layout(
#     title_text=f"$f_{{mw}}={fmw/GHz:.3f}GHz, V_s={Vs/uV:.2f}\mu{{V}}, R={R}\Omega, L={L/uH:.3f}\mu{{H}}$",
#     title_font_size=40,
#     width=1000,
#     height=1000,
#     margin_pad=8,
#     scene=dict(
#         camera={
#             "eye": {"x": -1, "y": -1.25, "z": 1.25},
#         },
#         xaxis={
#             "tickfont": {"size": 15},
#             "title": {"text": "Vb (μV)", "font": {"size": 40}},
#         },
#         yaxis={
#             "tickfont": {"size": 15},
#             "title": {"text": "Vmw (μV)", "font": {"size": 40}},
#         },
#         zaxis={
#             "tickfont": {"size": 15},
#             "title": {
#                 "text": "I/2efmw",
#                 "font": {
#                     "size": 40,
#                 },
#             },
#         },
#     ),
# )
# fig.show()

# Surface plot ################################################################
# Surface with plotly
print("Generating surface plot")
surface_plot = go.Surface(y=vb_list / uV, x=vmw_list / uV, z=i_average_mesh)
surface_plot_layout = go.Layout(
    title_text=f"$f_{{mw}}={fmw/GHz:.3f}GHz, V_s={Vs/uV:.2f}\mu{{V}}, R={R}\Omega, L={L/uH:.3f}\mu{{H}}$",
    width=1000,
    height=1000,
    margin_pad=8,
    scene=dict(
        camera={
            "eye": {"x": 1.25, "y": -1.25, "z": 1.25},
        },
        yaxis={
            "tickfont": {"size": 15},
            "title": {"text": "Vb (μV)", "font": {"size": 40}},
        },
        xaxis={
            "tickfont": {"size": 15},
            "title": {"text": "Vmw (μV)", "font": {"size": 40}},
        },
        zaxis={
            "tickfont": {"size": 15},
            "title": {
                "text": "I/2efmw",
                "font": {
                    "size": 40,
                },
            },
        },
    ),
)
fig = go.Figure(data=[surface_plot], layout=surface_plot_layout)
fig.show()
