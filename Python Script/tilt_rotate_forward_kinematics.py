# -*- coding: utf-8 -*-
"""
SHER-3.0 Tilt-Roll Forward Kinematics Visualization

Created on Thu Jun 19 13:26:00 2025

@author: Botao Zhao
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

def ROTz(theta):
    return np.array([[np.cos(theta), -np.sin(theta), 0],
                     [np.sin(theta),  np.cos(theta), 0],
                     [0,              0,             1]])

def Pos(theta, L):
    return np.array([[L*np.cos(theta)],
                     [L*np.sin(theta)],
                     [0]])

def Trans(R, p):
    T = np.vstack((np.hstack((R, p)), [0, 0, 0, 1]))
    return T

def apply_roll_rotation(points, angle_deg, offset_z):
    angle_rad = np.deg2rad(angle_deg)
    R = np.array([[1, 0, 0],
                  [0, np.cos(angle_rad), -np.sin(angle_rad)],
                  [0, np.sin(angle_rad),  np.cos(angle_rad)]])
    rotated = []
    for p in points:
        p_shifted = p - np.array([0, 0, offset_z])
        p_rotated = R @ p_shifted + np.array([0, 0, offset_z])
        rotated.append(p_rotated)
    return rotated

# Fixed Parameters
x_AR_l = -60.48
y_AR = -10
L_AQ = 32.5
L_QR = 42
L_AB = 28
L_BC = 78.5
L_CD = 15
L_DA = 73.5
phi1 = np.pi / 2
phi2 = np.pi * 138 / 180
L_DP = 94.49
phi3 = np.pi * 103.93 / 180
phi4 = np.pi * 31.06 / 180
L_T = 50
d1 = -25.5 # Distance between the center frame of the delta platform to point A

# Create figure and axes
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
fig.subplots_adjust(bottom=0.15)

# Slider axes
slider_ax_s = fig.add_axes([0.2, 0.02, 0.6, 0.03])
slider_ax_roll = fig.add_axes([0.2, 0.06, 0.6, 0.03])

# Create sliders
s_slider = Slider(slider_ax_s, 's', 0, 50, valinit=0)
roll_slider = Slider(slider_ax_roll, 'Roll (°)', -90, 90, valinit=0)

def update_plot(_=None):
    s = s_slider.val
    roll = roll_slider.val
    ax.cla()

    x_AR = x_AR_l + s
    L_AR = np.sqrt(x_AR**2 + y_AR**2)
    alpha3 = np.arccos((L_AR**2 + L_AQ**2 - L_QR**2) / (2 * L_AR * L_AQ))
    alpha4 = np.pi - np.arctan2(y_AR, x_AR)

    theta1 = np.pi - (phi2 + alpha3 + alpha4)
    alpha1 = phi1 - theta1
    L_BD = np.sqrt(L_AB**2 + L_DA**2 - 2 * L_AB * L_DA * np.cos(alpha1))
    beta1 = np.arccos((L_AB**2 + L_BD**2 - L_DA**2) / (2 * L_AB * L_BD))

    delta = np.arccos((L_CD**2 + L_BD**2 - L_BC**2) / (2 * L_CD * L_BD))
    theta2 = phi1 + beta1 - delta

    T_AD = Trans(ROTz(theta1), Pos(theta1, L_DA))
    T_DP = Trans(ROTz(theta2 - phi3 - theta1), Pos(theta2 - phi3 - theta1, L_DP))
    T_PT = Trans(ROTz(np.pi / 2 - phi4 - (theta2 - phi3 - theta1)), Pos(np.pi / 2 - phi4 - (theta2 - phi3 - theta1), 0))
    T_AT = T_AD @ T_DP @ T_PT

    thetat = -np.arccos(T_AT[0][0])
    pt = np.array([T_AT[0][3], 0, T_AT[1][3]])
    p_tip = pt + np.array([np.cos(thetat)*L_T, 0, np.sin(thetat)*L_T])
    p_end = pt + np.array([-np.cos(thetat)*L_T*2, 0, -np.sin(thetat)*L_T*2])

    # Define original points
    A = np.array([0, 0, 0])
    B = A + np.array([0, 0, L_AB])
    D = A + np.array([np.cos(theta1) * L_DA, 0, np.sin(theta1) * L_DA])
    C = D + np.array([np.cos(theta2) * L_CD, 0, np.sin(theta2) * L_CD])
    Q = A + np.array([-np.cos(alpha3 + alpha4) * L_AQ, 0, np.sin(alpha3 + alpha4) * L_AQ])
    R = np.array([x_AR, 0, y_AR])
    P = D + np.array([np.cos(phi3 - theta2) * L_DP, 0, -np.sin(phi3 - theta2) * L_DP])

    # Apply roll rotation to all points
    raw_points = [A, B, C, D, Q, R, P, p_tip, p_end]
    A, B, C, D, Q, R, P, p_tip, p_end = apply_roll_rotation(raw_points, roll, d1)

    # Plot links
    ax.plot([A[0], B[0]], [A[1], B[1]], [A[2], B[2]], 'r', label='AB')
    ax.plot([A[0], D[0]], [A[1], D[1]], [A[2], D[2]], 'g', label='AD')
    ax.plot([B[0], C[0]], [B[1], C[1]], [B[2], C[2]], 'b', label='BC')
    ax.plot([C[0], D[0]], [C[1], D[1]], [C[2], D[2]], 'k', label='CD')
    ax.plot([A[0], Q[0]], [A[1], Q[1]], [A[2], Q[2]], 'g', label='AQ')
    ax.plot([Q[0], R[0]], [Q[1], R[1]], [Q[2], R[2]], 'm', label='QR')
    ax.plot([D[0], P[0]], [D[1], P[1]], [D[2], P[2]], 'c', label='DP')
    ax.plot([C[0], P[0]], [C[1], P[1]], [C[2], P[2]], 'c', label='CP')
    ax.plot([p_tip[0], p_end[0]], [p_tip[1], p_end[1]], [p_tip[2], p_end[2]], 'k', label='Tool')

    # Plot the roll axis (dashed black line)
    x_range = np.array([-100, 200])
    y_fixed = np.array([0, 0])
    z_fixed = np.array([-60, -60])
    ax.plot(x_range, y_fixed, z_fixed, 'k--', linewidth=1, label='Roll Axis')

    # Axis settings
    ax.set_xlim([-100, 200])
    ax.set_ylim([-150, 150])
    ax.set_zlim([-100, 200])
    ax.set_box_aspect([1, 1, 1])
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title(f"Forward Kinematic Plot \n P = ({P[0]:.1f}, {P[1]:.1f}, {P[2]:.1f})")
    ax.legend()

# Initial draw
update_plot()

# Connect sliders
s_slider.on_changed(update_plot)
roll_slider.on_changed(update_plot)

plt.show()
