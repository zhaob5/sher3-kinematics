# -*- coding: utf-8 -*-
"""
SHER-3.0 Tilt-Roll Inverse Kinematics Visualization

Created on Thu Jun 19 11:27:24 2025

@author: Botao Zhao
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# offset_z is the distance between the center frame of the delta platform to point A on the tilt mechanism
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

# Constants
AB = 28
BC = 78.5
CD = 15
DA = 73.5
DP = 94.49
AQ = 32.5
QR = 42
x_AR_l = -60.48
y_AR = -10
L_T = 50
d1 = -25.5 # Distance between the center frame of the delta platform to point A

# Angles in radians
phi1 = np.pi / 2
phi2 = np.pi * 138 / 180
phi3 = np.pi * 103.93 / 180
phi4 = np.pi * 31.06 / 180

#print("theta2 = ", float(theta2*180/np.pi))

theat_tool_max = 150 # deg
theat_tool_min = 110 # deg
theat_tool_init = 130 # deg

# Create figure and axes
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
fig.subplots_adjust(bottom=0.15)

# Slider axes
slider_ax_theta_tool = fig.add_axes([0.2, 0.02, 0.6, 0.03])
slider_ax_roll = fig.add_axes([0.2, 0.06, 0.6, 0.03])

# Create sliders
theta_tool_slider = Slider(slider_ax_theta_tool, 'Tilt (°)', theat_tool_min, theat_tool_max, valinit=theat_tool_init)
roll_slider = Slider(slider_ax_roll, 'Roll (°)', -90, 90, valinit=0)

def update_plot(_=None):
    theta_tool = theta_tool_slider.val * np.pi/180
    roll = roll_slider.val
    ax.cla()
    
    theta2 = phi3 - np.pi + phi4 + theta_tool
    a = CD*np.cos(theta2) - AB*np.cos(phi1)
    b = - CD*np.sin(theta2) + AB*np.sin(phi1)
    gamma = np.arctan2(a, b)
    angle_adc = np.arccos((DA**2+a**2+b**2-BC**2)/(2*DA*np.sqrt(a**2+b**2)))
    theta1 = np.pi/2 - angle_adc + gamma
    
    Px = np.cos(theta1)*DA + np.cos(phi3 - theta2)*DP
    Py = np.sin(theta1)*DA - np.sin(phi3 - theta2)*DP
    
    # Now find the slide bar position s:
    alpha = np.pi - theta1 - phi2 # alpha = alpha3 + alpha4
    
    x_AR = -np.cos(alpha)*AQ - np.sqrt(QR**2 - (np.sin(alpha)*AQ - y_AR)**2)
    
    s = x_AR - x_AR_l
    
    # Define points
    A = np.array([0, 0, 0])
    B = A + np.array([0, 0, AB])
    D = A + np.array([np.cos(theta1) * DA, 0, np.sin(theta1) * DA])
    C = D + np.array([np.cos(theta2) * CD, 0, np.sin(theta2) * CD])
    Q = A + np.array([-np.cos(alpha) * AQ, 0, np.sin(alpha) * AQ])
    R = np.array([x_AR, 0, y_AR])
    P = np.array([Px, 0, Py])
    
    # Tool shaft
    p_tip = P + np.array([np.cos(theta_tool)*L_T*2, 0, np.sin(theta_tool)*L_T*2])
    p_end = P + np.array([-np.cos(theta_tool)*L_T, 0, -np.sin(theta_tool)*L_T])
    
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
    ax.set_title(f"Inverse Kinematic \n s = {s:.1f}, P = ({P[0]:.1f}, {P[1]:.1f}, {P[2]:.1f})")
    ax.legend()
    
# Initial draw
update_plot()

# Connect sliders
theta_tool_slider.on_changed(update_plot)
roll_slider.on_changed(update_plot)

plt.show()