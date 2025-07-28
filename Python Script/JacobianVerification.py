# -*- coding: utf-8 -*-
"""
Created on Sun Jul 27 15:24:01 2025

@author: Botao
"""

import numpy as np
from sympy import Matrix, zeros, pprint

################### Joint Parameters: ##############################
q1 = 1  # Prismatic
q2 = 2  # Prismatic
q3 = 3  # Prismatic
q4 = 30*np.pi/180  # Revolute
q5 = 10  # Prismatic

################### Joint Velocities: ##############################
q1_dot = 4  # mm/s
q2_dot = 5  # mm/s
q3_dot = 6  # mm/s
q4_dot = 6*np.pi/180  # rad/s
q5_dot = 10  # mm/s
q_dot = Matrix([
        [q1_dot],
        [q2_dot],
        [q3_dot],
        [q4_dot],
        [q5_dot]
])

################### Delta Platform: ###################################

# Robot parameters
rp = 34.4773 # mm platform radius
rb = 60.6927 # mm base radius
l = 68  # mm rod length

# Base joint angles
theta1 = np.pi / 3
theta2 = np.pi
theta3 = 5 * np.pi / 3

# Base joint positions (fixed)
rb1 = np.array([np.cos(theta1)*rb, np.sin(theta1)*rb, 0])
rb2 = np.array([np.cos(theta2)*rb, np.sin(theta2)*rb, 0])
rb3 = np.array([np.cos(theta3)*rb, np.sin(theta3)*rb, 0])

# Platform offsets (relative to center)
rp1 = np.array([np.cos(theta1)*rp, np.sin(theta1)*rp, 0])
rp2 = np.array([np.cos(theta2)*rp, np.sin(theta2)*rp, 0])
rp3 = np.array([np.cos(theta3)*rp, np.sin(theta3)*rp, 0])


L1 = np.array([0,0,q1])
L2 = np.array([0,0,q2])
L3 = np.array([0,0,q3])

e1 = L1 + rb1 - rp1
e2 = L2 + rb2 - rp2
e3 = L3 + rb3 - rp3

h1 = ((e2[0]**2 + e2[1]**2 + e2[2]**2) - (e1[0]**2 + e1[1]**2 + e1[2]**2))/2
h2 = ((e3[0]**2 + e3[1]**2 + e3[2]**2) - (e1[0]**2 + e1[1]**2 + e1[2]**2))/2

x21 = (e2 - e1)[0]
y21 = (e2 - e1)[1]
z21 = (e2 - e1)[2]

x31 = (e3 - e1)[0]
y31 = (e3 - e1)[1]
z31 = (e3 - e1)[2]    

k1 = (y31*z21/y21 - z31)/(x31 - y31*x21/y21)
k2 = (h2 - y31*h1/y21)/(x31 - y31*x21/y21)
k3 = (x31*z21/x21 - z31)/(y31 - x31*y21/x21)
k4 = (h2 - x31*h1/x21)/(y31 - x31*y21/x21)

T1 = k1**2 + k3**2 + 1
T2 = k1*k2 + k3*k4 - e1[0]*k1 - e1[1]*k3 - e1[2]
T3 = k2**2 + k4**2 -2*e1[0]*k2 - 2*e1[1]*k4 - l**2 + e1[0]**2 + e1[1]**2 + e1[2]**2

z = (-T2 + np.sqrt(T2**2 - T1*T3))/T1
x = k1*z + k2
y = k3*z + k4

r = np.array([x, y, z])

l1 = -L1 - rb1 + r + rp1
l2 = -L2 - rb2 + r + rp2
l3 = -L3 - rb3 + r + rp3


################ Roll-Tilt Mechanism: ###########################
    
s = q5
s_dot = q5_dot
psi = q4 # roll angle

# Fixed Parameters
x_AR_max = 60.48
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
d1 = 56.01 # Distance between the center frame of the delta platform to point A
d2 = 192.7
d3 = 25.5

# Angles in radians
phi1 = np.pi / 2
phi2 = np.pi * 138 / 180
phi3 = np.pi * 103.93 / 180
phi4 = np.pi * 31.06 / 180

# Position Analysis
x_AR = x_AR_max - s
L_AR = np.sqrt(x_AR**2 + y_AR**2)
alpha3 = np.arccos((L_AR**2 + L_AQ**2 - L_QR**2) / (2 * L_AR * L_AQ))
alpha4 = np.pi - np.arctan2(y_AR, -x_AR)
alpha = alpha3 + alpha4

theta1 = np.pi - (phi2 + alpha3 + alpha4)
alpha_1 = phi1 - theta1
L_BD = np.sqrt(L_AB**2 + L_DA**2 - 2 * L_AB * L_DA * np.cos(alpha_1))
beta_1 = np.arccos((L_AB**2 + L_BD**2 - L_DA**2) / (2 * L_AB * L_BD))

delta = np.arccos((L_CD**2 + L_BD**2 - L_BC**2) / (2 * L_CD * L_BD))
theta2 = phi1 + beta_1 - delta

# Velocity Analysis
L_JR = np.tan(alpha)*x_AR - y_AR
L_JQ = np.sqrt((L_JR+y_AR)**2 + x_AR**2) - L_AQ

rho_dot = s_dot / L_JR
theta1_dot = -L_JQ / L_AQ * rho_dot

beta_2 = np.arccos((L_BC**2 + L_BD**2 - L_CD**2) / (2 * L_BC * L_BD))
L_IA = L_AB * np.sin(beta_1 + beta_2) / np.sin(alpha_1 + beta_1 + beta_2)
L_ID = L_IA - L_DA

theta2_dot = -L_DA / L_ID * theta1_dot

# Position in {a} frame:
px_a = np.cos(theta1)*L_DA + np.cos(theta2 - phi3)*L_DP
py_a = np.sin(theta1)*L_DA + np.sin(theta2 - phi3)*L_DP

# Position in {p} frame:
px_p = px_a + d2
py_p = -np.sin(psi)*py_a - np.sin(psi)*d3
pz_p = np.cos(psi)*py_a + np.cos(psi)*d3 + d1

# Position in {b} frame:
px_b = px_a + d2 + x
py_b = -np.sin(psi)*py_a - np.sin(psi)*d3 + y
pz_b = np.cos(psi)*py_a + np.cos(psi)*d3 + d1 + z


################## Spacial Jacobian ################

J_l = zeros(3,3)

J_l[0,0] = l1[0]
J_l[0,1] = l1[1]
J_l[0,2] = l1[2]
J_l[1,0] = l2[0]
J_l[1,1] = l2[1]
J_l[1,2] = l2[2]
J_l[2,0] = l3[0]
J_l[2,1] = l3[1]
J_l[2,2] = l3[2]

J_z = zeros(3,3)

J_z[0,0] = l1[2]
J_z[1,1] = l2[2]
J_z[2,2] = l3[2]

J_d = J_l.inv() * J_z


A = np.sin(theta1)*L_DA*L_JQ/(L_AQ*L_JR) - np.sin(theta2-phi3)*L_DP*L_DA*L_JQ/(L_ID*L_AQ*L_JR)
B = -np.cos(theta1)*L_DA*L_JQ/(L_AQ*L_JR) + np.cos(theta2-phi3)*L_DP*L_DA*L_JQ/(L_ID*L_AQ*L_JR)

JacobianTip = zeros(6,5)

JacobianTip[:3, :3] = J_d

JacobianTip[1,3] = -np.cos(q4)*(py_a + d3)
JacobianTip[2,3] = -np.sin(q4)*(py_a + d3)
JacobianTip[3,3] = 1

JacobianTip[0,4] = A
JacobianTip[1,4] = -np.sin(q4)*B
JacobianTip[2,4] = np.cos(q4)*B
JacobianTip[4,4] = -np.cos(q4)*L_DA*L_JQ/(L_ID*L_AQ*L_JR) # Should there be a "-" in front?
JacobianTip[5,4] = -np.sin(q4)*L_DA*L_JQ/(L_ID*L_AQ*L_JR)

Js = JacobianTip  # 6x5 spatial Jacobian in base frame


################# Body Jacobian ###################

def adjoint_pseudo(T):
    R = T[:3, :3]
    return Matrix.vstack(
        Matrix.hstack(R, Matrix.zeros(3)),
        Matrix.hstack(Matrix.zeros(3), R)
    )

# Tbt is known
theta = np.pi/2 - phi4 - phi3 + theta2
Tbt = Matrix([
    [np.cos(theta), -np.sin(theta), 0, px_a+d2+x],
    [-np.sin(theta)*np.sin(psi), -np.cos(theta)*np.sin(psi), -np.cos(psi), -np.sin(psi)*py_a-np.sin(psi)*d3+y],
    [np.sin(theta)*np.cos(psi), np.cos(theta)*np.cos(psi), -np.sin(psi), np.cos(psi)*py_a+np.cos(psi)*d3+d1+z],
    [0, 0, 0, 1]
    ])  # 4x4 transformation from base to tool
Rbt = Tbt[:3, :3]

# Compute Ad_Tbt⁻¹
Ttb = Tbt.inv()
Ad_Ttb = adjoint_pseudo(Ttb)

# Compute body Jacobian
Jb = Ad_Ttb * Js

V_tool_t = Jb*q_dot
V_tool_b = Js*q_dot

# print("\nIn base frame:\n")
# print("x_dot_b = ", V_tool_b[0], "mm/s")
# print("y_dot_b = ", V_tool_b[1], "mm/s")
# print("z_dot_b = ", V_tool_b[2], "mm/s")
# print("theta_x_dot_b = ", V_tool_b[3]*180/np.pi, "deg/s")
# print("theta_y_dot_b = ", V_tool_b[4]*180/np.pi, "deg/s")
# print("theta_z_dot_b = ", V_tool_b[5]*180/np.pi, "deg/s\n")


print("\nIn tool frame:\n")
print("x_dot_t = ", V_tool_t[0], "mm/s")
print("y_dot_t = ", V_tool_t[1], "mm/s")
print("z_dot_t = ", V_tool_t[2], "mm/s")
print("theta_x_dot_t = ", V_tool_t[3]*180/np.pi, "deg/s")
print("theta_y_dot_t = ", V_tool_t[4]*180/np.pi, "deg/s")
print("theta_z_dot_t = ", V_tool_t[5]*180/np.pi, "deg/s\n")



################ Jacobian Inverse ##################

Js_pseudo_inv = (Js.T * Js).inv() * Js.T
Jb_pseudo_inv = (Jb.T * Jb).inv() * Jb.T

q_dot_verify_s = Js_pseudo_inv * V_tool_b
q_dot_verify_b = Jb_pseudo_inv * V_tool_t

print("Joint Velocity from Jacobian Pseudo-Inverse:")
print("q_dot =")
pprint(q_dot_verify_s) # The results should match your q_dot values