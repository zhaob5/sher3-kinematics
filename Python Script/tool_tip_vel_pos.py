import numpy as np
import matplotlib.pyplot as plt

# Constants
s_dot = 1
s_values = np.linspace(0, 50, 500)  # np array from 0 to 50
psi = 30*np.pi/180
psi_dot = 1
velocity_ratios = []
position_x = []
position_y = []
velocity_x = []
velocity_y = []
velocity_xp = []
velocity_yp = []
velocity_zp = []

# Link lengths and angles
x_AR_l = 60.48
# x_AR_l = 63
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
# L_DP = 74.8
phi3 = np.pi * 103.93 / 180
phi4 = np.pi * 31.06 / 180
L_T = 50

d1 = 56.01 # mm
d2 = 192.700 # mm
d3 = 25.500 # mm

for s in s_values:
    try:
        # Position Analysis
        x_AR = x_AR_l - s
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
        
        px = np.cos(theta1)*L_DA + np.cos(theta2 - phi3)*L_DP
        py = np.sin(theta1)*L_DA + np.sin(theta2 - phi3)*L_DP
        
        px_dot = -np.sin(theta1)*L_DA*theta1_dot - np.sin(theta2 - phi3)*L_DP*theta2_dot
        py_dot = np.cos(theta1)*L_DA*theta1_dot + np.cos(theta2 - phi3)*L_DP*theta2_dot

        vel_ratio = (theta2_dot*180/np.pi) / s_dot
        velocity_ratios.append(vel_ratio)
        
        position_x.append(px)
        position_y.append(py)
        
        velocity_x.append(px_dot)
        velocity_y.append(py_dot)
        
        v_xp = px_dot
        v_yp = -np.sin(psi)*py_dot - np.cos(psi)*psi_dot*(py + d3)
        v_zp = np.cos(psi)*py_dot - np.sin(psi)*psi_dot*(py + d3)
        
        velocity_xp.append(v_xp)
        velocity_yp.append(v_yp)
        velocity_zp.append(v_zp)
        
    except (ValueError, ZeroDivisionError):  # Handle domain errors in arccos or division
        velocity_ratios.append(np.nan)
        velocity_x.append(np.nan)
        velocity_y.append(np.nan)

# Plotting

plt.subplot(1, 3, 1)
plt.plot(s_values,  velocity_xp)
plt.title('vx')

plt.subplot(1, 3, 2)
plt.plot(s_values,  velocity_yp)
plt.title('vy')

plt.subplot(1, 3, 3)
plt.plot(s_values,  velocity_zp)
plt.title('vz')

plt.tight_layout()
plt.show()

# plt.figure()
# plt.plot(s_values, velocity_ratios)
# plt.xlabel('Stroke Position (mm)')
# plt.ylabel('Velocity ratio (deg/mm)')
# plt.title('Velocity Ratio')
# plt.grid(True)
# plt.show()

# plt.figure()
# plt.axis('equal')
# plt.plot(position_x, position_y)
# plt.xlabel('x (mm)')
# plt.ylabel('y (mm)')
# plt.title('Tool Tip Position in frame {a}')
# plt.grid(True)
# plt.show()

# plt.figure()
# plt.plot(s_values, velocity_x)
# plt.xlabel('s')
# plt.ylabel('velocity_x')
# plt.title('Velocity_x')
# plt.grid(True)
# plt.show()

# plt.figure()
# plt.plot(s_values, velocity_y)
# plt.xlabel('s')
# plt.ylabel('velocity_y')
# plt.title('Velocity_y')
# plt.grid(True)
# plt.show()
