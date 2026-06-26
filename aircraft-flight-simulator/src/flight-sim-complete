import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# ─────────────────────────────────────────────────────────────────────────────
# PARAMETERS  — edit values here
# ─────────────────────────────────────────────────────────────────────────────

# Physical constants
G   = 9.81      # m/s²

# Aircraft
M   = 1200.0    # kg
IYY = 4500.0    # kg·m²  (pitch moment of inertia)
S   = 16.0      # m²     (wing area)
C   = 1.5       # m      (mean aerodynamic chord)

# Aerodynamic coefficients
CL0      =  0.2
CL_ALPHA =  5.5
CL_MAX   =  1.5

CD0 =  0.02
K   =  0.04

CM0        =  0.0
CM_ALPHA   = -0.8
CM_Q       = -14.0
CM_DELTA_E = -1.2

# Altitude hold — outer PID loop
H_COMMAND = 800.0   # m

KPH =  0.005
KIH =  0.0
KDH =  0.0005   # positive → correct anticipatory damping on approach

# Pitch hold — inner PID loop
KP = -1.0
KI = -0.05      # small integral removes steady-state pitch error
KD = -0.8

# Speed hold — throttle PI
KP_T = 300.0    # N / (m/s)
KI_T =  30.0    # N / (m/s²)

# Simulation
DT      = 0.01
T_FINAL = 150.0

# Initial conditions
H0     = 1000.0   # m
V0     = 50.0     # m/s
GAMMA0 = 0.0      # rad
Q0     = 0.0      # rad/s

# ─────────────────────────────────────────────────────────────────────────────
# ISA ATMOSPHERE
# ─────────────────────────────────────────────────────────────────────────────

def isa_density(altitude):
    """Air density via ISA troposphere model (valid 0–11 km)."""
    T0, RHO0, L, R = 288.15, 1.225, 0.0065, 287.05
    alt = np.clip(altitude, 0.0, 11000.0)
    T   = T0 - L * alt
    return RHO0 * (T / T0) ** (G / (R * L) - 1.0)

# ─────────────────────────────────────────────────────────────────────────────
# AERODYNAMICS
# ─────────────────────────────────────────────────────────────────────────────

def compute_forces(V, rho, alpha, q, delta_e):
    """Return lift L (N), drag D (N), pitching moment M_aero (N·m)."""
    alpha = np.clip(alpha, np.radians(-15), np.radians(15))
    q_dyn = 0.5 * rho * V**2 * S

    CL = np.clip(CL0 + CL_ALPHA * alpha, -CL_MAX, CL_MAX)
    CD = CD0 + K * CL**2
    L  = q_dyn * CL
    D  = q_dyn * CD

    q_hat  = q * C / (2.0 * V)
    Cm     = CM0 + CM_ALPHA * alpha + CM_Q * q_hat + CM_DELTA_E * delta_e
    M_aero = q_dyn * C * Cm

    return L, D, M_aero

# ─────────────────────────────────────────────────────────────────────────────
# TRIM SOLVER
# ─────────────────────────────────────────────────────────────────────────────

def compute_trim(V, h):
    """
    Solve V_dot=0, gamma_dot=0, q_dot=0 simultaneously.
    Returns alpha_trim (rad), delta_e_trim (rad), T_trim (N).
    """
    rho   = isa_density(h)
    q_dyn = 0.5 * rho * V**2 * S

    def equations(vars):
        alpha, delta_e, T = vars
        CL = CL0 + CL_ALPHA * alpha
        CD = CD0 + K * CL**2
        Cm = CM0 + CM_ALPHA * alpha + CM_DELTA_E * delta_e

        L  = q_dyn * CL
        D  = q_dyn * CD
        Mq = q_dyn * C * Cm

        V_dot     = (T * np.cos(alpha) - D) / M - G * np.sin(0.0)
        gamma_dot = (L + T * np.sin(alpha)) / (M * V) - G / V
        q_dot     = Mq / IYY
        return [V_dot, gamma_dot, q_dot]

    CL_g  = M * G / q_dyn
    a_g   = (CL_g - CL0) / CL_ALPHA
    de_g  = -(CM0 + CM_ALPHA * a_g) / CM_DELTA_E
    T_g   = q_dyn * (CD0 + K * CL_g**2)

    sol, _, ier, msg = fsolve(equations, [a_g, de_g, T_g], full_output=True)

    if ier != 1:
        print(f"WARNING: Trim did not converge — {msg}")

    alpha_trim, delta_e_trim, T_trim = sol
    res = equations(sol)

    print("Trim solution:")
    print(f"  alpha_trim   = {np.degrees(alpha_trim):.6f} deg")
    print(f"  delta_e_trim = {np.degrees(delta_e_trim):.6f} deg")
    print(f"  T_trim       = {T_trim:.6f} N")
    print(f"  rho at h={h} m = {rho:.5f} kg/m³  (ISA)")
    print("Trim residuals (all should be ~1e-14):")
    print(f"  V_dot={res[0]:.2e}  gamma_dot={res[1]:.2e}  q_dot={res[2]:.2e}")
    print(f"  delta_e_trim {'OK' if abs(np.degrees(delta_e_trim)) <= 30 else 'OUTSIDE ±30 deg'}\n")

    return alpha_trim, delta_e_trim, T_trim

# ─────────────────────────────────────────────────────────────────────────────
# CONTROLLERS
# ─────────────────────────────────────────────────────────────────────────────

class AltitudeController:
    """Outer loop: altitude error → theta_command (rad). Output clipped ±10 deg."""

    def __init__(self):
        self.integral   = 0.0
        self.prev_error = None

    def warm_start(self, h, h_command):
        self.prev_error = h_command - h

    def update(self, h, h_command, dt):
        error           = h_command - h
        self.integral   = np.clip(self.integral + error * dt,
                                  np.radians(-10), np.radians(10))
        if self.prev_error is None:
            self.prev_error = error
        derivative      = (error - self.prev_error) / dt
        self.prev_error = error

        cmd = KPH * error + KIH * self.integral + KDH * derivative
        return np.clip(cmd, np.radians(-10), np.radians(10))


class PitchController:
    """Inner loop: pitch error → elevator deflection (rad). Output clipped ±30 deg."""

    def __init__(self, delta_e_trim):
        self.delta_e_trim = delta_e_trim
        self.integral     = 0.0
        self.prev_error   = None

    def warm_start(self, theta, theta_command_init):
        self.prev_error = theta_command_init - theta

    def update(self, theta, theta_command, dt):
        error           = theta_command - theta
        self.integral   = np.clip(self.integral + error * dt,
                                  np.radians(-20), np.radians(20))
        if self.prev_error is None:
            self.prev_error = error
        derivative      = (error - self.prev_error) / dt
        self.prev_error = error

        de = self.delta_e_trim + KP * error + KI * self.integral + KD * derivative
        return np.clip(de, np.radians(-30), np.radians(30))


class ThrottleController:
    """Speed hold: airspeed error → thrust (N). Capped at 2.5 × T_trim."""

    def __init__(self, T_trim, V_command):
        self.T_trim    = T_trim
        self.V_command = V_command
        self.integral  = 0.0

    def update(self, V, dt):
        error          = self.V_command - V
        self.integral += error * dt
        T = self.T_trim + KP_T * error + KI_T * self.integral
        return np.clip(T, 0.0, 2.5 * self.T_trim)

# ─────────────────────────────────────────────────────────────────────────────
# TRIM & INITIAL STATE
# ─────────────────────────────────────────────────────────────────────────────

alpha_trim, delta_e_trim, T_trim = compute_trim(V0, H0)

h     = H0
V     = V0
gamma = GAMMA0
q     = Q0
theta = alpha_trim + gamma
x     = 0.0

# ─────────────────────────────────────────────────────────────────────────────
# WARM-START CONTROLLERS
# ─────────────────────────────────────────────────────────────────────────────

alt_ctrl = AltitudeController()
alt_ctrl.warm_start(h, H_COMMAND)
theta_command_init = alt_ctrl.update(h, H_COMMAND, DT)

pitch_ctrl = PitchController(delta_e_trim)
pitch_ctrl.warm_start(theta, theta_command_init)

# Re-init alt_ctrl so the first real loop step sees the correct prev_error
alt_ctrl = AltitudeController()
alt_ctrl.warm_start(h, H_COMMAND)

throttle_ctrl = ThrottleController(T_trim, V_command=V0)

# ─────────────────────────────────────────────────────────────────────────────
# SIMULATION LOOP
# ─────────────────────────────────────────────────────────────────────────────

time_arr  = []
theta_arr = []
tcmd_arr  = []
de_arr    = []
h_arr     = []
V_arr     = []
T_arr     = []
gamma_arr = []

t = 0.0

while t < T_FINAL:
    V = max(V, 1.0)

    rho = isa_density(h)

    theta_command = alt_ctrl.update(h, H_COMMAND, DT)
    delta_e       = pitch_ctrl.update(theta, theta_command, DT)
    T             = throttle_ctrl.update(V, DT)

    alpha        = np.clip(theta - gamma, np.radians(-15), np.radians(15))
    L, D, M_aero = compute_forces(V, rho, alpha, q, delta_e)

    V_dot     = (T * np.cos(alpha) - D) / M - G * np.sin(gamma)
    gamma_dot = (L + T * np.sin(alpha)) / (M * V) - G * np.cos(gamma) / V
    q_dot     = M_aero / IYY
    theta_dot = q
    h_dot     = V * np.sin(gamma)
    x_dot     = V * np.cos(gamma)

    V     += V_dot     * DT
    gamma += gamma_dot * DT
    q     += q_dot     * DT
    theta += theta_dot * DT
    h     += h_dot     * DT
    x     += x_dot     * DT

    time_arr.append(t)
    theta_arr.append(np.degrees(theta))
    tcmd_arr.append(np.degrees(theta_command))
    de_arr.append(np.degrees(delta_e))
    h_arr.append(h)
    V_arr.append(V)
    T_arr.append(T)
    gamma_arr.append(np.degrees(gamma))

    t += DT

# ─────────────────────────────────────────────────────────────────────────────
# PLOTS
# ─────────────────────────────────────────────────────────────────────────────

fig, axes = plt.subplots(5, 1, figsize=(11, 18), sharex=True)
fig.suptitle("Aircraft Altitude Hold Simulation", fontsize=13, y=0.998)

axes[0].plot(time_arr, theta_arr, label="Pitch Angle θ")
axes[0].plot(time_arr, tcmd_arr, "--", label="Pitch Command")
axes[0].set_ylabel("Pitch (deg)"); axes[0].legend(); axes[0].grid(True)

axes[1].plot(time_arr, de_arr, color="tab:orange")
axes[1].axhline(np.degrees(delta_e_trim), color="k", linestyle=":",
                linewidth=0.8, label=f"δe_trim = {np.degrees(delta_e_trim):.2f}°")
axes[1].set_ylabel("Elevator (deg)"); axes[1].legend(); axes[1].grid(True)

axes[2].plot(time_arr, h_arr, label="Altitude")
axes[2].axhline(H_COMMAND, color="k", linestyle="--",
                linewidth=0.8, label=f"Command ({H_COMMAND} m)")
axes[2].set_ylabel("Altitude (m)"); axes[2].legend(); axes[2].grid(True)

axes[3].plot(time_arr, V_arr, color="tab:green", label="Airspeed V")
axes[3].axhline(V0, color="k", linestyle="--",
                linewidth=0.8, label=f"V_cmd = {V0} m/s")
axes[3].set_ylabel("Airspeed (m/s)"); axes[3].legend(); axes[3].grid(True)

axes[4].plot(time_arr, T_arr, color="tab:red", label="Thrust T")
axes[4].axhline(T_trim, color="k", linestyle="--",
                linewidth=0.8, label=f"T_trim = {T_trim:.1f} N")
axes[4].set_ylabel("Thrust (N)")
axes[4].set_xlabel("Time (s)"); axes[4].legend(); axes[4].grid(True)

plt.tight_layout()
plt.show()