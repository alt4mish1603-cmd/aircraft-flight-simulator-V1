# Aircraft Flight Simulation & Autopilot Control System

Built as an independent project in Year 1 to apply flight mechanics theory before formally studying control systems — no template, no tutorial, designed from first principles.

The simulator models longitudinal aircraft dynamics, aerodynamic forces and moments, trim conditions, and a cascaded autopilot system controlling pitch, altitude, and airspeed. It was extended in V2 to include real-time Arduino hardware integration → [aircraft-flight-simulator-V2-arduino](https://github.com/alt4mish1603-cmd/aircraft-flight-simulator-V2-arduino)

---

## Demo

> Altitude climbing from 1000m to 1200m in response to autopilot command.

![Altitude_response_1000-1200](https://github.com/user-attachments/assets/9b6ba876-e813-4f40-9bdb-28edf06ac990)

> Altitude descending from 1000m to 800m.

<img width="1919" height="1079" alt="Altitude_response_1000-800" src="https://github.com/user-attachments/assets/e45e1151-a92c-4d07-b643-5634992a3635" />


---

## Project Structure

```
aircraft-flight-simulator-V1/
│
├── src/
│   ├── aero.py           # Aerodynamic forces and moments
│   ├── atmosphere.py     # ISA atmospheric model
│   ├── config.py         # Aircraft parameters and controller gains
│   ├── controllers.py    # PID autopilot controllers
│   ├── trim.py           # Trim condition solver
│   └── main.py           # Simulation entry point
│
├── plots/                # Simulation output screenshots
└── README.md
```

---

## Physics — What the Simulator Models

### State Vector

The longitudinal aircraft model tracks six state variables at each timestep:

| Variable | Symbol | Description |
|----------|--------|-------------|
| Airspeed | V | Total velocity magnitude |
| Flight Path Angle | γ | Angle of velocity vector to horizontal |
| Pitch Rate | q | Angular velocity about pitch axis |
| Pitch Angle | θ | Nose angle relative to horizontal |
| Altitude | h | Height above sea level |
| Horizontal Position | x | Downrange distance |

The angle of attack follows directly from geometry: `α = θ − γ`

### Equations of Motion

The simulator solves four coupled nonlinear ODEs at each timestep, derived from Newton's second law applied to the aircraft:

```
V̇ = (T·cos α − D) / m  −  g·sin γ
γ̇ = (L + T·sin α) / (mV)  −  g·cos γ / V
q̇ = M_aero / Iyy
θ̇ = q
ḣ = V·sin γ
ẋ = V·cos γ
```

These are integrated using SciPy's `solve_ivp` at each timestep to propagate the aircraft state forward in time.

### Aerodynamic Model

**Lift** uses a linear lift curve slope:

```
CL = CL0 + CL_α · α
L  = ½ρV²S · CL
```

**Drag** uses the parabolic drag polar, combining parasite and induced drag:

```
CD = CD0 + k · CL²
D  = ½ρV²S · CD
```

**Pitching Moment** includes four contributions — zero-lift moment, static stability, pitch damping, and elevator control effectiveness:

```
Cm = Cm0 + Cm_α·α + Cm_q·q̂ + Cm_δe·δe
M  = ½ρV²Sc · Cm
```

Where `q̂ = qc/2V` is the non-dimensionalised pitch rate and `c` is the mean aerodynamic chord.

### ISA Atmosphere

Air density varies with altitude using the International Standard Atmosphere troposphere model:

```
T(h) = 288.15 − 0.0065·h        [K]
p(h) = 101325·(T/288.15)^5.256  [Pa]
ρ(h) = p / (287.05·T)           [kg/m³]
```

At the simulation's starting altitude of 1000m, ρ = 1.112 kg/m³ versus the sea-level value of 1.225 — a 9% reduction that directly affects lift and drag throughout the flight.

### Trim Analysis

Before the simulation begins, the aircraft is initialised at **trim** — the equilibrium condition where all accelerations are zero. Three equations are solved simultaneously using `scipy.optimize.fsolve` for three unknowns:

| Unknown | Physical meaning |
|---------|-----------------|
| α_trim | Angle of attack for lift = weight |
| δe_trim | Elevator deflection for moment balance |
| T_trim | Thrust for thrust = drag |

This ensures the aircraft starts in steady level flight before any autopilot commands are applied. Trim residuals converge to machine precision:

```
V_dot     = 6.63e-14
gamma_dot = 2.78e-17
q_dot     = 5.14e-17
```

These values are effectively zero and confirm a valid trim solution.

---

## Autopilot — Cascaded PID Control

Three nested control loops operate simultaneously:

```
Altitude error → [Altitude PID] → θ_command
Pitch error    → [Pitch PID]    → δe_command  → aircraft
Speed error    → [Throttle PI]  → T_command
```

The cascaded structure exists because altitude and pitch operate on fundamentally different timescales — pitch responds in roughly 1 second, altitude in roughly 100 seconds. Separate loops let each work at its natural frequency without interference.

### Controller Gains

| Controller | Kp | Ki | Kd |
|------------|----|----|----|
| Altitude (outer) | 0.005 | 0.0 | 0.0005 |
| Pitch (inner) | -1.0 | -0.05 | -0.8 |
| Throttle | 300.0 | 30.0 | — |

Negative pitch gains reflect the aerodynamic sign convention — a positive pitch error (nose too low) requires a negative elevator deflection (trailing edge down) to pitch the nose up.

---

## Aircraft Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| m | 1200 kg | Mass |
| Iyy | 4500 kg·m² | Pitch moment of inertia |
| S | 16 m² | Wing area |
| c | 1.5 m | Mean aerodynamic chord |
| CL_α | 5.5 /rad | Lift curve slope |
| CD0 | 0.02 | Zero-lift drag coefficient |
| k | 0.04 | Induced drag factor |
| Cm_α | -0.8 /rad | Static pitch stability derivative |
| Cm_q | -14 | Pitch damping derivative |

---

## How to Run

### Requirements

```
pip install numpy scipy matplotlib
```

### Run the simulation

```
python src/flight-sim-complete.py
```

The simulation initialises at trim, applies the configured altitude command, and generates plots of the aircraft state response on completion.

---

## What I Learned

- Deriving and implementing nonlinear equations of motion from Newton's second law
- Aerodynamic modelling — lift curve slope, parabolic drag polar, pitching moment contributions
- Trim analysis — simultaneous nonlinear solve for equilibrium flight condition using `fsolve`
- Cascaded PID control design and tuning for systems with separated timescales
- ISA atmosphere modelling and its effect on aerodynamic forces across altitude
- Scientific programming in Python with NumPy, SciPy, and Matplotlib

---

## Future Improvements

The natural next steps for this simulator are full 6-DOF dynamics, state-space linearisation around the trim point, eigenvalue-based stability analysis, and LQR optimal control design to replace the hand-tuned PID gains. Wind disturbance modelling and sensor/actuator dynamics would bring it closer to a realistic flight simulation environment.

---


## Author

**Altamish Shohed**
Queen Mary University of London — BEng Aerospace Engineering (Year 1)
