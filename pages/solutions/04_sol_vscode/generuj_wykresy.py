from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "wykresy"
OUT.mkdir(parents=True, exist_ok=True)

v0 = 37.0
alpha_deg = 85.0
g = 9.81
mass = 0.145
diameter = 0.074
air_density = 1.225
drag_coefficient = 0.47
area = np.pi * diameter**2 / 4
c_drag = 0.5 * air_density * drag_coefficient * area
k_drag = c_drag / mass
dt = 0.001
alpha = np.deg2rad(alpha_deg)
initial = np.array([0.0, 0.0, v0 * np.cos(alpha), v0 * np.sin(alpha)])

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.25,
})


def pochodne(stan, z_oporem):
    x, y, vx, vy = stan
    szybkosc = np.hypot(vx, vy)
    ax = -k_drag * szybkosc * vx if z_oporem else 0.0
    ay = -g - k_drag * szybkosc * vy if z_oporem else -g
    return np.array([vx, vy, ax, ay])


def krok_rk4(stan, krok, z_oporem):
    k1 = pochodne(stan, z_oporem)
    k2 = pochodne(stan + krok * k1 / 2, z_oporem)
    k3 = pochodne(stan + krok * k2 / 2, z_oporem)
    k4 = pochodne(stan + krok * k3, z_oporem)
    return stan + krok * (k1 + 2 * k2 + 2 * k3 + k4) / 6


def symuluj(z_oporem):
    czasy = [0.0]
    stany = [initial.copy()]
    czas = 0.0
    stan = initial.copy()

    while czas < 20.0:
        poprzedni = stan
        nastepny = krok_rk4(stan, dt, z_oporem)
        if nastepny[1] < 0.0:
            udzial = poprzedni[1] / (poprzedni[1] - nastepny[1])
            czasy.append(czas + udzial * dt)
            stany.append(poprzedni + udzial * (nastepny - poprzedni))
            break
        czas += dt
        stan = nastepny
        czasy.append(czas)
        stany.append(stan.copy())
    else:
        raise RuntimeError("Nie znaleziono chwili lądowania w ciągu 20 s.")

    return np.asarray(czasy), np.asarray(stany)


def zapisz(fig, nazwa):
    fig.tight_layout()
    fig.savefig(OUT / f"{nazwa}.svg", bbox_inches="tight")
    fig.savefig(OUT / f"{nazwa}.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


t_drag, s_drag = symuluj(True)
t_vacuum = np.linspace(0.0, 2 * initial[3] / g, 1000)
s_vacuum = np.column_stack((
    initial[2] * t_vacuum,
    initial[3] * t_vacuum - 0.5 * g * t_vacuum**2,
    np.full_like(t_vacuum, initial[2]),
    initial[3] - g * t_vacuum,
))

blue, orange = "#2364aa", "#d97732"

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(s_vacuum[:, 0], s_vacuum[:, 1], color=blue, ls="--", lw=2,
        label="Bez oporu")
ax.plot(s_drag[:, 0], s_drag[:, 1], color=orange, lw=2.2,
        label="Opór kwadratowy")
ax.set(title="Tor rzutu ukośnego", xlabel="Położenie poziome x (m)",
       ylabel="Wysokość y (m)")
ax.legend()
zapisz(fig, "tor_rzutu")

fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
axes[0].plot(t_vacuum, s_vacuum[:, 0], color=blue, ls="--", label="Bez oporu")
axes[0].plot(t_drag, s_drag[:, 0], color=orange, label="Opór kwadratowy")
axes[0].set_ylabel("Położenie x (m)")
axes[0].legend()
axes[1].plot(t_vacuum, s_vacuum[:, 1], color=blue, ls="--", label="Bez oporu")
axes[1].plot(t_drag, s_drag[:, 1], color=orange, label="Opór kwadratowy")
axes[1].set(xlabel="Czas t (s)", ylabel="Wysokość y (m)")
axes[1].legend()
fig.suptitle("Położenie w funkcji czasu")
zapisz(fig, "polozenie_czas")

fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
axes[0].plot(t_vacuum, s_vacuum[:, 2], color=blue, ls="--", label="Bez oporu")
axes[0].plot(t_drag, s_drag[:, 2], color=orange, label="Opór kwadratowy")
axes[0].set_ylabel("Prędkość vx (m/s)")
axes[0].legend()
axes[1].plot(t_vacuum, s_vacuum[:, 3], color=blue, ls="--", label="Bez oporu")
axes[1].plot(t_drag, s_drag[:, 3], color=orange, label="Opór kwadratowy")
axes[1].set(xlabel="Czas t (s)", ylabel="Prędkość vy (m/s)")
axes[1].legend()
fig.suptitle("Składowe prędkości w funkcji czasu")
zapisz(fig, "predkosc_czas")

energy_vacuum = 0.5 * mass * (s_vacuum[:, 2]**2 + s_vacuum[:, 3]**2) + mass * g * s_vacuum[:, 1]
energy_drag = 0.5 * mass * (s_drag[:, 2]**2 + s_drag[:, 3]**2) + mass * g * s_drag[:, 1]
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(t_vacuum, energy_vacuum, color=blue, ls="--", lw=2,
        label="Energia mechaniczna bez oporu")
ax.plot(t_drag, energy_drag, color=orange, lw=2,
        label="Energia mechaniczna z oporem")
ax.set(title="Zmiana energii mechanicznej", xlabel="Czas t (s)",
       ylabel="Energia (J)")
ax.legend()
zapisz(fig, "energia")

apex = int(np.argmax(s_drag[:, 1]))
print(f"v0 = {v0:.3f} m/s, alpha = {alpha_deg:.1f} deg")
print(f"Pole przekroju = {area:.6f} m^2, c_drag = {c_drag:.6f} kg/m")
print(f"Bez oporu: T = {t_vacuum[-1]:.3f} s, h_max = {s_vacuum[:, 1].max():.3f} m, R = {s_vacuum[-1, 0]:.3f} m")
print(f"Z oporem: T = {t_drag[-1]:.3f} s, h_max = {s_drag[apex, 1]:.3f} m, R = {s_drag[-1, 0]:.3f} m")
print(f"Wykresy zapisano w: {OUT}")
