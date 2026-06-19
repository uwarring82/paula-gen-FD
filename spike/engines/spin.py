"""
Qubit spin STATE and its coherent operations — the object that flows through the
integrated twin (prepare -> drive -> detect).

The state is the Bloch vector (x, y, z) of the |down>=|F=3,mF=+3> / |up>=|F=2,mF=+2>
qubit, with the convention z = +1 -> |down> (P_up = 0), z = -1 -> |up> (P_up = 1),
P_up = (1 - z)/2. Operations rotate it:

  * prepare(eps)        — prepared |down> with residual |up> population eps (prep error)
  * pulse(.., Omega, delta, t, phase) — a (possibly detuned) drive pulse: rotation by
    Omega_eff*t about (Omega cos phi, Omega sin phi, delta)/Omega_eff
  * free(.., delta, t)  — free precession about z by 2 pi delta t

All frequencies in Hz, times in s (the 2pi enters the rotation angle). The drive
detuning is supplied by the caller as delta = delta_set + delta_ACZeeman + delta_noise,
so AC-Zeeman shifts and (quasi-static) dephasing are just contributions to delta.
"""
from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass
class Bloch:
    x: float = 0.0
    y: float = 0.0
    z: float = 1.0          # |down> by default


def prepare(eps_prep: float = 0.0) -> Bloch:
    """Prepared |down> with residual |up> population eps_prep (P_up = eps_prep)."""
    return Bloch(0.0, 0.0, 1.0 - 2.0 * eps_prep)


def p_up(s: Bloch) -> float:
    """Excited-state probability P(|up>) = (1 - z)/2."""
    return 0.5 * (1.0 - s.z)


def _rodrigues(v, k, theta):
    """Rotate vector v about unit axis k by angle theta (Rodrigues' formula)."""
    vx, vy, vz = v
    kx, ky, kz = k
    ct, st = math.cos(theta), math.sin(theta)
    kdv = kx * vx + ky * vy + kz * vz
    cx, cy, cz = ky * vz - kz * vy, kz * vx - kx * vz, kx * vy - ky * vx   # k x v
    return (vx * ct + cx * st + kx * kdv * (1.0 - ct),
            vy * ct + cy * st + ky * kdv * (1.0 - ct),
            vz * ct + cz * st + kz * kdv * (1.0 - ct))


def pulse(s: Bloch, rabi_hz: float, detuning_hz: float, t_s: float, phase: float = 0.0) -> Bloch:
    """Apply a drive pulse of Rabi frequency rabi_hz (= Omega/2pi), detuning_hz and
    duration t_s, with drive phase `phase` [rad]. Rotates the Bloch vector by
    Omega_eff*t about the generalised-Rabi axis."""
    om, dl = rabi_hz, detuning_hz
    eff = math.hypot(om, dl)
    if eff == 0.0:
        return Bloch(s.x, s.y, s.z)
    k = (om * math.cos(phase) / eff, om * math.sin(phase) / eff, dl / eff)
    theta = 2.0 * math.pi * eff * t_s
    x, y, z = _rodrigues((s.x, s.y, s.z), k, theta)
    return Bloch(x, y, z)


def free(s: Bloch, detuning_hz: float, t_s: float) -> Bloch:
    """Free precession about z by 2 pi * detuning * t."""
    phi = 2.0 * math.pi * detuning_hz * t_s
    c, sn = math.cos(phi), math.sin(phi)
    return Bloch(c * s.x - sn * s.y, sn * s.x + c * s.y, s.z)


def pi_pulse_time(rabi_hz: float) -> float:
    """On-resonance pi-pulse duration 1/(2 rabi_hz)."""
    return 1.0 / (2.0 * rabi_hz)
