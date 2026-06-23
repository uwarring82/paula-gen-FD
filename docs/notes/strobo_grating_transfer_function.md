# The stroboscopic phase grating as a phase-space probe — transfer function

*Technical note, 2026-06-23. Self-contained derivation of the measurement transfer
function of the Strobo2.0 "active phase grating" and its relation to motional-state
(Wigner) tomography. All leading-order formulas are validated against
[`spike/engines/strobo_sim.py`](../../spike/engines/strobo_sim.py); see §8.*

## 0. Summary

In the impulsive weak-pulse limit the grating maps the motional state onto the qubit
through the **Wigner characteristic function** $\chi(\beta)=\mathrm{Tr}[\rho\,D(\beta)]$.
Which functional of $\chi$ you measure depends on the qubit observable:

| observable | leading-order signal | reconstructs |
|---|---|---|
| spin-flip **probability** $P_\downarrow=\langle A^\dagger A\rangle$ | **double** sum $\sum_{k,k'}\dots\chi(\beta_{k'}-\beta_k)$ | $W$ via a kernel deconvolution; motion-**blind** on the exact strobe |
| spin **coherence** $\langle A\rangle$ (Ramsey w/ reference) | **single** sum $\sum_k e^{-ik\delta\Delta_t}\chi(\beta_k)$ — a discrete FT of $\chi$ | $\chi\!\to\!W$ by Fourier transform |

The clean, exactly-invertible scheme uses the **coherence** with a spin-*dependent*
accumulated displacement; the spin-flip probability (what the bare carrier grating
gives) is the diagnostic/comb side and is motion-blind exactly on resonance.

## 1. Setup and assumptions

- Qubit $\{|\!\uparrow\rangle,|\!\downarrow\rangle\}$, single motional mode of frequency
  $\omega_{\rm lf}$, $a,a^\dagger$. Lamb–Dicke parameter $\eta$.
- The grating applies $N$ identical OC pulses, one per strobe period $\Delta_t$, each of
  area $\theta=\Omega_{\rm strobo}\,\delta t$. **Impulsive weak-pulse limit:**
  $\theta\ll1$ and $\omega_{\rm lf}\delta t\ll1$ (each pulse is short compared to the
  motional period).
- One pulse is, in the Lamb–Dicke approximation, the spin-flip + kick
  $\;\exp[-i\tfrac{\theta}{2}(\sigma_+\!\otimes\! D(i\eta)+\mathrm{h.c.})]$, with
  $D(\alpha)=\exp(\alpha a^\dagger-\alpha^* a)$ the displacement operator.
- A **grating phase** $\phi_g$ rotates the kick: $D(i\eta)\to D(i\eta\,e^{i\phi_g})$.
- A **drive detuning** $\delta$ gives the qubit a phase $e^{-ik\delta\Delta_t}$ between
  cycles. We start in $|\!\uparrow\rangle\otimes\rho$.

In the interaction picture w.r.t. $H_0=\omega_{\rm lf}a^\dagger a$, a kick applied at
cycle $k$ (time $t_k=k\Delta_t$) is rotated in phase space by the free motional
evolution, $D(\alpha)\to D(\alpha\,e^{-i\omega_{\rm lf}t_k})$. Define the **per-cycle
phase slip**

$$\Phi \;=\; \omega_{\rm lf}\Delta_t \;=\; 2\pi\,(f_{\rm lf}\Delta_t),$$

so the sampled displacement at cycle $k$ is

$$\boxed{\;\beta_k \;=\; i\eta\,e^{\,i(\phi_g-k\Phi)}\;}$$

On the **exact strobe** $\Delta_t=2\pi/\omega_{\rm lf}$ ($f_{\rm lf}\Delta_t=1$) we have
$\Phi=2\pi$ and every $\beta_k=i\eta\,e^{i\phi_g}$ — the kick points the same way each
cycle. Off the strobe ($f_{\rm lf}\Delta_t\neq$ integer) the $\beta_k$ fan out around a
circle of radius $\eta$.

## 2. The accumulated flip amplitude

To first order in $\theta$ the qubit picks up at most one flip, summed coherently over
the $N$ cycles. The amplitude operator to reach $|\!\downarrow\rangle$ (acting on the
motion) is

$$A \;=\; -\,i\,\frac{\theta}{2}\sum_{k=0}^{N-1} e^{-ik\delta\Delta_t}\,D(\beta_k).$$

This single operator generates **both** qubit observables.

## 3. Two observables

**(a) Spin-flip probability** $P_\downarrow=\langle A^\dagger A\rangle$. Using the Weyl
identity $D(\mu)^\dagger D(\nu)=D(\nu-\mu)\,e^{\,i\,\mathrm{Im}(\nu\mu^*)}$ and
$\langle D(\beta)\rangle=\chi(\beta)$,

$$\boxed{\;
P_\downarrow(\delta,\phi_g)=\Big(\tfrac{\theta}{2}\Big)^2
\sum_{k,k'=0}^{N-1} e^{\,i(k-k')\delta\Delta_t}\;
e^{\,i\,\mathrm{Im}(\beta_{k'}\beta_k^*)}\;
\chi(\beta_{k'}-\beta_k)\;}
$$

A real, incoherent **double** sum: a fixed linear functional of $\chi$ (hence of $W$),
sampling $\chi$ at the kick *differences* $\beta_{k'}-\beta_k$.

**(b) Spin coherence** $\langle A\rangle$ (measured as Ramsey fringes against a reference
pathway, giving $\mathrm{Re}\langle A\rangle,\ \mathrm{Im}\langle A\rangle$):

$$\boxed{\;
\langle A\rangle(\delta,\phi_g)=-\,i\,\frac{\theta}{2}
\sum_{k=0}^{N-1} e^{-ik\delta\Delta_t}\,\chi(\beta_k)\;}
$$

A **single** sum — a *discrete Fourier transform of $\chi$* sampled along the spiral
$\beta_k$. This is the formula that makes "Fourier-type reconstruction" precise.

Note $\langle A^\dagger A\rangle\neq|\langle A\rangle|^2$ in general (operator variance),
which is exactly why the two rows of the table behave differently.

## 4. Exact strobe — the universal comb

With $\Phi=2\pi$ all $\beta_k=i\eta e^{i\phi_g}$. Then $\beta_{k'}-\beta_k=0$ and
$\mathrm{Im}(\beta_{k'}\beta_k^*)=0$, so the **probability** collapses to

$$P_\downarrow(\delta)=\Big(\tfrac{\theta}{2}\Big)^2
\left|\sum_{k=0}^{N-1}e^{-ik\delta\Delta_t}\right|^2
=\Big(\tfrac{\theta}{2}\Big)^2
\left|\frac{\sin(N\delta\Delta_t/2)}{\sin(\delta\Delta_t/2)}\right|^2 .$$

This **Dirichlet/Fejér comb** is independent of $\eta$ and of the motional state — the
spin-flip probability is *motion-blind on resonance* (teeth at $\delta=k/\Delta_t=k\,
f_{\rm lf}$, width $\sim1/(N\Delta_t)$). It is the origin of the heterodyne beat: detune
by $f_{\rm IF}$ and the cycle-domain population nutates with half-beat
$1/(2 f_{\rm IF}\Delta_t)$ ([`heterodyne_beat`](../../spike/twin_strobo.py)).

The **coherence**, by contrast, does *not* collapse:
$\langle A\rangle=-i\tfrac{\theta}{2}\,\chi(i\eta e^{i\phi_g})\,S_N(\delta)$ — it still
reads $\chi$ at one phase-space point. So the motion is invisible to $P_\downarrow$ but
visible to the coherence, even on the exact strobe.

## 5. Off-strobe / phase-varying grating — sampling $\chi(\beta)$

Off the exact strobe the $\beta_k$ rotate and both kernels become motion-sensitive:

- **Coherence:** $\langle A\rangle=-i\tfrac{\theta}{2}\sum_k e^{-ik\delta\Delta_t}
  \chi\!\big(i\eta e^{i(\phi_g-k\Phi)}\big)$ — a discrete FT of $\chi$ over the spiral.
  The grating phase $\phi_g$ rotates the sampling ring; the detuning $\delta$ supplies the
  Fourier phase.
- **Probability:** the §3(a) double sum; its $\eta$-dependence comes entirely from the
  off-diagonal $k\neq k'$ terms (zero on the strobe). This is the ~13% leading-order
  signal seen when the strobe is mistuned 3% (§8).

**Reach.** A single carrier grating only samples $|\beta|=\eta$ (a thin ring — $\eta$ is
small). To probe larger $|\beta|$ (finer phase-space structure) the per-cycle kicks must
**add coherently**, which requires a spin-*dependent* displacement (the MS/SDF / sideband
configuration), giving a net $\beta_{\rm tot}=\sum_k\delta\beta_k$ up to $\sim N\eta$ and
the exact coherence $\langle\sigma_+\rangle=\tfrac12\chi(\beta_{\rm tot})$. Rastering
$\phi_g$ (direction) and $N$ or $\delta$ (magnitude) maps $\chi$ over the plane.

## 6. From $\chi$ to the Wigner function

$\chi$ and $W$ are a 2-D Fourier pair,

$$\chi(\beta)=\int d^2\alpha\,W(\alpha)\,e^{\beta\alpha^*-\beta^*\alpha},\qquad
W(\alpha)=\frac{1}{\pi^2}\int d^2\beta\,\chi(\beta)\,e^{\alpha\beta^*-\alpha^*\beta}.$$

Hence three equivalent reconstruction routes from the measured kernel:

1. **Characteristic-function tomography** — sample $\chi(\beta)$ (coherence/SDF) and 2-D
   FFT $\to W$. *Direct route for this grating* (Banaszek–Wódkiewicz / Leibfried family).
2. **Quadrature / inverse Radon** — by the Fourier-slice theorem a radial slice
   $\chi(s e^{i\phi})$ is the 1-D FT of the marginal $\mathrm{pr}(x_\phi)$; 1-D-invert
   each slice to marginals, then filtered back-projection $\to W$ (Vogel–Risken). The
   *dual* of route 1; "P_flip = Radon of $W$" holds only approximately (small-$\eta$ gives
   quadrature *moments*, not the full marginal; large $\eta$ makes the kernel nonlinear).
3. **Displaced parity (no transform)** — $\langle\Pi(\beta)\rangle=\tfrac{\pi}{2}W(\beta)$
   gives $W$ *directly* if the readout maps onto motional parity (a full sideband π map)
   rather than the coherence (Leibfried; Lutterbach–Davidovich).

## 7. Finite-pulse and non-ideal corrections

Beyond the impulsive limit, replace the pulse area by its detuned, finite-width form

$$\frac{\theta}{2}\;\to\;\frac{\Omega_{\rm strobo}}{\Omega_{\rm eff}}
\sin\!\Big(\frac{\Omega_{\rm eff}\,\delta t}{2}\Big),\qquad
\Omega_{\rm eff}=\sqrt{\Omega_{\rm strobo}^2+\delta_{\rm AC}^2},$$

and include the motion *during* the pulse through a sampling factor
$\mathrm{sinc}(\omega_{\rm lf}\delta t/2)$ that smears the effective $\eta$
($\omega_{\rm lf}\delta t\approx0.16$ rad at $\delta t=0.02\,\mu$s — modest but nonzero).
Finite $N$ and decoherence cap the reachable $|\beta|_{\max}$, convolving the
reconstructed $W$ with a point-spread of width $\sim1/|\beta|_{\max}$ — so the
**noise/filter-function** analysis sets the tomographic resolution.

## 8. Non-perturbative regime and numerical validation

When $\theta$ is not small (or many phonons participate) the weak-pulse expansion fails
and one diagonalizes the full Floquet step $U_{\rm cycle}=e^{-iH_0\Delta_t}\,
e^{-iH_{\rm pulse}\delta t}$ — exactly what `strobo_sim` does. Checks performed:

- **Comb (§4):** `strobo_detuning_scan` gives full-contrast teeth at $\delta=k f_{\rm lf}$,
  width $\approx1/(N\Delta_t)=26$ kHz; **$\eta$-independent** on the exact strobe.
- **Double-sum kernel (§3a):** the analytical $P_\downarrow$ for $|0\rangle$
  ($\chi(\beta)=e^{-|\beta|^2/2}$) matches `strobo_detuning_scan` to **0.05–0.2 %** in the
  weak/off-strobe regime, the residual $\propto\theta$ (validates the leading order).
- **Motion-blindness:** on the exact strobe $P_\downarrow$ is $\eta$-independent; a 3 %
  strobe mistune makes it $\eta$-sensitive at $\sim$13 % (leading order, not saturation).
- **Heterodyne (§4):** `strobo_population_vs_cycles` nutates with half-beat
  $1/(2f_{\rm IF}\Delta_t)$ (50/100/200 kHz $\to$ 13/7/3 cycles).

## 9. References (standard motional-state tomography — verify before citing)

- K. Vogel, H. Risken, *Phys. Rev. A* **40**, 2847 (1989) — quadrature/Radon tomography.
- D. Leibfried *et al.*, *Phys. Rev. Lett.* **77**, 4281 (1996) — trapped-ion motional
  Wigner function via displaced parity.
- L. G. Lutterbach, L. Davidovich, *Phys. Rev. Lett.* **78**, 2547 (1997) — direct Wigner
  measurement from interference.
- K. Banaszek, K. Wódkiewicz, *Phys. Rev. Lett.* **76**, 4344 (1996) — direct probing of
  the characteristic function / quasiprobabilities.
- S. Wallentowitz, W. Vogel, *Phys. Rev. Lett.* **75**, 2932 (1995); P. J. Bardroff
  *et al.*, *Phys. Rev. Lett.* **77**, 2198 (1996) — pattern-function reconstruction.

See also [`docs/LOGBOOK.md`](../LOGBOOK.md) (2026-06-23 entries) and the memory note
`project-strobo-grating-receiver`.
