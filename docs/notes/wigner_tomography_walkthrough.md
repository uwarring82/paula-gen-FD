# Wigner tomography of a displaced state — a step-by-step walkthrough

*Companion to the theory note [`strobo_grating_transfer_function.md`](strobo_grating_transfer_function.md)
(which has the full **Notation and terminology** glossary) and the worked-example data in
[`docs/examples/wigner_tomography/`](../examples/wigner_tomography/README.md). This note
describes, step by step, how the digital twin reconstructs the **Wigner function** of a
known displaced (coherent) motional state from **qubit population** measurements. Generated /
reproduced by [`spike/twin_wigner_tomography.py`](../../spike/twin_wigner_tomography.py).*

> ⚠️ **Scope — this does NOT run the stroboscopic grating.** This walkthrough demonstrates the
> *reconstruction mathematics* (populations → $\chi$ → $W$) and the *measurement model* of the
> conditional-displacement / spin-dependent-force readout. It uses the **analytic** $\chi$ of
> the prepared state and an **idealized** readout $P_\downarrow=\tfrac12(1+C\,\mathrm{Re}/\mathrm{Im}\,\chi)$;
> it does **not** run the stroboscopic phase-grating propagator (`engines.strobo_sim`), so the
> grating's own physics — the Floquet comb, the per-cycle kicks, the finite reach
> $\lvert\Delta\beta\rvert\le2\eta$, the multipulse phase-programming criterion, the
> finite-pulse corrections — is *bypassed* (see §"Honesty and scope" below).

## What we are doing, and why

The trapped ion's *motion* is a quantum harmonic oscillator; its state $\rho$ is fully
described by a **Wigner function** $W(\alpha)$ — a quasi-probability map in **phase space**
(the position–momentum plane, dimensionless coordinate $\alpha$). We want to *measure* that
$W$ for a known displaced state, to demonstrate motional-state tomography end to end.

The qubit cannot read $W$ directly. But a population measurement can read the **characteristic
function** $\chi(\beta)=\mathrm{Tr}[\rho\,D(\beta)]$ — the 2-D Fourier transform of $W$ — where
$D(\beta)=e^{\beta a^\dagger-\beta^* a}$ is the displacement operator. So the plan is: *measure
$\chi$, then Fourier-transform to $W$.* (Plain-language definitions of qubit, phase space,
$\chi$, $W$, $D$, coherent state, etc. are in the glossary of the theory note.)

## The measurement sequence

![Wigner-tomography measurement sequence: prepare the state, qubit pi/2, conditional
displacement D(beta), analysis pi/2 with phase phi, detect.](../figures/wigner_tomography_sequence.png)

For each probe displacement $\beta$, **one experimental run** is: (i) cool the ion and
**prepare** the state $|\gamma\rangle$; (ii) a qubit $\pi/2$ pulse — put the qubit in a
superposition; (iii) a **conditional displacement** $D(\beta)$ — a *spin-dependent* kick that
entangles the qubit with the motion (on the apparatus this would be realised by the
spin-dependent-force / grating; **here it is an idealized model — the grating propagator is
not run**, see the scope note); (iv) an
**analysis** $\pi/2$ pulse of phase $\varphi$ — choose which quadrature of $\chi$ to read;
(v) **detect** the qubit (bright/dark). Repeating $M$ times gives the spin-down probability
$P_\downarrow$; scanning the probe $\beta$ and the readout phase $\varphi$ maps out $\chi$.

## Step 1 — Prepare the state

The known input is a **coherent state** $|\gamma\rangle$ with $\gamma = 1.3\,e^{i\pi/4}$
(amplitude $1.3$, phase $45^\circ$). A coherent state is the most "classical" motional state
— a displaced ground state — and its Wigner function is a Gaussian blob centred at
$\alpha=\gamma$. In a real experiment you make it by cooling to the ground state and applying
a calibrated displacement; in the twin it is the known input we will try to recover.

## Step 2 — Measure the characteristic function

The conditional-displacement sequence makes the qubit population a *direct* readout of $\chi$:

$$P_\downarrow(\beta) = \tfrac12\big(1 + C\,\mathrm{Re}\,\chi(\beta)\big)\ \ (\varphi=0),
\qquad
P_\downarrow(\beta) = \tfrac12\big(1 + C\,\mathrm{Im}\,\chi(\beta)\big)\ \ (\varphi=\pi/2),$$

with $C$ the interferometer fringe contrast. We **scan** the probe $\beta$ over a phase-space
grid — the parameters varied are the displacement *magnitude* $\lvert\beta\rvert$ and the
grating optical *phase* $\phi_g=\arg\beta$ (so $\beta=\lvert\beta\rvert e^{i\phi_g}$) — and the
two readout phases $\varphi\in\{0,\pi/2\}$. Each setting is repeated $M$ times to estimate
$P_\downarrow$. The **raw data** is therefore a scan log of *(scanned parameters → measured
population)*; see [`coherent_tomography_raw.dat`](../examples/wigner_tomography/coherent_tomography_raw.dat).

![Raw measured populations: spin-down probability over the scanned displacement plane (two
readout quadratures) and a 1-D cut with shot error bars.](../figures/twin_wigner_raw_data.png)

The fringes in $P_\downarrow$ *are* the displacement information; the 1-D cut shows the
measured population (with binomial shot error bars) against the ideal curve. Here $C=0.90$,
$M=500$ shots per point.

## Step 3 — Reconstruct $\chi$, then the Wigner function

Invert the population model to recover the measured characteristic function (one complex
number per probe $\beta$, from the two quadratures):

$$\hat\chi(\beta) = \frac{(2P_\downarrow^{(\varphi=0)}-1) + i\,(2P_\downarrow^{(\varphi=\pi/2)}-1)}{C}.$$

Then the Wigner function is the 2-D Fourier transform of $\chi$:

$$W(\alpha) = \frac{1}{\pi^2}\int d^2\beta\,\hat\chi(\beta)\,e^{\alpha\beta^*-\alpha^*\beta}
\ \approx\ \frac{1}{\pi^2}\sum_{\rm grid}\hat\chi(\beta)\,e^{\alpha\beta^*-\alpha^*\beta}\,\Delta^2\beta,$$

evaluated over the measured $\beta$-grid (engine function `wigner_from_samples`).

## Step 4 — Validate

The reconstructed $W$ should be a Gaussian blob at $\alpha=\gamma$. The displacement is read
off as the centroid $\langle\alpha\rangle=\gamma$.

![Analytic vs reconstructed Wigner function, and the residual; + marks the input gamma, x the
recovered displacement.](../figures/twin_wigner_reconstruction.png)

**Result:** recovered $\lvert\gamma\rvert = 1.304$ (input $1.300$) and
$\arg\gamma = 0.790$ rad (input $0.785$) — amplitude and phase to better than $1\%$; the peak
Wigner error $\max\lvert W_{\rm rec}-W_{\rm ana}\rvert = 0.035$ is set by the $M=500$ shot
noise. With no shot noise the pipeline is exact (machine precision).

## Honesty and scope

- **The stroboscopic grating is NOT run here.** This validates the reconstruction pipeline
  and the measurement *model*, but $\chi(\beta)$ is the **analytic** characteristic function of
  the prepared coherent state and the readout is the **idealized** $P_\downarrow=\tfrac12(1+C\,\mathrm{Re}/\mathrm{Im}\,\chi)$.
  The grating propagator (`engines.strobo_sim`) is not invoked, so the actual grating physics —
  the Floquet comb, the per-cycle kicks $\beta_k$, the finite reach
  $\lvert\Delta\beta\rvert\le2\eta\approx0.78$, the multipulse phase-programming criterion (§6
  of the theory note) and the finite-pulse corrections (§8) — is bypassed. A **grating-faithful**
  version would instead generate the $\chi$ samples from the propagator running in the SDF /
  conditional-displacement mode (an engine extension not yet built), inheriting those limits.
- **Twin-simulated, not a real measurement.** The raw data is generated by the engine and is
  labelled as such (provenance header + the example README); it lives under `docs/examples/`,
  separate from `sources/data/` (the real DAQ path).
- **Reach (in the real route).** A single two-pulse Ramsey reaches only
  $\lvert\Delta\beta\rvert\le2\eta\approx0.78$ — too small to sample $\chi$ out to its support;
  the grating-faithful route would use the **spin-dependent-force / extended-reach** scheme (a
  calibrated conditional displacement up to $\sim N\eta$), per the theory note (§5–7).
- **Ideal scalings.** The recovered fidelity degrades with lower contrast $C$ and motional
  decoherence (heating, dephasing); the sensitivity / bandwidth / dynamic-range numbers for
  our apparatus are worked out in the tutorial notebook
  [`strobo_grating_tomography.ipynb`](../notebooks/strobo_grating_tomography.ipynb) §7.

## Reproduce

```
python -m spike.twin_wigner_tomography
```

Driver [`spike/twin_wigner_tomography.py`](../../spike/twin_wigner_tomography.py); theory in
the [transfer-function note](strobo_grating_transfer_function.md) (§6–7); data, figures and a
short index in [`docs/examples/wigner_tomography/`](../examples/wigner_tomography/README.md).
Tested in `spike/test_twin_wigner_tomography.py`.
