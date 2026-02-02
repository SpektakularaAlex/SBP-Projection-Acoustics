# SBP-Projection Acoustics

Stable high-order finite-difference solvers for the **acoustic wave equation** using the **SBP-Projection method**.  
This repository collects the code, experiments, and selected figures from our course project in scientific computing for PDEs.

📄 **Project report (PDF):** `report/Projekt_1_Bervet_PDE.pdf`

---

## Overview

We study and implement provably stable SBP-Projection discretizations for the acoustic wave equation, focusing on:

- **Well-posed boundary conditions** (Dirichlet & characteristic)
- **High-order SBP finite differences** (central and upwind variants)
- **Energy stability via the energy method**
- **Eigenvalue spectra and RK4 CFL estimates**
- Time integration using **4th-order Runge–Kutta (RK4)**

The implementation is primarily in **MATLAB**.

---

## Model problem

### 1D acoustic wave equation (bounded domain)

We consider the 1D acoustic system on a bounded interval with an absorption term:

- State: \(u = (p, v)^T\), where \(p\) is pressure and \(v\) is particle velocity  
- Material parameters: density \(\rho(x)\), sound speed \(c(x)\), absorption \(\beta(x)\ge 0\)

Boundary conditions are written in the general form
\[
p + a_\ell v = 0 \ \text{at } x=x_\ell,
\qquad
p + a_r v = 0 \ \text{at } x=x_r.
\]

A key result from the energy method is the well-posedness condition:
- \(a_\ell \ge 0\)
- \(a_r \le 0\)

Characteristic boundary conditions correspond to choosing incoming characteristic variables to be zero, which yields:
\[
a_\ell = \rho c, \qquad a_r = -\rho c.
\]

(See report for full derivations and energy estimates.)

---

## Method: SBP-Projection

Boundary conditions are imposed strongly using a projection operator of the form
\[
P = I - \bar{H}L^T (L\bar{H}L^T)^{-1} L,
\]
so the semi-discrete system becomes
\[
u_t = -P\,C^{-1}(D_x + D)\,P u.
\]

We test both:
- **central SBP** derivative operators, and
- **upwind SBP** operators \(D^\pm\),

and verify that the resulting schemes satisfy a discrete energy estimate.

---

## Numerical stability & CFL (RK4)

To assess time-step restrictions, we rewrite the semi-discrete system as:
\[
u_t = M u,
\]
compute eigenvalues of \(M\), and use the RK4 stability radius \(R \approx 2.78\) to estimate a stable time step:
\[
k|\lambda| \le R \quad \Rightarrow \quad \alpha = \frac{k}{h} \approx \frac{2.78}{|\lambda|h}.
\]

In the report, this yields a CFL estimate of approximately **α ≈ 1.593** for the tested SBP-Projection setup.

---

## Results (figures)

All figures shown in this README are expected to live in `figures/`.

> Add your images to `figures/` and update the filenames below.

### Eigenvalues of the semi-discrete operator

<p align="center">
  <img src="figures/eigs_dirichlet_bc1.png" width="31%">
  <img src="figures/eigs_dirichlet_bc2.png" width="31%">
  <img src="figures/eigs_characteristic.png" width="31%">
</p>
<p align="center"><em>Eigenvalue plots of hM for different boundary conditions.</em></p>

### Example wave propagation / solution snapshots (optional)

<p align="center">
  <img src="figures/wave_solution.png" width="75%">
</p>

---

