# Robust stability of nonlinear biofuel models via µ / γ bounds

MATLAB code for Sections 5.2–5.9 of the research report *Integrated Structural and Probabilistic Approaches for Biological and Epidemiological Systems* (ERC project INSPIRE).

The scripts certify robust stability of two microbial biofuel production models under ±20% kinetic-parameter uncertainty, using the structured singular value µ (box uncertainty, ∞-norm) and its Euclidean-norm generalisation γ (ball uncertainty, 2-norm). The four upper bounds come from Álamo & Dormido (Int. J. Robust Nonlinear Control, 2001): Property 1, Property 2, the recursive bound and the Perron bound, with MATLAB's `mussv` as an external reference. Section 4 of the report reproduces the crane example from that paper and has no code here.

## Models

- **Model A** (open loop): states *(n, p, b)* — cell density, pump, biofuel. Six uncertain parameters: `alpha_n, delta_n, gamma_p, beta_p, alpha_b, delta_b`.
- **Model B** (feedback): adds the repressor *R*. Nine uncertain parameters: the six above plus `k_p, k_b, gamma_R`. The repressor steady state *R̄* depends only on fixed parameters, so it is a constant and is not swept.

Every parameter is perturbed by ±20% (`rho = 0.20`), i.e. `theta = theta_nom (1 + rho*delta)` with `delta` in [-1, 1].

## Four ways to build the LFT

| Method | Idea | Script(s) |
|---|---|---|
| 2 — BDC | Write the Jacobian as `J = B diag(D) C` and treat each entry `D_i` as an uncertain channel, with ranges swept over the parameter box. | `method2_bdc/` |
| 3 — Equilibrium as parameters | Kinetic parameters plus the equilibrium coordinates `n̄, p̄, b̄` as extra channels, anchored at the midpoint of the swept equilibrium box. | `method3_eqparam/` |
| 4 — Total derivative | Only kinetic channels; the motion of the equilibrium is folded into each channel through the implicit function theorem. | baseline run of `taylor_remainder/` |

The ±20% box (`delta` in [-1, 1]^m) is *covered* in the ∞-norm when `1/µ ≥ 1`, and in the 2-norm when `1/γ ≥ sqrt(m)`.

## Repository layout

```
method2_bdc/                 BDC_newrange_B.m
method3_eqparam/             biofuelA.m, biofuelB.m
vertex_montecarlo_checks/    modelA_openloop_check.m, modelB_closedloop_check.m
physicality_radius/          empirical_radius_l2.m
taylor_remainder/            biofuel_lft_remainder.m, log_uncertainty.m
```

## Requirements

- MATLAB with the Optimization Toolbox (`fsolve`), Robust Control Toolbox (`mussv`) and Symbolic Math Toolbox (the `taylor_remainder` scripts).
- YALMIP with the SeDuMi solver, used by Property 2.
- The four bound functions on the MATLAB path: `mu_upper_bound_property1`, `mu_upper_bound_property2`, `mu_upper_bound_recursive`, `gamma_perron_bound`.

## Experiments

### 1. Vertex and Monte Carlo checks — `vertex_montecarlo_checks/`
*Report: equilibrium intervals (Table 12), BDC channel ranges (Table 15), vertex / Monte Carlo paragraph in Section 5.6.*

`modelA_openloop_check.m` (Model A, 2⁶ = 64 vertices) and `modelB_closedloop_check.m` (Model B, 2⁹ = 512 vertices) do the groundwork for the other scripts. They re-solve the equilibrium at the corners of the parameter box and on a dense grid (A) or a random interior sample (B), and report:

- the **equilibrium box** for `n̄, p̄, b̄` (and whether the extremes sit at the vertices);
- the **BDC D-box**, i.e. the range of every Jacobian channel `D_i` (7 channels for A, 8 for B);
- a **stability comparison** on the exact Jacobian `J = B diag(D) C`: all vertices, uniform Monte Carlo over the physical box, Monte Carlo over the independent D-box, and a radius sweep that scales the box by α (both sampled and vertex-based).

### 2. Method 2 — BDC decomposition — `method2_bdc/BDC_newrange_B.m`
*Report: Section 5.2.*

Builds the scaled `B, C` matrices, anchors the Jacobian at the midpoint of the D-box and forms `M(s) = C (sI − J0)⁻¹ B diag(W)` with one real scalar block per channel (q = 7 for A, 8 for B). It sweeps ω over [10⁻³, 10³] and prints `sup µ` (Property 1) and `sup γ` (Property 2, recursive, Perron) with the corresponding radii. Select the model with `model_selection` at the top.

The D-box (`D_lo`, `D_hi`) is hard-coded from the "Final D-box" printed by the check scripts.

### 3. Method 3 — equilibrium as parameters — `method3_eqparam/`
*Report: Section 5.3.*

`biofuelA.m` and `biofuelB.m` add the equilibrium coordinates `n̄, p̄, b̄` as extra channels, each with its own radius taken from the swept equilibrium box, and anchor the Jacobian at the box midpoint. Each sensitivity matrix `J_k` is factored by SVD into `L_k R_k` (parameters that move two Jacobian entries, such as `p̄`, become repeated blocks of size 2), giving m = 9 channels / q = 10 for Model A and m = 12 / q = 14 for Model B. The scripts print the factorisation error (machine precision), the pre-sweep check at ω = 0, and the µ / γ results over 15 frequencies.

The equilibrium intervals (`n_bar_lo`, `n_bar_hi`, …) are hard-coded from the check scripts.

### 4. Method 4 — total derivative, and the Taylor-form remainder — `taylor_remainder/`
*Report: Sections 5.4–5.6 (total derivative), 5.8 (remainder), 5.9 (logarithmic uncertainty).*

`biofuel_lft_remainder.m` builds the total-derivative LFT symbolically (Jacobians, implicit-function-theorem sensitivities, SVD rank factorisation) and then bounds the error of the first-order approximation. Set `model` to `'A'`, `'B'` or `'both'`.

The remainder is estimated two ways and combined entrywise:

- **Way A — residual sampling:** the true Jacobian (equilibrium re-solved with `fsolve`) is compared with the affine surrogate over box corners and Monte Carlo draws.
- **Way B — analytic curvature:** the total-derivative operator is applied a second time to obtain a Hessian-type bound.

The remainder matrix is folded back into the LFT as one rank-1 real-scalar block per nonzero entry, and `run_mu` evaluates µ and the γ bounds twice. **The baseline run (`A (baseline)`, `B (baseline)`) is the plain total-derivative LFT of Method 4**; the second run (`with remainder`) shows how much of that margin survives once the linearisation error is included. The frequency grid is set in `run_mu` (`w = logspace(-3,3,N)`); the report's total-derivative table uses 25 points.

`log_uncertainty.m` is the same pipeline with a switch `unc_model`: `'mult'` for `theta = thetabar (1 + rho delta)` and `'log'` for `theta = thetabar exp(rho delta)`, which keeps every parameter positive for any `rho`. The first-order blocks are identical for the two maps; only the curvature differs, and the logarithmic map adds a diagonal correction to the Hessian bound (Section 5.9).

### 5. Physicality-aware robustness radius — `physicality_radius/empirical_radius_l2.m`
*Report: Section 5.7 and Table in 5.7.1.*

Checks the certified radii against brute force and adds a physicality constraint: a certificate is algebraic and does not know that parameters, equilibrium coordinates and Jacobian entries must stay positive. Pick one of six experiments with `experiment`:

`A_total`, `A_eqparam`, `A_bdc`, `B_total`, `B_eqparam`, `B_bdc` (Model × method; "total" = Method 4).

For each one the script sweeps a dilation α of the box (all 2ᵐ corners, uniform interior draws with local refinement of any unstable hit) and a sphere sweep for the 2-norm. It then reports two physicality walls:

- **closed-form wall** `α_phys = min_i(mid_i / hw_i)` for the BDC and equilibrium-as-parameter methods, where each quantity is built as `mid + hw·delta`;
- **empirical wall**, the first α on the grid at which a corner becomes non-physical (the only option for the total-derivative method, whose equilibrium is re-solved).

The final answer is `min(1/µ, α_phys)`, and the script says whether the certificate or physicality is the binding constraint.

The certified values `muInv` and `gamInv` for each experiment are hard-coded in `parse_experiment`, and the equilibrium boxes in the model constants; update them if the upstream scripts change.

