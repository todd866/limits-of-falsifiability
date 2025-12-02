# The Limits of Falsifiability: Simulation Suite

**A computational companion to "The Limits of Falsifiability" (BioSystems 258, 2025)**

**Author:** Ian Todd, University of Sydney
**DOI:** [10.1016/j.biosystems.2025.105608](https://doi.org/10.1016/j.biosystems.2025.105608)
**Timeline:** Received 18 Aug 2025 → Accepted 2 Oct 2025 → Published 16 Oct 2025

This repository demonstrates the physical and mathematical constraints on Popperian falsifiability in biological systems. It provides **executable proofs** that binary hypothesis testing fails when systems operate in high-dimensional phase spaces near thermodynamic limits.

---

## The Core Argument: System vs Shadow

The central thesis of this work is **ontological**, not merely statistical. We distinguish between:

| Concept | Definition | Example |
|---------|------------|---------|
| **The System (D_sys)** | The high-dimensional causal manifold where biological function actually occurs | The protein folding landscape, the neural state space |
| **The Shadow (D_obs)** | The low-dimensional projection accessible to measurement | A time series, a binary assay, electrode recordings |

Standard falsifiability assumes we can project D_sys → D_obs via a binary cut without destroying causal structure. These simulations demonstrate that when **D_sys >> D_obs** and energy scales are low (Sub-Landauer), this assumption is **mathematically incoherent**.

---

## Papers

| File | Description |
|------|-------------|
| `biosystems_2025_published.pdf` | **Paper 1:** The original paper as published in BioSystems 258 (October 2025) |
| `paper2_shadow_geometry.pdf` | **Paper 2:** "The Geometry of Biological Shadows" — computational companion quantifying topological aliasing |

**Paper 1** (The Warning) was drafted with Claude 4 (web) and establishes the philosophical argument.

**Paper 2** (The Ruler) was developed with Claude 4.5 Opus, GPT-5.1 Pro, and Gemini 3 Pro. It provides:
- Quantitative metrics for topological aliasing (47% misclassification, 199 "teleportations")
- The √N scaling law for sub-Landauer detection
- The falsifiability regime diagram

Paper 2 incorporates feedback from the Fulcher Lab at USyd ("define dimensionality") and operationalizes the D_sys/D_obs distinction into measurable quantities.

---

## Modules & Philosophical Claims

### 1. The Shadow Box (`05_shadow_box.py`) ⭐ FLAGSHIP

**"Plato's Cave for Falsifiability"**

- **The Concept:** Demonstrates the ontological gap between a chaotic system (Lorenz attractor, D_sys ≈ 2.06) and its observation (D_obs = 2).
- **The Proof:** Shows that a "clean" binary falsification line in the observation space groups together topologically disconnected causal states.
- **Key Metrics:**
  - **Aliasing Rate:** Percentage of states where shadow truth contradicts system truth (~47%)
  - **Topological Violations:** Times the shadow "teleports" while the system flows continuously

### 2. Binary Projection Problem (`01_binary_projection.py`)

**"The 1/k^n Problem" — Eq. 1**

- **The Concept:** Models the information destruction under binary projection.
- **The Proof:** As system dimensionality (n) rises, information preserved by any single binary test drops to zero (1/k^n).
- **Result:** Binary tests approach chance accuracy (50%) in high-D, even when perfectly distinguishable via multivariate methods.

### 3. Sub-Landauer Dynamics (`02_sub_landauer_sr.py`)

**"Noise is Signal" — Eq. 8-9**

- **The Concept:** Models the Sub-Landauer Domain where signal energy E < k_B T ln 2.
- **The Proof:** Demonstrates **Stochastic Resonance**—a signal invisible to single sensors becomes clear when pooled across a population.
- **Key Equation:** SNR ∝ √N. Falsifiability is an emergent property of the ensemble, not the unit.

### 4. Predictability Horizon (`03_predictability_horizon.py`)

**"The Hard Limit on Prediction" — Eq. 3-4**

- **The Concept:** Models Lyapunov exponents in chaotic systems.
- **The Proof:** Measurement precision (Δx) only buys logarithmic time: T_pred ∝ ln(1/Δx).
- **Implication:** Even with quantum-limited sensors, specific biological trajectories are strictly unfalsifiable beyond T_pred.

### 5. Scale-Dependent Regimes (`04_scale_dependent.py`)

**"The Regime Diagram" — Principle 1**

- **The Concept:** Maps the boundary where Popperian logic breaks down.
- **The Proof:** Generates a phase diagram of Signal Strength × Dimensionality.
- **Result:** Defines the "Ensemble Regime" (High-D, Low-Signal) where we must switch from single-case falsification to pattern-matching epistemology.

### 6. Non-Ergodic Memory (`06_nonergodic_memory.py`) 🆕

**"The System Remembers in Dimensions We Cannot See"**

- **The Concept:** When hidden states carry memory, time averages ≠ ensemble averages.
- **The Proof:** Trajectories with different hidden states converge to DIFFERENT values (0.25 vs 0.75). The ensemble mean (0.5) is achieved by NO individual trajectory.
- **Key Insight:** "Just measure longer" doesn't help when the system is non-ergodic.

### 7. Sample Complexity (`07_sample_complexity.py`) 🆕

**"The Required Data Doesn't Exist"**

- **The Concept:** Coverage of high-D spaces requires exponentially many samples.
- **The Proof:** With 1,000 samples and 3 bins/dimension:
  - n=5: 99% coverage (fine)
  - n=10: 1.7% coverage (terrible)
  - n=15: 0.01% coverage (essentially zero)
- **Key Insight:** Rare events in high-D spaces are invisible—you will NEVER sample them.

---

## The Inference Trilemma

Paper 2 establishes that **all three classical escape routes from measurement uncertainty are blocked**:

| Escape Route | Why It Fails | Simulation |
|---|---|---|
| **Time averaging** | Non-ergodicity: hidden memory makes time avg ≠ ensemble avg | `06_nonergodic_memory.py` |
| **Ensemble averaging** | Curse of dimensionality: N ~ k^n samples required | `07_sample_complexity.py` |
| **Direct measurement** | Perturbation: energy injection destroys the phenomenon | (physical principle) |

> "This is not a technological limitation awaiting better instruments; it is a structural feature of high-dimensional systems operating near thermodynamic limits."

---

## Definitions

To address common critiques regarding dimensionality:

- **Nominal Dimensionality (N):** The number of variables measured (e.g., 100 neurons).
- **Intrinsic Dimensionality (D_int):** The degrees of freedom of the generating manifold.
  - *Note:* Falsifiability fails even when D_int is low (e.g., Lorenz attractor), provided the projection is orthogonal to causal flow.
- **Sub-Landauer Pattern:** A structure where E_pattern < k_B T ln 2, making it thermodynamically impossible to resolve as a discrete bit without ensemble averaging.
- **Participation Ratio:** D_PR = (Σλ_i)² / Σ(λ_i²) — operational measure of intrinsic dimensionality from covariance eigenvalues.

---

## Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Run the flagship simulation (The Shadow Box)
python sims/05_shadow_box.py

# Run the full suite
for script in sims/*.py; do python "$script"; done
```

All figures are saved to `figures/`.

---

## Key Equations

| Equation | Description |
|----------|-------------|
| Ω_preserved / Ω_total = 1/k^n | Information preserved under binary projection (Eq. 1) |
| E_Landauer = k_B T ln 2 ≈ 3.0 × 10⁻²¹ J | Landauer limit at 310 K (Eq. 2) |
| T_pred ≲ (1/λ) ln(L/Δx) | Predictability horizon for chaotic systems (Eq. 3-4) |
| SNR(Y) ∝ √N | Stochastic resonance scaling (Eq. 9) |

---

## The Argument in Brief

Karl Popper's falsifiability criterion assumes scientific hypotheses can be reduced to binary tests. We show this assumption is **scale-dependent** and saturates in high-dimensional biological systems operating near physical measurement limits.

Three constraints compound to limit falsifiability:

1. **High dimensionality** → binary projection destroys almost all information
2. **Thermodynamic limits** → the Landauer bound sets a floor on bit recording
3. **Chaotic dynamics** → predictability horizons limit deterministic specification

**Central thesis:** Popperian falsification is a special case applicable to low-dimensional systems with strong signals. High-dimensional biological systems require ensemble-based, multi-scale inference.

---

## Related Work

This is part of a research program on dimensional constraints in biology:

- **Todd (2025a):** "The limits of falsifiability" — Paper 1 (BioSystems 258)
- **Todd (2025b):** "The geometry of biological shadows" — Paper 2 (this repo, in prep)
- **Todd (2025c):** "Timing inaccessibility and the projection bound" (BioSystems)
- **Todd (2025d):** "The physics of immune cooperation" (submitted)

---

## License

Code: MIT License
Paper: © 2025 Ian Todd (Open Access CC-BY)
