# Changelog

The version of record is the published paper (BioSystems 258, 105608; DOI in the
README), tagged `v1.0.0`. This repository tracks ongoing development; `main` carries
the maintained version.

## 2026-07-12 — corrections

- **Landauer framed as a detection-energy floor.** v1.0 stated that measuring a
  sub-Landauer phenomenon "would inject more energy than the phenomenon contains,
  fundamentally altering the system under study" (and similar phrasings). This conflates
  two quantities: Landauer's *k*<sub>B</sub>*T* ln 2 bounds the dissipation of an
  irreversible *erasure*, not the energy a *detector* must inject to read a signal — and the
  paper's own stochastic-resonance pooling mechanism shows sub-Landauer signals are
  collectively detectable without per-unit energy injection exceeding the signal. v2.0
  removes the detection-energy-floor phrasing, states plainly that "Landauer's principle
  concerns the thermodynamic cost of erasing information, not detecting it," reframes the
  limit as a *reference scale* for reliable single-event discrimination against thermal noise
  (the substantive constraint being reliability, not energy injection), and hedges the
  single-channel gating figure as a "work-scale proxy." The qualitative claim — a
  sub-Landauer regime of causally potent, not-reliably-binary-readable patterns — is
  unchanged.

## 2026-06-11 — corrections

- **Stochastic-resonance scaling.** v1.0 reported pooled signal-to-noise figures
  (N=100 → ~30, N=500 → ~135) and described them as confirming SNR ∝ √N. Those
  numbers are a *power* SNR (∝ N; note 135/30 = 4.5 = 500/100), whereas the analytic
  prediction is the *amplitude* SNR (∝ √N). v2.0 reports the amplitude form. The
  qualitative result (ensemble pooling exposes sub-threshold signals) is unchanged.

## v2.0 (December 2025)
- Added the Framework-Dependence section; extended the Duhem–Quine discussion;
  reports the amplitude-SNR scaling (see correction above).

## v1.0 (October 2025)
- Published in BioSystems (tag `v1.0.0`).
