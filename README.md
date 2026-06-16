<div align="center">

<a href="https://github.com/effjy/krakken-cryptanalysis/"><img src="titles/krakken-2048-cryptanalysis-title.svg" height="52" alt="Krakken-2048 Cryptanalysis"></a>

Krakken-2048 — Cryptanalysis & Verification Suite

**A 2048-bit SPN-ARX wide-trail permutation, its sponge hash, and the empirical + provable analysis behind it.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](#license)
[![Language: C](https://img.shields.io/badge/language-C-00599C.svg?logo=c)](#)
[![SIMD: AVX2](https://img.shields.io/badge/SIMD-AVX2-orange.svg)](#)
[![Build: make](https://img.shields.io/badge/build-make-success.svg)](#building--running-the-tests)
[![Permutation: 8 rounds](https://img.shields.io/badge/rounds-8-blue.svg)](#)
[![State: 2048-bit](https://img.shields.io/badge/state-2048--bit-9cf.svg)](#)
[![Active S-boxes: 229](https://img.shields.io/badge/min%20active%20S--boxes-229-brightgreen.svg)](#what-the-papers-cover)
[![Char. bound: 2⁻¹³⁷⁴](https://img.shields.io/badge/characteristic%20bound-2%E2%81%BB%C2%B9%C2%B3%E2%81%B7%E2%81%B4-brightgreen.svg)](#what-the-papers-cover)
[![NIST SP 800-22: PASS](https://img.shields.io/badge/NIST%20SP%20800--22-PASS-success.svg)](#3-nist-sp-800-22-on-the-sponge-keystream)
[![Division property: MILP](https://img.shields.io/badge/division%20property-MILP%20(SCIP)-blueviolet.svg)](#4-division-property-integral-search--division)
[![Papers: 2](https://img.shields.io/badge/papers-2-blueviolet.svg)](#-download-the-papers)

</div>

---

## Overview

**Krakken-2048** is a 2048-bit cryptographic permutation built in the wide-trail
tradition (AES / Whirlpool / Keccak lineage). Its state is 32 × 64-bit words
viewed as an 8 × 4 matrix, and one round applies nine layers in order:

| # | Layer | Role |
|---|-------|------|
| 1 | **Theta** | column-parity mixing |
| 2 | **Tentacle (MDS)** | circulant GF(2⁸) Whirlpool matrix, branch number **9** |
| 3 | **Rho** | per-lane bit rotations |
| 4 | **Pi** | lane permutation |
| 5 | **Chi** | nonlinear layer, ABYSSAL 8-bit S-box (diff. uniformity 4) |
| 6 | **XRBD** | XOR-Rotation Butterfly Diffusion — the headline diffusion primitive |
| 7 | **Pressure** | 64-bit ARX intra-column step |
| 8 | **Iota** | SHAKE128-derived round constants |
| 9 | **InkCloud** | lane shuffle |

The permutation iterates **8 rounds**. It is used in a sponge construction
(rate 160 bytes, capacity 96 bytes) to produce the **Krakken-2048 hash**, and as
a keystream generator in a stream-cipher / volume-encryption mode.

This repository collects the **analysis** of that design: two companion papers
and the **reproducible test harnesses** that produced the empirical results in
them.

---

## 📄 Download the papers

Two companion notes accompany the main (separately published) construction paper.
Click to download the LaTeX source; compiled PDFs are also included in this repo.

> ### ⬇️ Paper 1 — Experimental Verification & Distinguisher Battery
> **[`krakken_experiments.tex`](./krakken_experiments.tex)** &nbsp;·&nbsp; [PDF](./krakken_experiments.pdf)
>
> Exact proofs of the cryptographic foundations + a broad empirical distinguisher
> battery against the compiled permutation, hash, and keystream.

> ### ⬇️ Paper 2 — Differential Active-S-box Analysis (Findings)
> **[`krakken_findings.tex`](./krakken_findings.tex)** &nbsp;·&nbsp; [PDF](./krakken_findings.pdf)
>
> The MILP active-S-box bound analysis (differential & linear), trail
> verification, the decisive solver finding, and the reduced-round + hash-level
> experiments.

To rebuild a PDF locally:

```bash
pdflatex krakken_experiments.tex
pdflatex krakken_findings.tex
```

---

## What the papers cover

### Provable bounds (MILP)

- **229 active S-boxes** proven over the full 8-round permutation, for **both
  differential and linear** trails, via an exact byte-lane MILP solved to a
  **closed optimality gap** (SCIP) — every per-round value `5, 37, 69, 101, 133,
  165, 197, 229` is proven directly, not extrapolated.
- This bounds any single differential/linear characteristic by **2⁻¹³⁷⁴**.
- **Isolation of the XRBD layer:** the SPN core alone admits a cheap lane-aligned
  trail (27 active S-boxes at 2 rounds); XRBD provably destroys it (→ 37 at 2
  rounds, 69 at 3).
- **Methodological finding:** SCIP closes these instances where CBC's lower
  bound stays frozen — a strong cut-generating solver is the decisive ingredient.

### Exhaustively verified foundations

- **MDS branch number 9** — all 12,869 square submatrices nonsingular over GF(2⁸).
- **ABYSSAL S-box** — bijection, differential uniformity 4 (DP = 2⁻⁶), max linear
  bias 16 (corr² = 2⁻⁶), BCT max 6, via full DDT / LAT / BCT enumeration.
- **Permutation layers** (Pi, InkCloud) bijective; verified inverse permutation
  round-trips at all round counts.

### Empirical distinguisher battery (no distinguisher survives 8 rounds)

Avalanche, rotational / RX-difference, integral / saturation, cube / algebraic
degree, exact algebraic degree (Möbius), forward & two-sided zero-sum,
differential clustering, impossible differentials, boomerang / BCT,
differential-linear, and invariant-subspace / S-box-coset structure.

The **analytical** complement to the integral test — a sound **bit-based
division-property MILP** — is now implemented (Paper 1 §2.7, harness in
`division/`); see [§4 below](#4-division-property-integral-search--division).

### New test suites in this repository

Three additional families were added and written into both papers:

| Area | Target | Result |
|------|--------|--------|
| **SAC / completeness / monobit** | bare permutation | ideal at every round, full bit-level completeness |
| **NIST SP 800-22** | sponge keystream | all 15 test families pass |
| **Collision / distribution** | sponge hash | birthday-ideal, no near-collisions, uniform |

---

## Repository layout

```
.
├── krakken_experiments.tex / .pdf   # Paper 1 — experimental verification & battery
├── krakken_findings.tex   / .pdf    # Paper 2 — MILP active-S-box findings
├── diffusion_tests/                 # SAC / completeness / monobit harness
│   ├── sac_test.c
│   ├── krakken_multi.c              # AVX2 permutation (linked)
│   └── Makefile
├── collision_tests/                 # sponge-hash collision / distribution battery
│   ├── collision_test.c
│   ├── krakken_multi.c              # AVX2 permutation + hash (linked)
│   └── Makefile
└── division/                        # bit-based division-property integral search (MILP)
    ├── krakken.c / .h               # reference permutation (built to libkrakken.so)
    ├── gen_facets.py                # compact EXACT S-box division model generator
    ├── krakken_divprop.py           # the search (SCIP MILP)
    ├── model.py / cbdp.py / layers.py / sbox_trails.py
    ├── validate_*.py                # soundness gates
    └── Makefile / runs.sh / README.md
```

> **Note:** `krakken_multi.c` is the AVX2 reference implementation of the
> permutation and sponge hash; both harnesses link against it. It is included in
> each test directory so each builds standalone.

---

## The test harnesses

### 1. SAC / completeness / monobit — `diffusion_tests/`

Refines the *mean* avalanche figure into the full per-bit dependency structure.
For every (input bit *i*, output bit *j*) pair it measures

```
SAC[i][j] = Pr[ output bit j flips | input bit i is flipped ]   (ideal ½)
```

across the 2048 × 2048 = 4.19 M cell matrix, and derives:

- **SAC matrix** — max & mean deviation from ½, cells beyond 5σ.
- **Completeness** — every cell strictly between 0 and N (every input bit
  influences every output bit, none deterministically): **0 dead / 0 stuck cells
  at every round**.
- **Monobit** — per-output-bit balance over random inputs.

A correct *mean* avalanche can hide a structured dependency matrix; this test
catches that, and finds nothing — the matrix is statistically ideal from round 1.

### 2. Sponge-hash collisions & distribution — `collision_tests/`

Three quick probes of the 256-bit sponge hash:

- **Birthday collisions** — hash N messages, truncate digests to *t* bits, count
  colliding pairs against the ideal `N²/2^(t+1)`. Observed tracks expected; **no
  collision excess** at any truncation width.
- **Near-collision diffusion** — 1-bit input difference → 256-bit digest
  difference: mean Hamming weight **127.98** (ideal 128), min **93** (far from 0).
- **Digest uniformity** — byte-value **χ² = 273.8** on 255 dof (z = +0.83).

### 3. NIST SP 800-22 on the sponge keystream

The permutation, used as a keystream generator
(`ciphertext = plaintext ⊕ Sponge(key ‖ nonce ‖ LE64(i))`) over **empty**
encrypted volume containers (so the ciphertext is effectively the raw keystream),
was assessed with the NIST Statistical Test Suite over 10 × 10⁶-bit sequences.
**All 15 test families pass** for all tested files (Frequency, Block Frequency,
Cumulative Sums, Runs, Longest Run, Rank, FFT, Non-Overlapping Template,
Overlapping Template, Universal, Approximate Entropy, Random Excursions &
Variant, Serial, Linear Complexity). A high-throughput monobit/runs cross-check
over 10⁸ bits is unremarkable, with byte entropy > 7.999 bits/byte.

### 4. Division-property integral search — `division/`

The analytical counterpart to the empirical integral test: the **bit-based
division property** (two-subset CBDP), modelled as a MILP and solved with **SCIP**.
It propagates a division trail through the round function and proves an output bit
**balanced** (zero-sum over the chosen input cube) exactly when its trail is
infeasible — a *proof*, not a statistical observation. The model is **sound**:
every bit it certifies balanced is a genuine zero-sum bit (the ARX and S-box
models are conservative, so it may miss some but never reports a false one).

The practical obstacle for a 2048-bit state is the S-box: the textbook *selector*
model adds ≈ 2015 binary variables per S-box (≈ 5 × 10⁵ per round), and SCIP could
not resolve a single output bit. This harness replaces it with a compact **exact**
inequality model (`gen_facets.py`: **661 inequalities, no auxiliary variables**,
verified over all 2¹⁶ points of the S-box division polytope), which drops the
one-round program from **937,952 → 422,112 variables** and makes it solvable —
round 1 now solves to optimality where the selector model never terminated. The
exhaustive multi-round sweep that would fix the integral-distinguisher length is a
long-running computation; this is the harness behind §2.7 of Paper 1.

```bash
cd division
make            # build libkrakken.so + generate the exact facet model (~10 s)
./runs.sh 1     # round 1, full word:0 cube, all output bits (logs to runs/)
./runs.sh validate   # correctness gates (wiring, ARX rule, S-box)
```

> **Requirements differ from the C harnesses:** a Python venv with
> **pyscipopt** (SCIP), **numpy** and **scipy** (default `~/venv`; override with
> `PY=/path/to/python`), plus `gcc`. See `division/README.md` for full usage.

---

## Building & running the tests

**Requirements:** a C compiler (`gcc`/`clang`), an **AVX2-capable** CPU,
`make`, and a math + pthread toolchain (standard on Linux).

```bash
# SAC / completeness / monobit
cd diffusion_tests
make            # build
make quick      # fast smoke test  (rounds=8, N=512)
make run        # default probe     (rounds=8, N=4096)
make scan       # all rounds 1..8   (~3 min)

# Sponge-hash collision / distribution
cd ../collision_tests
make            # build
make quick      # fast smoke test  (N=100,000)
make run        # full battery      (N=1,000,000, ~3 s)
```

Both harnesses use a deterministic `splitmix64` PRNG (seeded by the round count)
so every run is reproducible.

---

## Scope & honest limitations

These results are offered as a **thorough self-conducted first-pass evaluation**,
not as a substitute for external cryptanalysis:

- The active-S-box bounds are **proven**; the distinguisher, SAC, NIST, and
  collision results are **empirical** — each rules out a *detectable* weakness of
  its kind at the sampled scale, which is necessary but not sufficient.
- The proven bounds concern **single characteristics**, not differential clusters
  or linear hulls.
- The linear bound coincides with the differential one by wide-trail symmetry,
  not via an independent bit-level mask model.
- The ARX layer is modelled conservatively (carries ignored), keeping the counts
  sound lower bounds.
- **Independent external cryptanalysis remains the essential next step**, and
  full reproducible artifacts are provided to support it.

---

## Citation

```bibtex
@misc{krakken2048,
  author = {Lachance-Caumartin, Jean-Fran\c{c}ois},
  title  = {Krakken-2048: Provable Active-S-box Bounds, Experimental
            Verification, and Distinguisher Battery},
  year   = {2026},
  note   = {ORCID: 0009-0005-6377-1675}
}
```

---

## License

Released under the **MIT License**. See `LICENSE` (add one before publishing if
not present).

---

<div align="center">
<sub>Built for scrutiny — independent cryptanalysis welcome.</sub>
</div>
