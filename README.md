# METTS.jl — Quick Guide & Example (extended \(t-t'-J\)  Model on a Cylinder)

This repository contains a minimal, reproducible example for running
**METTS** (Minimally Entangled Typical Thermal States) with **ITensors**
on the 2D extended **\(t-t'-J\) ** model wrapped as a cylinder. The script
computes finite-temperature observables, selected four-fermion correlators,
and writes results to HDF5.

---

## Contents

```
examples/ttJ_metts.jl           <-- runnable METTS script for the t-t'-J cylinder
src/                                <-- package code
outfiles.metts/                     <-- created at runtime; outputs written here
```

---

## Requirements

- Julia **1.11** or newer
- Packages installed into **this** project:
  - `ITensors`
  - `ITensorCorrelators`
  - `METTS`
  - `HDF5`
  - `Dumper`
  - `LinearAlgebra`
  - `Random`
  - `Serialization`
  - `Printf`

> **Note:** when working inside this repository, run Julia with `--project`
> so that `using METTS` and the other local dependencies resolve correctly.

---

## Model

The script simulates the extended **\(t-t'-J\)** Hamiltonian on a square-lattice cylinder:

- open boundary conditions along the long direction \(x\)
- periodic boundary conditions along the wrapped direction \(y\)

The Hamiltonian contains:

- nearest-neighbor hopping with amplitude `t`
- next-nearest-neighbor diagonal hopping with amplitude `t_prime`
- nearest-neighbor spin exchange `J`

The local Hilbert space is the standard no-double-occupancy **tJ** basis:
`Emp`, `Up`, `Dn`.

---

## Quick Start

```bash
# 1) Enter the project
cd /path/to/METTS.jl

# 2) Instantiate and add required packages into this project
julia --project -e 'using Pkg; Pkg.instantiate()'

# 3) Run the t-t'-J METTS script
#    Arguments:
#    L W J filling t t_prime T tau maxD cutoff seed NMETTS Nwarm
julia --project ttJ_metts.jl \
  24 4 0.4 0.9375 1.0 0.2 0.10 0.2 2000 1e-6 1 100 5

# 4) Locate the produced output
find /metts.cylinder.tj.tjp/outfiles.metts -name outfile.h5
```

---

## Arguments

| Argument | Meaning |
|---|---|
| `L` | Lattice length |
| `W` | Lattice width (cylinder circumference) |
| `J` | Spin-exchange coupling |
| `filling` | Electron filling per site, e.g. `0.9375` |
| `t` | Nearest-neighbor hopping |
| `t_prime` | Next-nearest-neighbor diagonal hopping |
| `T` | Temperature, with `β = 1/T` |
| `tau` | Imaginary-time TDVP step |
| `maxD` | Maximum MPS bond dimension |
| `cutoff` | Truncation cutoff |
| `seed` | Random seed |
| `NMETTS` | Total METTS iterations including warm-up |
| `Nwarm` | Number of warm-up iterations not stored |

### Usage

```bash
julia --project ttJ_metts.jl \
  L W J filling t t_prime T tau maxD cutoff seed NMETTS Nwarm
```

## What the Script Does

For a given parameter set, the script:

1. builds the \(t-t'-J\)  Hamiltonian as an MPO on a cylinder
2. initializes a product state at the target filling
3. optionally performs a **DMRG warm-up** if the run starts from scratch
4. runs the METTS loop using imaginary-time evolution via TDVP
5. measures local observables and correlation functions after warm-up
6. periodically measures selected four-fermion correlators
7. stores collapse samples so the run can be resumed automatically

---


### Fresh run

If no previous `samples.txt` exists, the script:

- builds a random product state at the requested filling
- enforces equal numbers of up and down electrons
- performs a DMRG warm-up before starting the METTS sampling

The DMRG warm-up depends on temperature:

- **low \(T\)**: longer warm-up and larger bond dimensions
- **intermediate \(T\)**: moderate warm-up
- **high \(T\)**: shorter warm-up

### Restarted run

If a previous `samples.txt` is found, the script:

- reads the last stored collapse sample
- reconstructs the product state from it
- skips DMRG warm-up
- resumes from the next total METTS step automatically

So restart is automatic as long as the run directory already contains `samples.txt`.

---

## Outputs

### Output root

This script currently writes to the fixed path:

```text
/metts.cylinder.tj.tjp/outfiles.metts
```

Each run is saved under:

```text
L.<L>.W.<W>/J.<J>/t.<t>/t_prime.<t_prime>/filling.<filling>/T.<T>/D.<maxD>/tau.<tau>.cutoff.<cutoff>.seed.<seed>/
```

Inside that directory you will typically find:

```text
outfile.h5
samples.txt
```

---

## HDF5 Contents

The main measured observables are written to `outfile.h5`.

For each measured METTS sample, the script dumps:

- `energy` — total energy
- `SZ` — local \( \langle S_i^z \rangle \)
- `Ntotal` — local density \( \langle n_i \rangle \)
- `SZSZ` — \( \langle S_i^z S_j^z \rangle \)
- `SPSM` — \( \langle S_i^+ S_j^- \rangle \)
- `SMSP` — \( \langle S_i^- S_j^+ \rangle \)
- `NtotNtot` — \( \langle n_i n_j \rangle \)
- `SS` — full spin correlator built as  
  \[
  \mathbf S_i \cdot \mathbf S_j
  =
  \langle S_i^z S_j^z \rangle
  + \tfrac12 \left(
  \langle S_i^- S_j^+ \rangle
  + \langle S_i^+ S_j^- \rangle
  \right)
  \]
- `cdagupcup` — one-body spin-up fermion correlator
- `cdagdncdn` — one-body spin-down fermion correlator
- `svn` — bipartite von Neumann entanglement entropy at the center cut

### Four-fermion correlators

The script also computes selected four-fermion correlators of the form

\[
\langle
c^\dagger_{\uparrow}
c^\dagger_{\downarrow}
c_{\uparrow}
c_{\downarrow}
\rangle
\]

for combinations of nearest-neighbor bond indices, and stores them under:

```text
fourfermioncorr/CdagupCdagdnCupCdn/step_<meas>/idx_i_j_k_l
```

These are scheduled less frequently depending on temperature:

- `T <= 0.04`: every measured step
- `0.04 < T <= 0.08`: every 2 measured steps
- `0.08 < T <= 0.2`: every 3 measured steps
- `T > 0.2`: every 5 measured steps

This keeps the heavy correlator measurements under control.

---

## Collapse Samples

The file `samples.txt` stores the collapsed product-basis sample after each measured METTS step.

Each line has the form:

```text
meas_index: [ ... integers ... ]
```

The integer codes are interpreted as:

- `1` → `Emp`
- `2` → `Up`
- `3` → `Dn`

These samples are used for restart.

---

## Minimal Runs

### Sanity check

```bash
julia --project ttJ_metts.jl \
  4 4 0.4 0.875 1.0 0.2 0.20 0.2 200 1e-6 1 10 5
```

### Typical production-style run

```bash
julia --project ttJ_metts.jl \
  24 4 0.4 0.9375 1.0 0.2 0.05 0.2 2500 1e-6 1 2000 10
```


---

## Performance Notes

A few important implementation details:

- imaginary-time evolution uses `timeevo_tdvp_extend`
- the script evolves by `-β/2` at each METTS step
- the post-evolution collapse basis is `"X"`
- the post-DMRG initial collapse basis is `"Z"`


---

## Practical Notes

- The total number of electrons is set by `round(Int, filling * N)`.
- The script requires this number to be even, since it initializes equal numbers of up and down electrons.
- Measurements are only stored for steps with `step_total > Nwarm`.
- Make sure `NMETTS > Nwarm`, otherwise you will get no stored measurements.

---

## Troubleshooting

### `Ne must be even`
The chosen `filling × (L*W)` produced an odd number of electrons.
Choose parameters so that the total particle number is even.

### No measured data appears
Check that:

```text
NMETTS > Nwarm
```

Otherwise the run never enters the measurement stage.

### Restart behaves strangely
Delete the previous `samples.txt` if you want to force a fully fresh run.

### Very slow execution
The most expensive parts are:

- TDVP imaginary-time evolution
- full correlation matrices
- four-fermion correlators

For testing, reduce:

- `L`, `W`
- `maxD`
- `NMETTS`

or increase `T`.

---

## Performance Tips

- A DMRG warm-up helps reduce initial transient behavior.
- Four-fermion correlators are expensive; the script already throttles them by temperature.
- For quick checks, start with small cylinders and small `NMETTS`.
- Keep `tau` moderate; values that are too small can be slow, while values that are too large can hurt accuracy.

---

## Summary

This script is a METTS driver for the **extended \(t-t'-J\) model on cylindrical lattices**. It is meant for finite-temperature studies of:

- local density and spin structure
- spin and density correlation functions
- one-body fermionic correlators
- entanglement entropy
- selected four-fermion observables relevant to pairing physics

with automatic restart through stored collapse samples.