# METTS.jl — Quick Guide & Example (Fermi–Hubbard on a Cylinder)

This repository contains a minimal, reproducible example for running
**METTS** (Minimally Entangled Typical Thermal States) with **ITensors**
on the 2D Fermi–Hubbard model wrapped as a cylinder. The example
computes finite‑temperature observables and writes results to HDF5.

---

## Contents

```
examples/hubbard_cylinder_metts.jl   <-- runnable example script
src/                                <-- package code
outfiles.metts/                     <-- created at runtime; outputs written here
```

---

## Requirements

- Julia **1.11** or newer
- Packages installed into **this** project: `ITensors`, `ITensorMPS`, `HDF5`

> **Note:** Do **not** add `METTS` from the registry when working inside this repo.
> Activating this project makes `using METTS` work.

---

## Quick Start

```bash
# 1) Enter the project
cd /path/to/METTS.jl

# 2) Instantiate and add required packages *into this project*
julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.add(["ITensors","ITensorMPS","HDF5"])'

# 3) Run an example (L W U filling t T tau maxD seed nmetts Nwarm)
julia --project examples/hubbard_cylinder_metts.jl   24 4 10.0 0.9375 1.0 0.10 0.3 2000 1 100 5

# 4) Find produced HDF5 outputs
find outfiles.metts -maxdepth 8 -name outfile.h5
```

---

## Arguments

| Argument | Meaning |
|---|---|
| `L` | Lattice length (total sites = `L × W`) |
| `W` | Lattice width (cylinder) |
| `U` | On‑site interaction |
| `filling` | Electrons per site (e.g. `0.9375`) |
| `t` | Nearest‑neighbor hopping (usually `1.0`) |
| `T` | Temperature (`β = 1/T`) |
| `tau` | TDVP step for imaginary‑time evolution |
| `maxD` | Maximum MPS bond dimension |
| `seed` | RNG seed |
| `nmetts` | Total METTS iterations |
| `Nwarm` | Warm‑up iterations (not stored) |

**Notes**

- Measurements are written only *after* warm‑up.
- Ensure `nmetts > Nwarm`.

---

## Outputs

**Default output root:** `./outfiles.metts` (or the directory set by `METTS_OUTDIR`).

Each run is saved under:
```
outfiles.metts/L.<L>.W.<W>/U.<U>/filling.<f>/T.<T>/D.<maxD>/tau.<tau>.seed.<seed>/
```

**Within a run directory you will find:**
- **`outfile.h5`** — parameters under `meta/`; groups `step_xxxxx/` for each measured step  
  (datasets include: `energy`, `Ntotal`, `DO`, `SZSZ`, `SPSM`, `SMSP`, `NtotNtot`, `SS`).
- **`samples.txt`** — product‑basis states per measured step. Codes `1..4` map to `Emp`, `Up`, `Dn`, `UpDn`.

Resuming a run is automatic if `samples.txt` already exists.

---

## Minimal Runs

**Sanity check:**
```bash
julia --project examples/hubbard_cylinder_metts.jl   4 4 10.0 0.9375 1.0 0.10 0.3 100 1 10 5
```

**Paper‑grade template:**
```bash
julia --project examples/hubbard_cylinder_metts.jl   32 4 10.0 0.9375 1.0 0.10 0.2 2000 1 2000 100
```

---

## Troubleshooting

- **ERROR: Package ITensors not found**  
  Use the project environment: run with `julia --project` (or `Pkg.activate(".")`).

- **No measured data appears**  
  Set `nmetts > Nwarm` so that measurements are taken after warm‑up.

- **Slow / memory‑heavy**  
  Decrease `maxD` or increase `tau` slightly.

---

## Performance Tips

- A short DMRG warm‑up avoids long initial transients.
- Choose `tau` and `maxD` carefully; typical `tau = 0.1–0.3`.