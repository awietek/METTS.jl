METTS.jl — Quick Guide and Example (Fermi–Hubbard on a cylinder)
================================================================

This repository contains a minimal, reproducible example for running
METTS (Minimally Entangled Typical Thermal States) with ITensors.jl
on the 2D Fermi–Hubbard model wrapped as a cylinder. The example
computes finite-temperature observables and writes results to HDF5.

Contents
--------
- examples/hubbard_cylinder_metts.jl  <-- runnable example script
- outfiles.metts/                     <-- outputs will be written here (created at runtime)
- src/         <-- package code and tests

Requirements
------------
- Julia 1.11 or newer
- Packages (installed into THIS project):
    ITensors
    ITensorMPS
    HDF5

Note: you do NOT add METTS from the registry here because you are already
inside the METTS repo. Activating this project makes `using METTS` work.

Quick Start (5 minutes)
-----------------------
1) Open a terminal and go to the repo root:
   cd /path/to/METTS.jl

2) Install dependencies into this project environment:
   julia --project -e 'using Pkg; Pkg.instantiate();
                       Pkg.add(["ITensors","ITensorMPS","HDF5"])'

3) Run the example (small test case):
   julia --project examples/hubbard_cylinder_metts.jl \
     24 4 10.0 0.9375 1.0 0.10 0.3 2000 1 100 5

   Arguments are positional:
     L W U filling t T tau maxD seed nmetts Nwarm

4) Check that outputs were written:
   find outfiles.metts -name outfile.h5 -maxdepth 8

   You should see a path like:
   outfiles.metts/L.24.W.4/U.10.000/filling.0.93750/T.0.10000/D.2000/tau.0.30.seed.1/outfile.h5
   The same folder will also contain samples.txt

What the example does
---------------------
- Builds the Fermi–Hubbard Hamiltonian on an L x W cylindrical square lattice
- Optional short DMRG warmup to get a reasonable initial state
- METTS iterations:
  - imaginary-time evolve by beta/2 using TDVP
  - measure standard observables
  - collapse to a product state in a chosen basis
  - repeat

Inputs (command-line arguments)
-------------------------------
L         Integer   lattice length along x (sites = L*W)
W         Integer   lattice width along y (cylinder)
U         Float64   on-site interaction
filling   Float64   electrons per site (e.g. 0.9375)
t         Float64   nearest-neighbor hopping (usually 1.0)
T         Float64   temperature (beta = 1/T)
tau       Float64   TDVP time step for imaginary time evolution
maxD      Int       max MPS bond dimension during time evolution
seed      Int       RNG seed
nmetts    Int       total number of METTS iterations
Nwarm     Int       number of warm-up iterations (no measurements stored)

Output files and layout
-----------------------
Base output directory is:
  - METTS_OUTDIR if the environment variable is set, else
  - ./outfiles.metts relative to the current working directory

Within that directory, runs are organized by parameters:
  outfiles.metts/L.<L>.W.<W>/U.<U>/filling.<filling>/T.<T>/D.<maxD>/\
  tau.<tau>.seed.<seed>/

Each run directory contains:
- outfile.h5
    meta/...
      Scalars recording L, W, U, filling, t, T, tau, maxD, seed, Nwarm
    step_000011/...
      energy         Float64
      Ntotal         Vector{Float64} length N
      DO             Vector{Float64} (double occupancy per site)
      SZSZ           Matrix{Float64} N x N
      SPSM           Matrix{Float64} N x N
      SMSP           Matrix{Float64} N x N
      NtotNtot       Matrix{Float64} N x N
      SS             Matrix{Float64} N x N

    A new group step_<index> is appended for every measured iteration
    (i.e. for steps > Nwarm).

- samples.txt
    Text file with one collapsed product-basis sample per measured step:
      <step>: [state codes]
    The state codes are integers 1..4 mapped in the script to
    "Emp", "Up", "Dn", "UpDn".

Resuming a run
--------------
The script auto-detects samples.txt. If present, it will:
- reconstruct the initial product state from the last saved sample
- skip the DMRG warmup
- continue METTS iterations, appending new step groups to outfile.h5
  and new rows to samples.txt

Minimal examples
----------------
Small sanity run:
  julia --project examples/hubbard_cylinder_metts.jl 4 4 10.0 0.9375 1.0 0.10 0.3 100 1 10 5

Bigger run template to get results in the paper (we reject about N warm samples initially; higher tempertuares require larger Nwarm; 
we also run 5 parallel seeds and average over them for each temperature):
  julia --project examples/hubbard_cylinder_metts.jl 32 4 10.0 0.9375 1.0 0.10 0.2 2000 1 2000 100

Notes:
- Measurements are stored only after warm-up. Ensure nmetts > Nwarm.

Troubleshooting
---------------
- ERROR: Package ITensors not found
  You are not in the project environment. Always run with:
    julia --project examples/...
  or add to the top of the script:
    import Pkg; Pkg.activate(joinpath(@__DIR__, "..")); Pkg.instantiate()
  and then run:
    julia examples/...

- No measured data appears
  Check that nmetts > Nwarm. Warm-up iterations do not write outputs.

- Slow or memory-heavy evolution
  Decrease maxD or increase tau slightly (with care)
  
Performance tips
----------------
- Start with a short DMRG warmup to avoid very long initial transients.
- Choose tau and maxD to balance accuracy and runtime. Typical tau is 0.1–0.3.
