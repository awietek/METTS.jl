using METTS
using Dumper
using HDF5
using Printf

# Usage:
#   julia chained_bar_analysis.jl metts_bc2.0.h5 metts_bc4.0.h5 metts_bc8.0.h5 ...
#
# Each file must have been produced by heisenberg_interval.jl and contain:
#   beta_collapse, beta_min, beta_max, tau  (scalars, written with write_data!)
#   betas       (1 x n_betas, written once with dump!)
#   energy      (nmetts x n_betas, one row per METTS sample)
#   log_norm    (nmetts x n_betas, one row per METTS sample)

"""
Reads one METTS HDF5 output file and returns the data needed for chained BAR.

Returns
-------
beta_collapse  : Float64
betas          : Vector{Float64} of measured beta values
measurements   : Vector of Dicts, one per METTS sample; each Dict maps
                 observable name -> Vector{Float64} over betas
"""
function load_metts_file(filename::AbstractString)
    dfile = DumpFile(filename)

    beta_collapse = read_data(dfile, "beta_collapse")
    betas         = vec(read_data(dfile, "betas"))   # (1, n_betas) -> n_betas

    energy_mat   = read_data(dfile, "energy")        # (nmetts, n_betas)
    log_norm_mat = read_data(dfile, "log_norm")      # (nmetts, n_betas)
    nmetts = size(energy_mat, 1)

    measurements = [Dict("energy"   => energy_mat[i, :],
                         "log_norm" => log_norm_mat[i, :])
                    for i in 1:nmetts]

    return beta_collapse, betas, measurements
end

# ---- Load all files -------------------------------------------------------

filenames = length(ARGS) > 0 ? ARGS : error("usage: julia chained_bar_analysis.jl file1.h5 file2.h5 ...")

data = [load_metts_file(f) for f in filenames]

# Sort by beta_collapse ascending
sort!(data; by = d -> d[1])

beta_collapses   = [d[1] for d in data]
all_betas        = [d[2] for d in data]
all_measurements = [d[3] for d in data]

println("Loaded $(length(filenames)) files:")
for (k, (bc, betas, meas)) in enumerate(zip(beta_collapses, all_betas, all_measurements))
    println("  state $k: beta_collapse = $bc, window = [$(betas[1]), $(betas[end])], N = $(length(meas))")
end

# ---- MBAR free energies and reweighted observables -----------------------

n_bootstrap = 20
betas_out, energy_mean, energy_err, energy_neff, free_energies, free_energy_err =
    bootstrap_mbar_reweight_observable(all_measurements, all_betas, beta_collapses,
                                       "energy"; n_bootstrap=n_bootstrap)

println("\nMBAR free energies (bootstrap n=$n_bootstrap):")
println("  beta_collapse    F            stderr")
for (k, bc) in enumerate(beta_collapses)
    @printf("  %10.4f  %12.6f  %12.6f\n", bc, free_energies[k], free_energy_err[k])
end

println("\nMBAR energy estimates (bootstrap n=$n_bootstrap):")
println("  beta       energy        stderr        N_eff")
for (beta, e, err, neff) in zip(betas_out, energy_mean, energy_err, energy_neff)
    @printf("  %6.3f  %12.6f  %12.6f  %8.1f\n", beta, e, err, neff)
end

# Save results
out = DumpFile("chained_bar_results.h5")
out["betas"]            = betas_out
out["energy_mean"]      = energy_mean
out["energy_err"]       = energy_err
out["energy_neff"]      = energy_neff
out["free_energies"]    = free_energies
out["free_energy_err"]  = free_energy_err
out["beta_collapses"]   = beta_collapses
println("\nResults written to chained_bar_results.h5")
