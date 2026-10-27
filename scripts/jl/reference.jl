# Rodas5P reference trace of a Julia model in the bimgrl trace contract.
#
#   julia --startup-file=no --project=<skill-dir>/scripts/jl reference.jl <model.jl> <out_stem> --t1 1000 [--bcl 1000 --beats 10]
#       [--stim-start 0 --stim-duration 1] [--reltol 1e-8 --abstol 1e-10] [--grid 0.1]
#
# The model file defines STATES, U0, P and rhs!(du,u,p,t) as for adapter.jl. If --bcl is given,
# the run is protocol B: `beats` beats, the state at the start of every beat is saved to
# <out_stem>_beats.csv and the last beat on the grid to <out_stem>.csv; the stimulus edges are
# tstops so the solver never steps across them. Without --bcl: protocol A, 0..t1.

using OrdinaryDiffEq
using DelimitedFiles
using JSON
using Dates

function getarg(name, default)
    i = findfirst(==(name), ARGS)
    i === nothing ? default : parse(typeof(default), ARGS[i + 1])
end

const MODEL = ARGS[1]
const OUT = ARGS[2]
include(abspath(MODEL))
t1 = getarg("--t1", 1000.0)
bcl = getarg("--bcl", 0.0)
beats = getarg("--beats", 10)
stim_start = getarg("--stim-start", 0.0)
stim_dur = getarg("--stim-duration", 1.0)
reltol = getarg("--reltol", 1e-8)
abstol = getarg("--abstol", 1e-10)
grid = getarg("--grid", 0.1)

f!(du, u, p, t) = rhs!(du, u, p, t)
if bcl > 0
    t1 = beats * bcl
end
tstops = bcl > 0 ? sort(vcat([stim_start + k * bcl for k in 0:beats-1], [stim_start + stim_dur + k * bcl for k in 0:beats-1])) : [stim_start, stim_start + stim_dur]
prob = ODEProblem(f!, Float64.(U0), (0.0, t1), P)
t0w = time()
sol = solve(prob, Rodas5P(); reltol=reltol, abstol=abstol, tstops=tstops, maxiters=10^8)
wall = time() - t0w
start = bcl > 0 ? t1 - bcl : 0.0
ts = collect(start:grid:t1)
Y = reduce(vcat, (sol(t)' for t in ts))
open(OUT * ".csv", "w") do io
    println(io, join(["t"; STATES], ","))
    writedlm(io, hcat(ts, Y), ',')
end
if bcl > 0
    tb = [k * bcl for k in 0:beats]
    Yb = reduce(vcat, (sol(t)' for t in tb))
    open(OUT * "_beats.csv", "w") do io
        println(io, join(["t"; STATES], ","))
        writedlm(io, hcat(tb, Yb), ',')
    end
end
meta = Dict("model" => basename(MODEL), "scheme" => "reference", "dt" => nothing, "reference" => true,
            "solver" => Dict("name" => "Rodas5P", "rtol" => reltol, "atol" => abstol), "language" => "julia",
            "t_unit" => "ms", "protocol" => bcl > 0 ? "B" : "A", "beats" => bcl > 0 ? beats : 1, "bcl" => bcl > 0 ? bcl : nothing,
            "stimulus" => Dict("start" => stim_start, "duration" => stim_dur, "period" => bcl, "amplitude" => nothing),
            "status" => "ok", "message" => "", "n_steps" => length(sol.t), "wall_s_total" => wall, "wall_s_per_beat" => nothing,
            "parameters" => Dict(), "float_type" => "float64", "blocks" => nothing, "created" => string(now()))
open(OUT * ".json", "w") do io
    JSON.print(io, meta, 2)
end
println("wrote $(OUT).csv ($(length(ts)) rows, $(length(sol.t)) accepted steps, $(round(wall; digits=1)) s)")
