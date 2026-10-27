# Jacobian sparsity of a Julia model by SparseConnectivityTracer (survives branches that
# Symbolics cannot trace), written as JSON for `bimgrl`'s structure analysis.
#
#   julia --startup-file=no --project=<skill-dir>/scripts/jl sparsity.jl <model.jl> <out.json>

using SparseConnectivityTracer
using JSON

const MODEL = ARGS[1]
include(abspath(MODEL))
f!(du, u) = rhs!(du, u, P, 0.0)
S = jacobian_sparsity(f!, zeros(length(U0)), Float64.(U0), TracerSparsityDetector())
rows, cols, _ = findnz(S)
deps = Dict{String,Vector{String}}()
for (r, c) in zip(rows, cols)
    push!(get!(deps, STATES[r], String[]), STATES[c])
end
open(ARGS[2], "w") do io
    JSON.print(io, Dict("states" => STATES, "deps" => deps, "nnz" => length(rows)), 2)
end
println("wrote $(ARGS[2]): $(length(rows)) structural entries")
