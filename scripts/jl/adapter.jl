# JSON-lines RHS adapter for a Julia model, serving `bimgrl`'s native path.
#
#   julia --project=<skill-dir>/scripts/jl adapter.jl <model.jl> [--params params.json]
#
# <model.jl> is a Julia file that defines, after being `include`d:
#   const STATES :: Vector{String}        state names, in the order of u
#   const U0     :: Vector{Float64}       initial state
#   const P                                parameters (any Julia object; a NamedTuple, Vector or struct)
#   rhs!(du, u, p, t)                      the right-hand side, in place
#   optional: const UNITS :: Dict{String,String}; set_params!(p, dict) returning the new p
# Protocol (references/contracts.md): one JSON request per stdin line, one JSON reply per line.
# The Jacobian is exact (ForwardDiff); `info` reports jac = true.

using JSON
using ForwardDiff

const MODEL = ARGS[1]
include(abspath(MODEL))

p = P
const N = length(U0)

function rhs_batch(ts, ys)
    out = Vector{Vector{Float64}}(undef, length(ys))
    for (i, (t, y)) in enumerate(zip(ts, ys))
        du = similar(Float64.(y))
        rhs!(du, Float64.(y), p, Float64(t))
        out[i] = du
    end
    out
end

function jac_batch(ts, ys)
    out = Vector{Vector{Vector{Float64}}}(undef, length(ys))
    for (i, (t, y)) in enumerate(zip(ts, ys))
        f = (du, u) -> rhs!(du, u, p, Float64(t))
        J = ForwardDiff.jacobian(f, zeros(N), Float64.(y))
        out[i] = [J[r, :] for r in 1:N]
    end
    out
end

params_dict(p::NamedTuple) = Dict(string(k) => v for (k, v) in pairs(p) if v isa Number)
params_dict(p::AbstractVector) = Dict("p[$i]" => v for (i, v) in enumerate(p) if v isa Number)
params_dict(p) = Dict(string(f) => getfield(p, f) for f in fieldnames(typeof(p)) if getfield(p, f) isa Number)

for line in eachline(stdin)
    isempty(strip(line)) && continue
    req = JSON.parse(line)
    op = req["op"]
    reply = try
        if op == "info"
            Dict("n" => N, "states" => STATES, "y0" => U0, "parameters" => params_dict(p),
                 "units" => (@isdefined(UNITS) ? UNITS : nothing), "jac" => true)
        elseif op == "rhs"
            Dict("dy" => rhs_batch(req["t"], req["y"]))
        elseif op == "jac"
            Dict("J" => jac_batch(req["t"], req["y"]))
        elseif op == "set_params"
            global p = (@isdefined(set_params!) ? set_params!(p, req["parameters"]) : merge(p, (; (Symbol(k) => v for (k, v) in req["parameters"])...)))
            Dict("ok" => true)
        elseif op == "quit"
            println(JSON.json(Dict("ok" => true))); flush(stdout); break
        else
            Dict("error" => "unknown op $op")
        end
    catch e
        Dict("error" => sprint(showerror, e))
    end
    println(JSON.json(reply))
    flush(stdout)
end
