# A tiny model in the adapter's form, for testing: a two-state Markov pair plus a pool-buffer
# binding pair, no membrane potential. Conserves A + B (occupancy) and c + b (mass).
const STATES = ["A", "B", "c", "b"]
const U0 = [0.9, 0.1, 1e-4, 0.01]
const P = (alpha = 2.0, beta = 0.5, kon = 100.0, koff = 0.1, bmax = 0.05)
function rhs!(du, u, p, t)
    A, B, c, b = u
    du[1] = -p.alpha * A + p.beta * B
    du[2] = p.alpha * A - p.beta * B
    bind = p.kon * c * (p.bmax - b) - p.koff * b
    du[3] = -bind
    du[4] = bind
    nothing
end
