Base.@kwdef struct SimParam
    # Reynolds number
    Re::Real = 250.
    # In-flow [m2/s]
    Q::Real = 1.
    # CFL
    cfl::Real = 0.5
end
getuc(pa::SimParam, gb::GridBox) = 0.75*pa.Q/gb.D
getnu(pa::SimParam, gb::GridBox) = getuc(pa, gb) * gb.D / pa.Re

struct Simulation
    u::UGrid
    v::VGrid
    p::PGrid
    params::SimParam
    dt::Real
end

function Simulation(params::SimParam, gridbox::GridBox)
    dt = params.cfl * gridbox.dx / getuc(params, gridbox)
    Simulation(
        UGrid(gridbox),
        VGrid(gridbox),
        PGrid(gridbox),
        params,
        dt
        )
end