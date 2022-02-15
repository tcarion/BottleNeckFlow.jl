Base.@kwdef mutable struct SimParam
    # Reynolds number
    Re::Real = 250.
    # In-flow [m2/s]
    Q::Real = 1.
    # CFL
    cfl::Real = 0.5
    # Relaxation coefficient
    alpha::Real = 1.97
    # Ramp time for bump
    dtdτ::Real = 1.
    # Tolerance for Gauss Seidel
    tol::Real = 1e-3
    # Ramp for khi
    ramp::Int = 200
end
getuc(pa::SimParam, gb::GridBox) = 0.75*pa.Q/gb.D
getum(pa::SimParam, gb::GridBox) = pa.Q / gb.D
getnu(pa::SimParam, gb::GridBox) = getuc(pa, gb) * gb.D / pa.Re

mutable struct Simulation
    u::UGrid
    v::VGrid
    p::PGrid
    params::SimParam
    dt::Real
    canal::CanalConfig
    logger::AbstractLogger
    step::Int
end

function Simulation(params::SimParam, gridbox::GridBox, canal::CanalConfig, logger = SimpleLogger())
    dt = params.cfl * gridbox.dx / getuc(params, gridbox)
    Simulation(
        UGrid(gridbox),
        VGrid(gridbox),
        PGrid(gridbox),
        params,
        dt,
        canal,
        logger,
        0
        )
end

function Base.copy(sim::Simulation)
    Simulation(
        copy(sim.u),
        copy(sim.v),
        copy(sim.p),
        sim.params,
        sim.dt,
        sim.canal,
        sim.logger,
        sim.step
    )
end