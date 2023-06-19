const DEFAULT_LOGGER = SimpleLogger

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
getuc(pa::SimParam, gb::GridBox) = 1.5*pa.Q/gb.D
getum(pa::SimParam, gb::GridBox) = 2. / 3. * getuc(pa, gb)
getnu(pa::SimParam, gb::GridBox) = getuc(pa, gb) * gb.D / pa.Re
function Base.show(io::IO, params::SimParam)
    println(io, "### SimParams:")
    println(io, "\t Reynolds number: $(params.Re)")
    println(io, "\t In-flow [m2/s]: $(params.Q)")
    println(io, "\t CFL: $(params.cfl)")
    println(io, "\t Relaxation coefficient: $(params.alpha)")
    println(io, "\t dt/dτ: $(params.dtdτ)")
    println(io, "\t Tol for GS: $(params.tol)")
    println(io, "\t Ramp for khi: $(params.ramp)")
end
from_ncf(::Type{SimParam}, fname::String) = _convert(SimParam, fname)


mutable struct Simulation
    u::UGrid
    v::VGrid
    p::PGrid
    params::SimParam
    dt::Real
    canal::CanalConfig
    logger::AbstractLogger
    step::Int
    uprev::Union{Nothing, UGrid}
    vprev::Union{Nothing, VGrid}
end

function Simulation(params::SimParam, gridbox::GridBox, canal::CanalConfig, logger = DEFAULT_LOGGER())
    dt = params.cfl * gridbox.dx / getuc(params, gridbox)
    Simulation(
        UGrid(gridbox),
        VGrid(gridbox),
        PGrid(gridbox),
        params,
        dt,
        canal,
        logger,
        0,
        nothing,
        nothing
        )
end
function Simulation(fname::String, logger = DEFAULT_LOGGER())
    gb = from_ncf(GridBox, fname)
    params = from_ncf(SimParam, fname)
    canal = from_ncf(CanalConfig, fname)
    # ugrid = UGrid(gb)
    # vgrid = VGrid(gb)
    # pgrid = PGrid(gb)
    Dataset(fname, "r") do ds
        dt = ds.attrib["dt"]
        step = ds["steps"][end]
        ugrid = UGrid(ds["u"][:,:,end], gb)
        vgrid = VGrid(ds["v"][:,:,end], gb)
        pgrid = PGrid(ds["p"][:,:,end], gb)
        Simulation(
            ugrid,
            vgrid,
            pgrid,
            params,
            dt,
            canal,
            logger,
            step,
            nothing,
            nothing
        )
    end
end

function Base.show(io::IO, sim::Simulation)
    println(io, "##### Simulation:")
    println(io, "\t dt [s]: $(sim.dt)")
    println(io, "\t nstep: $(sim.step)")
    println(io, "\t Max speed [m/s]: $(getuc(sim.params, sim.u.grid))")
    println(io, "\t nu [m^2/s]: $(getnu(sim.params, sim.u.grid))")
    println(io, "\t Convetive time [s]: $(sim.u.grid.L / getum(sim.params, sim.u.grid))")

    show(io, sim.params)
    show(io, sim.u.grid)
    show(io, sim.canal)
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
        sim.step,
        nothing,
        nothing
    )
end