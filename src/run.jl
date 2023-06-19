function runsim(gb::GridBox, canal::CanalConfig, n, fname = ""; save = 50)
    sim = Simulation(SimParam(), gb, canal)
    runsim!(sim, n, fname; save = save)
    sim
end
function runsim!(sim::Simulation, n, fname = ""; save = 50)
    lograte = 50
    sim.step !== 0 && error("Simulation is not at initial state.") 
    
    poiseuille!(sim)
    boundaries!(sim)
    !isempty(fname) && create(fname, sim)
    
    # un_1 = copy(sim.u)
    # vn_1 = copy(sim.v)
    # ti = @elapsed begin
    #     make_step!(sim)
    # end
    # sim.uprev = copy(sim.u)
    # sim.vprev = copy(sim.v)
    # with_logger(sim.logger) do 
    #     @info "Step number $step done in $(round(ti * 1.0e3, digits = 2))ms\nSimulation time : $(sim.step * sim.dt)"
    # end
    # saveif(fname, sim, save)

    for step in 1:n
        ti = @elapsed begin
            make_step!(sim)
        end
        sim.uprev = copy(sim.u)
        sim.vprev = copy(sim.v)

        if sim.step % lograte == 0
            with_logger(sim.logger) do 
                @info """Step number $step done in $(round(ti * 1.0e3, digits = 2))ms\nSimulation time : $(sim.step * sim.dt)
                Diagnosis :
                Re_h = $(reynh(sim))
                Re_h_omega = $(reynhvort(sim))
                """
            end
        end
        # ti = @elapsed begin
        #     make_step!(sim, un_1, vn_1)
        # end

        saveif(fname, sim, save)

        # un_1 = copy(sim.u)
        # vn_1 = copy(sim.v)
    end
    !isempty(fname) && add_to_nc(fname, sim)
end

function saveif(fname, sim, save)
    !isempty(fname) && (sim.step % save == 0) && add_to_nc(fname, sim)
end

function make_step!(sim::Simulation)
    # if sim.step == 0
    # step_euler!(sim)
    step!(sim)
    # else
    boundaries!(sim)
    step_poisson!(sim, sim.params.tol)
    poisson_project!(sim)
    boundaries!(sim)
    sim.step = sim.step + 1
end

# function make_step!(sim::Simulation, un_1::UGrid, vn_1::VGrid)
#     # if sim.step == 0
#     step_AB2!(sim, un_1, vn_1)
#     # else
#     boundaries!(sim)
#     step_poisson!(sim, sim.params.tol)
#     poisson_project!(sim)
#     boundaries!(sim)
#     sim.step = sim.step + 1
# end