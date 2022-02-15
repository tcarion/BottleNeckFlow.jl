function runsim!(sim::Simulation, n)
    sim.step !== 0 && error("Simulation is not at initial state.") 
    poiseuille!(sim)
    boundaries!(sim)

    un_1 = copy(sim.u)
    vn_1 = copy(sim.v)
    ti = @elapsed begin
        make_step!(sim)
    end
    with_logger(sim.logger) do 
        @info "Step number $step done in $(round(ti * 1.0e3, digits = 2))ms\nSimulation time : $(sim.step * sim.dt)"
    end
    for step in 2:n
        ti = @elapsed begin
            make_step!(sim, un_1, vn_1)
        end
        with_logger(sim.logger) do 
            @info """Step number $step done in $(round(ti * 1.0e3, digits = 2))ms\nSimulation time : $(sim.step * sim.dt)
            Diagnosis :
            Re_h = $(reynh(sim))
            Re_h_omega = $(reynhvort(sim))
            """
        end
        un_1 = copy(sim.u)
        vn_1 = copy(sim.v)
    end
end

function make_step!(sim::Simulation)
    # if sim.step == 0
    step_euler!(sim)
    # else

    step_poisson!(sim, sim.params.tol)
    poisson_project!(sim)
    boundaries!(sim)
    sim.step = sim.step + 1
end

function make_step!(sim::Simulation, un_1::UGrid, vn_1::VGrid)
    # if sim.step == 0
    step_AB2!(sim, un_1, vn_1)
    # else

    step_poisson!(sim, sim.params.tol)
    poisson_project!(sim)
    boundaries!(sim)
    sim.step = sim.step + 1
end