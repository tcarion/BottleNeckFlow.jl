@recipe(CanalPlot) do scene
    Attributes(
        color = :red
    )
end

@recipe(GriddedPlot) do scene
    Attributes(
        colormap = :jet
    )
end

@recipe(PlotGrid) do scene
    Attributes(
        marker = :square
    )
end

const CanalPlotArg = CanalPlot{<:Tuple{AbstractVector, CanalConfig}}
const GriddedPlotArg = GriddedPlot{<:Tuple{AbstractGrid}}
const PlotGridArg = PlotGrid{<:Tuple{AbstractGrid}}
# Makie.argument_names(::Type{<: CanalPlotArg}) = (:x, :canal,)

function Makie.plot!(cp::CanalPlotArg)
    xs = cp[1]
    canal = cp[2]
    newy = δ(xs.val, canal.val)

    h = canal.val.H
    hobs = Observable(h)
    ytop = Observable(h .- newy)

    ybot = Observable(-h .+ newy)
    # ys[] = newy
    # d[] = canal.val.d
    function update_plot(xs, canal)
        # colors[]

        # clear tde vectors inside tde observables
        empty!(ytop[])
        empty!(ybot[])
        newy = δ(xs, canal)
        hobs[] = canal.H
        ytop[] = hobs[] .- newy
        ysab[] = -hobs[] .+ newy
        # tden refill tdem witd our updated values
        # for (t, s) in zip(times, stockvalues)
        #     pusd!(ys[], Point2f(t, s.low))
        # end
        # append!(colors[], [x.close > x.open for x in stockvalues])
        # colors[] = colors[]
    end

    Makie.Observables.onany(update_plot, xs, canal)

    lines!(cp, xs, ytop, color = cp.color)
    lines!(cp, xs, ybot, color = cp.color)
    cp
end

function Makie.plot!(gp::GriddedPlotArg)
    grid = gp[1]

    # mat = Observable(trueview(grid.val))

    # function update_plot(grid)
    #     empty!(mat[])
    #     mat[] = trueview(grid)
    # end

    # Makie.Observables.onany(update_plot, grid)

    # z2nan = map(grid) do g
    #     replace(x -> isapprox(0., x) ? NaN : x, g)
    # end
    contourf!(gp, getxs(grid.val), getys(grid.val), grid, colormap = gp.colormap)
end

function Makie.plot!(gp::PlotGridArg)
    grid = gp[1]

    scatter!(gp, [(x, y) for x in getxs(grid.val) for y in getys(grid.val)], marker = gp.marker)
end