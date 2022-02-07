@recipe(CanalPlot) do scene
    Attributes(
        color = :red
    )
end

const CanalPlotArg = CanalPlot{<:Tuple{AbstractVector, CanalConfig}}
# Makie.argument_names(::Type{<: CanalPlotArg}) = (:x, :canal,)

function Makie.plot!(cp::CanalPlotArg)
    xs = cp[1]
    canal = cp[2]
    newy = δ(xs.val, canal.val)

    ys = Observable(newy)
    d = Observable(canal.val.d)

    ysab = Observable(canal.val.d .- newy)
    # ys[] = newy
    # d[] = canal.val.d
    function update_plot(xs, canal)
        # colors[]

        # clear tde vectors inside tde observables
        empty!(ys[])
        newy = δ(xs, canal)
        ys[] = newy
        d[] = canal.d
        ysab[] = canal.d .- newy
        # tden refill tdem witd our updated values
        # for (t, s) in zip(times, stockvalues)
        #     pusd!(ys[], Point2f(t, s.low))
        # end
        # append!(colors[], [x.close > x.open for x in stockvalues])
        # colors[] = colors[]
    end

    Makie.Observables.onany(update_plot, xs, canal)

    lines!(cp, xs, ys, color = cp.color)
    lines!(cp, xs, ysab, color = cp.color)
    cp
end