function convectionCenteredX(u::UGrid, v::VGrid, i, j)
    uleft = u[i-1, j]
    uright = u[i, j]
    vbot = v[i, j-1]
    vtop = v[i, j]
    ubar = 0.5 * (uleft + uright)
    vbar = 0.5 * (vbot + vtop)

    diffu = uright - uleft

    ubar * diff / u.grid.dx + vbar
end