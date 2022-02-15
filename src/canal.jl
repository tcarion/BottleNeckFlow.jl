struct CanalConfig
    H::Real
    hd::Real
    Ld::Real
    x0::Real
    L::Real
    d::Real
end

CanalConfig(H::Real, hd::Real, Ld::Real, x0::Real, L::Real) = CanalConfig(H, hd, Ld, x0, L, 2*H)
CanalConfig(H::Real) = CanalConfig(H, 0.55*H, 2*H, 4*H, 16*H)

function δ(x::Real, c::CanalConfig)
    diff = x - c.x0
    l = c.Ld * 0.5
    if (-l <= diff) && (diff <= l)
        tr = cos(2pi * (x - c.x0) / c.Ld)
        c.hd * ( (1 + tr) * 0.5 )^2
    else
        0.
    end
end

δ(x::AbstractVector, c::CanalConfig) = δ.(x, Ref(c))

function inbump(c::CanalConfig, x, y)
    diff = abs(x - c.x0)
    l = c.Ld * 0.5
    if diff <= l
        tr = cos(2pi * (x - c.x0) / c.Ld)
        yc = c.hd * ( (1 + tr) * 0.5 )^2
        if c.H - abs(y) <= yc
            true
        else
            false
        end
    else
        false
    end
end