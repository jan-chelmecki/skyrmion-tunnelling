@inline function map_index(k, l, nx, ny, ::FreeBoundary)
    if 1 <= k <= nx && 1 <= l <= ny
        return k::Int, l::Int, true::Bool
    else
        return 0::Int, 0::Int, false::Bool
    end
end

@inline function map_index(k, l, nx, ny, ::PeriodicBoundary)
    return mod1(k, nx), mod1(l, ny), true
end