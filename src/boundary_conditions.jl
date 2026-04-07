@inline function map_index(k, l, nx, ny, ::FreeBoundary)
    if 1 <= k <= nx && 1 <= l <= ny
        return k::Int, l::Int, true::Bool
    else
        return 0::Int, 0::Int, false::Bool
    end
end

@inline function map_index(k, l, nx, ny, ::PeriodicBoundary)
    # I write long branch statements because they should be faster than return mod1(k, nx), mod1(l, ny)
    if k<1
        kk = k + nx
    elseif k>nx
        kk = k-nx
    else
        kk = k
    end

    if l<1
        ll = l+ny
    elseif l>ny
        ll = l-ny
    else
        ll = l
    end

    return kk, ll, true
end