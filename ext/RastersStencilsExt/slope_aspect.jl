
abstract type SlopeFilter end
abstract type SlopeConvolution <: SlopeFilter end
# Not all slope algorithms can provide aspect
abstract type SlopeAspectConvolution <: SlopeConvolution end

struct FD2 <: SlopeAspectConvolution end
struct FD3Reciprocal <: SlopeAspectConvolution end
struct FD3ReciprocalSquared <: SlopeAspectConvolution end
struct FD3Linear <: SlopeAspectConvolution end
struct FDFrame <: SlopeAspectConvolution end
struct SimpleDifference <: SlopeConvolution end

@inline function aspect_filter(method::SlopeConvolution, n::Window, d)
    fx, fy = _slope_conv(method, n, d)
    return _aspect(fx, fy)
end

@inline function slope_filter(method::SlopeConvolution, n::Window, d)
    fx, fy = _slope_conv(method, n, d)
    return _slope(fx, fy)
end

@inline function slopeaspect_filter(method::SlopeAspectConvolution, n::Window, d)
    fx, fy = _slope_conv(method, n, d)
    return _slope(fx, fy), _aspect(fx, fy)
end

_slope(fx, fy) = atan(√(fx^2 + fy^2))
function _aspect(fx, fy)
    (ismissing(fx) || ismissing(fy)) && return missing
    # Rotate - we want high Y (north) as the origin
    # TODO: pass through the Order for X/Y dims 
    # So the result always has zero at North
    -atan(fx, fy) 
end

@inline function _slope_conv(::FD2, n::Window, d)
    fx = (n[6] - n[4]) / 2d
    fy = (n[8] - n[2]) / 2d
    return fx, fy
end

@inline function _slope_conv(::FD3Reciprocal, n::Window, d)
    fx = (n[3] -n[1] + √(2(n[6] - n[4])) + n[9] - n[7]) / (4 + 2 * √(2)) * d
    fy = (n[7] -n[1] + √(2(n[8] - n[2])) + n[9] - n[3]) / (4 + 2 * √(2)) * d
    return fx, fy
end

@inline function _slope_conv(::FD3Linear, n::Window, d)
    fx = (n[3] - n[1] + n[6] - n[4] + n[9] - n[7]) / 6d
    fy = (n[7] - n[1] + n[8] - n[2] + n[9] - n[3]) / 6d
    return fx, fy
end

@inline function _slope_conv(::FD3ReciprocalSquared, n::Window, d)
    fx = (n[3] - n[1] + 2(n[6] - n[4]) + n[9] - n[7]) / 8d
    fy = (n[7] - n[1] + 2(n[8] - n[2]) + n[9] - n[3]) / 8d
    return fx, fy
end

@inline function _slope_conv(::FDFrame, n::Window, d)
    fx = (n[3] - n[1] + n[9] - n[7]) / 4d
    fy = (n[7] - n[1] + n[9] - n[3]) / 4d
    return fx, fy
end

@inline function _slope_conv(::SimpleDifference, n::Window, d)
    fy = (n[5] - n[2]) / d
    return fy, fy # No aspect
end

struct MaxSlope <: SlopeFilter end

@inline function slope_filter(method::MaxSlope, n::Window, g)
    # slopes = (
    #     abs((g - n[2]) / g), 
    #     abs((g - n[4]) / g),
    #     abs((g - n[6]) / g), 
    #     abs((g - n[8]) / g), 
    #     abs((g - n[1]) / (√(2)*g)), 
    #     abs((g - n[3]) / (√(2)*g)), 
    #     abs((g - n[7]) / (√(2)*g)), 
    #     abs((g - n[9]) / (√(2)*g)), 
    # )
    # xmissing = map(ismissing, slopes)
    # if all(xmissing)
    #     return missing
    # elseif any(xmissing)
    #     return maximum(skipmissing(slopes))
    # else
    #     return maximum(slopes)
    # end
end

function slope end

function aspect end

"""
    slope_aspect(elevation::Raster, method; [cellsize=1])
    
Calculate both slope and aspect. Since these are calculated together anyway, 
it's more efficient to use a single function call if you do want both 
slope and aspect.

Returns a `RasterStack` with `:slope` and `:aspect` layers.
"""
function slope_aspect end


for (f, filt) in (:slope => :slope_filter, :aspect => :aspect_filter, :_slopeaspect => :slopeaspect_filter)
    @eval begin 
        function $(f)(elevation::AbstractRaster, method=FD2(); cellsize=1)
            stencil = Window{1}(); 
            padval = missingval(elevation)
            newdata = Neighborhoods.broadcast_stencil(stencil, parent(elevation); boundary=Remove(padval)) do w
                $(filt)(method, w, cellsize)
            end
            rebuild(elevation; data=newdata, name=$(QuoteNode(f))) 
        end
    end
end

function slopeaspect(elevation, method=FD2(); cellsize=1)
    sa = _slopeaspect(elevation, method; cellsize)
    slope = first.(sa)
    aspect = last.(sa)
    nt = (; slope, aspect)
    return RasterStack(nt)
end
