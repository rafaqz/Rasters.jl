abstract type AbstractModifications end
struct NoMod{T,Mi} <: AbstractModifications
     missingval::Mi
end
NoMod{T}(missingval::Mi) where {T,Mi} = NoMod{T,Mi}(missingval)
NoMod() = NoMod{Any}()
NoMod{T}() where T = NoMod{T}(nothing)
NoMod{T}(::NoKW) where T = NoMod{T}(nothing)

Base.eltype(::NoMod{T}) where T = T
source_eltype(::NoMod{T}) where T = T

struct Mod{T1,T2,Mi,Ma,S,O,F} <: AbstractModifications
     missingval::Mi
     coalesceval::Ma
     scale::S
     offset::O
     coerce::F
     function Mod(::Type{T}, missingval, coalesceval, scale, offset, coerce) where T
         coalesceval = coalesceval === missingval ? nothing : coalesceval
         if isnokw(coerce) || isnothing(coerce)
             coerce = convert
         end
         vals = map(_nokw2nothing, (missingval, coalesceval, scale, offset))
         T1 = _resolve_mod_eltype(T, vals...)
         new{T1,T,map(typeof, vals)...,typeof(coerce)}(vals..., coerce)
     end
end

Base.eltype(::Mod{T1}) where T1 = T1
source_eltype(::Mod{<:Any,T2}) where T2 = T2


function _resolve_mod_eltype(::Type{T}, missingval, coalesceval, scale, offset) where T
    T1 = isnothing(coalesceval) ? T : promote_type(T, typeof(coalesceval))
    T2 = isnothing(scale) ? T1 : promote_type(T1, typeof(scale))
    T3 = isnothing(offset) ? T2 : promote_type(T2, typeof(offset))
    return T3
end

missingval(m::Mod) = m.missingval
coalesceval(m::Mod) = isnothing(m.coalesceval) ? m.missingval : m.coalesceval
missingval(m::NoMod) = m.missingval
coalesceval(m::NoMod) = missingval(m)

struct ModifiedDiskArray{I,T,N,V,M} <: DiskArrays.AbstractDiskArray{T,N}
    var::V
    mod::M
end
function ModifiedDiskArray(v::V, m::M; invert=false) where {V<:AbstractArray{<:Any,N},M} where N
    T = invert ? source_eltype(m) : eltype(m)
    return ModifiedDiskArray{invert,T,N,V,M}(v, m)
end

Base.parent(A::ModifiedDiskArray) = A.var
Base.size(A::ModifiedDiskArray, args...) = size(parent(A), args...)
filename(A::ModifiedDiskArray) = filename(parent(A))
missingval(A::ModifiedDiskArray) = A.missingval
coalesceval(A::ModifiedDiskArray) = A.coalesceval
DiskArrays.haschunks(A::ModifiedDiskArray) = DiskArrays.haschunks(parent(A))
DiskArrays.eachchunk(A::ModifiedDiskArray) = DiskArrays.eachchunk(parent(A))

function DiskArrays.readblock!(
    A::ModifiedDiskArray{false,<:Any,0}, out_block, I::AbstractVector...
)
    out_block[] = _applymod(parent(A)[I...][], A.mod)
    return out_block
end
function DiskArrays.readblock!(
    A::ModifiedDiskArray{true,T,<:Any,0}, out_block, I::AbstractVector...
) where T
    out_block[] = _invertmod(Val{T}(), parent(A)[I...], A.mod)
    return out_block
end
function DiskArrays.readblock!(
    A::ModifiedDiskArray{false}, out_block, I::AbstractVector...
)
    out_block .= _applymod.(parent(A)[I...], (A.mod,))
    return out_block
end
function DiskArrays.readblock!(
    A::ModifiedDiskArray{true,T}, out_block, I::AbstractVector...
) where T
    out_block .= _invertmod.((Val{T}(),), parent(A)[I...], (A.mod,))
    return out_block
end

function DiskArrays.writeblock!(
    A::ModifiedDiskArray{false,<:Any,0,<:AbstractArray{T}}, block, I::AbstractVector...
) where T

    parent(A)[I...] = _invertmod(Val{source_eltype(A.mod)}(), block[], A.mod)
    return nothing
end
function DiskArrays.writeblock!(
    A::ModifiedDiskArray{true,<:Any,0,<:AbstractArray{T}}, _block, I::AbstractVector...
) where T
    parent(A)[I...] = _applymod(Val{eltype(A.mod)}(), block[], A.mod)
    return nothing
end
function DiskArrays.writeblock!(
    A::ModifiedDiskArray{<:Any,<:Any,<:Any,<:AbstractArray{T}}, block, I::AbstractVector...
) where T
    parent(A)[I...] = _invertmod.((Val{source_eltype(A.mod)}(),), block, (A.mod,))
    return nothing
end
function DiskArrays.writeblock!(
    A::ModifiedDiskArray{true,<:Any,<:Any,<:AbstractArray{T}}, _block, I::AbstractVector...
) where T
    parent(A)[I...] = _applymod.((Val{eltype(A.mod)}(),), block, (A.mod,))
    return nothing
end

Base.@assume_effects :foldable function _applymod(x, m::Mod)
    if _ismissing(x, missingval(m))
        coalesceval(m)
    else
        _scaleoffset(x, m)
    end
end
Base.@assume_effects :foldable _applymod(x, m::NoMod) = x

_ismissing(x, mv) = isequal(x, mv)
_ismissing(_, ::Nothing) = false

_scaleoffset(x, m::Mod) = _scaleoffset(x, m.scale, m.offset)
_scaleoffset(x, scale, offset) = x * scale + offset
_scaleoffset(x, ::Nothing, offset) = x + offset
_scaleoffset(x, scale, ::Nothing) = x * scale
_scaleoffset(x, ::Nothing, ::Nothing) = x

Base.@assume_effects :foldable function _invertmod(::Val{T}, x, m::Mod) where T
    tm = if isnothing(m.missingval)
        x
    else
        if _ismissing(x, m.coalesceval)
            return m.missingval
        else
            x
        end
    end
    return _scaleoffset_inv(T, tm, m)
end
Base.@assume_effects :foldable _invertmod(v, x, m::NoMod) = x

Base.@assume_effects :foldable _scaleoffset_inv(::Type{T}, x, m::Mod) where T = 
    _scaleoffset_inv(m.coerce, T, x, m)::T
Base.@assume_effects :foldable _scaleoffset_inv(coerce::Base.Callable, ::Type{T}, x, m::Mod) where T =
    coerce(T, _scaleoffset_inv1(x, m.scale, m.offset))::T

Base.@assume_effects :foldable _scaleoffset_inv1(x, scale, offset) = (x - offset) / scale
Base.@assume_effects :foldable _scaleoffset_inv1(x, scale, ::Nothing) = x / scale
Base.@assume_effects :foldable _scaleoffset_inv1(x, ::Nothing, offset) = x - offset
Base.@assume_effects :foldable _scaleoffset_inv1(x, ::Nothing, ::Nothing) = x


function _stack_mods(
    eltypes::Vector, metadata::Vector, missingval::Vector, coalesceval;
    scaled, coerce
)
    map(eltypes, metadata, missingval) do T, md, mv
        scale, offset = get_scale(md, scaled)
        _mod(T, mv, coalesceval, scale, offset, coerce)
    end
end
function _stack_mods(
    eltypes::Vector, metadata::Vector, missingval, coalesceval::Vector;
    scaled::Bool, coerce
)
    map(eltypes, metadata, coalesceval) do T, md, mk
        scale, offset = get_scale(md, scaled)
        _mod(T, missingval, mk, scale, offset, coerce)
    end
end
function _stack_mods(
    eltypes::Vector, metadata::Vector, missingval::Vector, coalesceval::Vector;
    scaled::Bool, coerce
)
    map(eltypes, metadata, missingval, coalesceval) do T, md, mv, mk
        scale, offset = get_scale(md, scaled)
        _mod(mv, mk, scale, offset, coerce)
    end
end
function _stack_mods(
    eltypes::Vector, metadata::Vector, missingval, coalesceval;
    scaled::Bool, coerce
)
    map(eltypes, metadata) do T, md
        scale, offset = get_scale(md, scaled)
        _mod(T, missingval, coalesceval, scale, offset, coerce)
    end
end

function _mod(::Type{T}, metadata, missingval, coalesceval; scaled::Bool, coerce) where T
    scale, offset = get_scale(metadata, scaled)
    _mod(T, missingval, coalesceval, scale, offset, coerce)
end
function _mod(::Type{T}, missingval, coalesceval, scale, offset, coerce) where T
    coalesceval = if isnokw(coalesceval)
        # If there is no missingval dont mask
        isnokwornothing(missingval) ? nothing : missing
    else
        # Unless coalesceval was passed explicitly
        coalesceval === missingval ? nothing : coalesceval
    end
    if isnokwornothing(coalesceval) && isnokwornothing(scale) && isnokwornothing(offset)
        return NoMod{T}(missingval)
    else
        return Mod(T, missingval, coalesceval, scale, offset, coerce)
    end
end

@inline get_scale(metadata::NoKW, scaled::Bool) = nothing, nothing
@inline function get_scale(metadata, scaled::Bool)
    scale = scaled ? get(metadata, "scale", nothing) : nothing
    offset = scaled ? get(metadata, "offset", nothing) : nothing
    return scale, offset
end

function _writer_mod(::Type{T}; missingval, coalesceval, scale, offset, coerce) where T
    missingval1 = if isnokw(missingval) || isnothing(missingval)
        if isnokw(coalesceval) || isnothing(coalesceval)
            nothing
        else
            _type_missingval(T)
        end
    elseif ismissing(missingval)
        _type_missingval(T)
    else
        missingval
    end
    coalesceval1 = if isnokw(coalesceval)
        if Missing <: T
            missing
        else
            nothing
        end
    else
        coalesceval
    end
    return _mod(T, missingval1, coalesceval1, scale, offset, coerce)
end

_mod_eltype(::AbstractArray{T}, ::NoMod) where T = T
_mod_eltype(::AbstractArray, m::Mod{T}) where T = T

_mod_inverse_eltype(::AbstractArray{T}, ::NoMod) where T = T
_mod_inverse_eltype(::AbstractArray{T}, m::Mod) where T =
    Base.promote_op(_invertmod, typeof(m.coerce), T, typeof(m))

_maybe_modify(var, m::Mod; kw...) = ModifiedDiskArray(var, m; kw...)
_maybe_modify(var, ::NoMod; kw...) = var
