abstract type AbstractModifications end
struct NoMod{Mi} <: AbstractModifications
     missingval::Mi
end
NoMod() = NoMod(nothing)
NoMod(::NoKW) = NoMod(nothing)
struct Mod{Mi,Ma,S,O,F} <: AbstractModifications
     missingval::Mi
     maskingval::Ma
     scale::S
     offset::O
     coerce::F
     function Mod(missingval, maskingval, scale, offset, coerce)
         if isnokw(coerce) || isnothing(coerce) 
             coerce = convert
         end
         vals = map(_nokw2nothing, (missingval, maskingval, scale, offset))
         new{map(typeof, vals)...,typeof(coerce)}(vals..., coerce)
     end
end

missingval(m::Mod) = m.missingval
maskingval(m::Mod) = m.maskingval
missingval(m::NoMod) = m.missingval
maskingval(m::NoMod) = nothing


struct ModifiedDiskArray{T,N,V,M} <: DiskArrays.AbstractDiskArray{T,N}
    var::V
    mod::M
end
function ModifiedDiskArray(v::V, m::M) where {V<:AbstractArray{<:Any,N},M} where N
    T = _mod_eltype(v, m)
    return ModifiedDiskArray{T,N,V,M}(v, m)
end

Base.parent(A::ModifiedDiskArray) = A.var
Base.size(A::ModifiedDiskArray, args...) = size(parent(A), args...)
filename(A::ModifiedDiskArray) = filename(parent(A))
missingval(A::ModifiedDiskArray) = A.missingval
maskingval(A::ModifiedDiskArray) = A.maskingval
DiskArrays.haschunks(A::ModifiedDiskArray) = DiskArrays.haschunks(parent(A))
DiskArrays.eachchunk(A::ModifiedDiskArray) = DiskArrays.eachchunk(parent(A))

function DiskArrays.readblock!(A::ModifiedDiskArray, out_block, I::AbstractVector...)
    broadcast!(_applymod, out_block, A.var[I...], (A.mod,))
    return nothing
end

function DiskArrays.writeblock!(
    A::ModifiedDiskArray{<:Any,<:Any,<:AbstractArray{T}}, in_block, I::AbstractVector...
) where T
    A.var[I...] = _invertmod.((Val{T}(),), in_block, (A.mod,))
    return nothing
end

Base.@assume_effects :foldable function _applymod(x, m::Mod)
    tm = if isnothing(m.maskingval)
        x
    else
        if _ismissing(x, m.missingval)
            return m.maskingval
        else
            x
        end
    end
    return _scaleoffset(tm, m)
end

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
        if _ismissing(x, m.maskingval)
            return m.missingval
        else
            x
        end
    end
    return _scaleoffset_inv(T, tm, m)
end

_scaleoffset_inv(::Type{T}, x, m::Mod) where T = _scaleoffset_inv(m.coerce, T, x, m)
_scaleoffset_inv(coerce::Base.Callable, ::Type{T}, x, m::Mod) where T = 
    coerce(T, _scaleoffset_inv(x, m.scale, m.offset))
_scaleoffset_inv(x, scale, offset) = (x - offset) / scale
_scaleoffset_inv(x, scale, ::Nothing) = x / scale
_scaleoffset_inv(x, ::Nothing, offset) = x - offset
_scaleoffset_inv(x, ::Nothing, ::Nothing) = x


function _stack_mods(metadata::Vector, missingval::Vector, maskingval; scaled, coerce)
    map(metadata, missingval) do md, mv
        scale, offset = _get_scale_offset(md, scaled)
        _mod(mv, maskingval, scale, offset, coerce)
    end
end
function _stack_mods(metadata::Vector, missingval, maskingval::Vector; scaled::Bool, coerce)
    map(metadata, maskingval) do md, mk
        scale, offset = _get_scale_offset(md, scaled)
        _mod(missingval, mk, scale, offset, coerce)
    end
end
function _stack_mods(metadata::Vector, missingval::Vector, maskingval::Vector; scaled::Bool, coerce)
    map(metadata, missingval, maskingval) do md, mv, mk
        scale, offset = _get_scale_offset(md, scaled)
        _mod(mv, mk, scale, offset, coerce)
    end
end
function _stack_mods(metadata::Vector, missingval, maskingval; scaled::Bool, coerce)
    map(metadata) do md
        scale, offset = _get_scale_offset(md, scaled)
        _mod(missingval, maskingval, scale, offset, coerce)
    end
end

function _mod(metadata, missingval, maskingval; scaled::Bool, coerce)
    scale, offset = _get_scale_offset(metadata, scaled)
    _mod(missingval, maskingval, scale, offset, coerce)
end
function _mod(missingval, maskingval, scale, offset, coerce)
    if isnothing(maskingval) && isnothing(scale) && isnothing(offset)
        return NoMod(missingval)
    else
        return Mod(missingval, maskingval, scale, offset, coerce)
    end
end

@inline _get_scale_offset(metadata::NoKW, scaled) = (nothing, nothing)
@inline function _get_scale_offset(metadata, scaled)
    scale = scaled ? get(metadata, "scale", nothing) : nothing
    offset = scaled ? get(metadata, "offset", nothing) : nothing
    return scale, offset
end

function _writer_mod(::Type{T}; missingval, maskingval, scale, offset, coerce) where T
    missingval1 = if isnokw(missingval) || isnothing(missingval)
        if isnokw(maskingval) || isnothing(maskingval)
            nothing
        else
            _type_missingval(T) 
        end
    elseif ismissing(missingval)
        _type_missingval(T) 
    else
        missingval
    end
    maskingval1 = if isnokw(maskingval)
        if Missing <: T
            missing
        else
            nothing
        end
    else
        maskingval
    end
    return _mod(missingval1, maskingval1, scale, offset, coerce)
end

_mod_eltype(::AbstractArray{T}, ::NoMod) where T = T
_mod_eltype(::AbstractArray{T}, m::Mod) where T =
    Base.promote_op(_applymod, T, typeof(m))

_mod_inverse_eltype(::AbstractArray{T}, ::NoMod) where T = T
_mod_inverse_eltype(::AbstractArray{T}, m::Mod) where T =
    Base.promote_op(_invertmod, typeof(m.coerce), T, typeof(m))

_maybe_modify(var, m::Mod) = ModifiedDiskArray(var, m)
_maybe_modify(var, ::NoMod) = var
