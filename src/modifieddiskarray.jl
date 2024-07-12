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

function _mod(cf::Bool, metadata; missingval, maskingval, coerce=convert)
    scale = cf ? get(metadata, "scale", nothing) : nothing
    offset = cf ? get(metadata, "offset", nothing) : nothing
    _mod(missingval, maskingval, scale, offset, coerce)
end
function _mod(missingval, maskingval, scale, offset, coerce=convert)
    if isnothing(maskingval) && isnothing(scale) && isnothing(offset)
        return NoMod(missingval)
    else
        return Mod(missingval, maskingval, scale, offset, coerce)
    end
end

_mod_eltype(::AbstractArray{T}, ::NoMod) where T = T
_mod_eltype(::AbstractArray{T}, m::Mod) where T =
    Base.promote_op(_applymod, T, typeof(m))

_mod_inverse_eltype(::AbstractArray{T}, ::NoMod) where T = T
_mod_inverse_eltype(::AbstractArray{T}, m::Mod) where T =
    Base.promote_op(_invertmod, typeof(m.coerce), T, typeof(m))

_maybe_modify(var, m::Mod) = ModifiedDiskArray(var, m)
_maybe_modify(var, ::NoMod) = var

struct ModifiedDiskArray{T,N,V,M} <: DiskArrays.AbstractDiskArray{T,N}
    var::V
    mod::M
end
function ModifiedDiskArray(v::V, m::M) where {V<:AbstractArray{<:Any,N},M} where N
    T = _mod_eltype(v, m)
    return ModifiedDiskArray{T,N,V,M}(v, m)
end

Base.parent(A::ModifiedDiskArray) = A.var
Base.size(A::ModifiedDiskArray, args...) = size(A.var, args...)
DiskArrays.haschunks(A::ModifiedDiskArray) = DiskArrays.haschunks(A.var)
DiskArrays.eachchunk(A::ModifiedDiskArray) = DiskArrays.eachchunk(A.var)

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
_scaleoffset(x, scale, offset) = muladd(x, scale, offset)
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
