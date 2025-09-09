
abstract type AbstractModifications{T} end

struct NoMod{T,Mi} <: AbstractModifications{T}
     missingval::Mi
end
NoMod{T}(missingval::Mi) where {T,Mi} = NoMod{T,Mi}(missingval)
NoMod() = NoMod{Any}()
NoMod{T}() where T = NoMod{T}(nothing)
NoMod{T}(::NoKW) where T = NoMod{T}(nothing)

Base.eltype(::NoMod{T}) where T = T
source_eltype(::NoMod{T}) where T = T

# Modifies array values by scale
struct Mod{T1,T2,Mi,S,O,F} <: AbstractModifications{T1}
     missingval::Mi
     scale::S
     offset::O
     coerce::F
     function Mod{T}(missingval, scale, offset, coerce) where T
         missingval = missingval isa Pair && missingval[1] === missingval[2] ? missingval[1] : missingval
         if isnokw(coerce) || isnothing(coerce)
             coerce = convert
         end
         vals = map(_nokw2nothing, (missingval, scale, offset))
         T1 = _resolve_mod_eltype(T, vals...)
         new{T1,T,map(typeof, vals)...,typeof(coerce)}(vals..., coerce)
     end
end

Base.eltype(::Mod{T1}) where T1 = T1
source_eltype(::Mod{<:Any,T2}) where T2 = T2

function _resolve_mod_eltype(::Type{T}, missingval, scale, offset) where T
    omv = _outer_missingval(missingval)
    T1 = isnothing(omv) ? T : promote_type(T, typeof(omv))
    T2 = isnothing(scale) ? T1 : promote_type(T1, typeof(scale))
    T3 = isnothing(offset) ? T2 : promote_type(T2, typeof(offset))
    return T3
end

missingval(m::Mod) = m.missingval
missingval(m::NoMod) = m.missingval

_inner_missingval(m::Mod) = _inner_missingval(m.missingval)
_inner_missingval(mv) = mv
_inner_missingval(mv::Pair) = mv[1]

_outer_missingval(m::AbstractModifications) = _outer_missingval(m.missingval)
_outer_missingval(mv::Pair) = mv[2]
_outer_missingval(mv) = mv

# A wrapper for disk arrays that modifieds values lazily:
# scale and offset or replacing missing values
struct ModifiedDiskArray{T,N,D,M} <: DiskArrays.AbstractDiskArray{T,N}
    data::D
    mod::M
end
function ModifiedDiskArray(
    data::D, mod::M
) where {D<:AbstractArray{<:Any,N},M<:AbstractModifications{T}} where {T,N}
    return ModifiedDiskArray{T,N,D,M}(data, mod)
end

_maybe_modify(A::AbstractArray, mod::AbstractModifications) = 
    ModifiedDiskArray(A, mod)
_maybe_modify(A::AbstractArray, ::Nothing) = A

filename(A::ModifiedDiskArray) = filename(parent(A))
missingval(A::ModifiedDiskArray) = missingval(A.mod)
_metadata(A::ModifiedDiskArray, args...) = _metadata(parent(A), args...)

Base.parent(A::ModifiedDiskArray) = A.data
Base.size(A::ModifiedDiskArray, args...) = size(parent(A), args...)

DiskArrays.haschunks(A::ModifiedDiskArray) = DiskArrays.haschunks(parent(A))
DiskArrays.eachchunk(A::ModifiedDiskArray) = DiskArrays.eachchunk(parent(A))

function DiskArrays.readblock!(
    A::ModifiedDiskArray{<:Any,<:Any,<:Any,<:Mod}, outer_block, I::AbstractVector...
)
    if isdisk(parent(A))
        inner_block = similar(outer_block, eltype(parent(A)))
        DiskArrays.readblock!(parent(A), inner_block, I...)
    elseif isdisk(parent(parent(A)))
        inner_block = similar(outer_block, eltype(parent(parent(A))))
        DiskArrays.readblock!(parent(parent(A)), inner_block, I...)
    else
        inner_block = view(parent(A), I...)
    end
    outer_block .= _applymod.(inner_block, (A.mod,))
    return outer_block
end
function DiskArrays.readblock!(
    A::ModifiedDiskArray{<:Any,<:Any,<:Any,<:NoMod}, out_block, I::AbstractVector...
)
    if isdisk(parent(A))
        DiskArrays.readblock!(parent(A), out_block, I...)
    else
        out_block .= view(parent(A), I...)
    end
end

function DiskArrays.writeblock!(
    A::ModifiedDiskArray{<:Any,<:Any,<:Any,<:Mod}, block, I::AbstractVector...
)
    if isdisk(parent(A))
        modblock = _invertmod.((Val{source_eltype(A.mod)}(),), block, (A.mod,))
        return DiskArrays.writeblock!(parent(A), modblock, I...) 
    else
        parent(A)[I...] = _invertmod.((Val{source_eltype(A.mod)}(),), block, (A.mod,))
    end
end
function DiskArrays.writeblock!(
    A::ModifiedDiskArray{<:Any,<:Any,<:Any,<:NoMod}, block, I::AbstractVector...
)
    if isdisk(parent(A))
        DiskArrays.writeblock!(parent(A), block, I...)
    else
        parent(A)[I...] = block
    end
end

Base.@assume_effects :foldable function _applymod(x, m::Mod)
    if _ismissing(x, _inner_missingval(m))
        _outer_missingval(m)
    else
        _scaleoffset(x, m)
    end
end
Base.@assume_effects :foldable _applymod(x, ::NoMod) = x

_ismissing(x, mv) = ismissing(x) || x === mv 
_ismissing(x, ::Nothing) = ismissing(x)

_scaleoffset(x, m::Mod) = _scaleoffset(x, m.scale, m.offset)
_scaleoffset(x, scale, offset) = x * scale + offset
_scaleoffset(x, ::Nothing, offset) = x + offset
_scaleoffset(x, scale, ::Nothing) = x * scale
_scaleoffset(x, ::Nothing, ::Nothing) = x

Base.@assume_effects :foldable function _invertmod(::Val{T}, x, m::Mod) where T
    tm = if !isnothing(m.missingval) && _ismissing(x, _outer_missingval(m))
        return _inner_missingval(m)
    else
        x
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

function _stack_mods(eltypes::Vector, metadata::Vector, missingval::AbstractVector;
    scaled::Bool, coerce
)
    map(eltypes, metadata, missingval) do T, md, mv
        scale, offset = get_scale(md, scaled)
        _mod(T, mv, scale, offset, coerce)
    end
end
function _stack_mods(eltypes::Vector, metadata::Vector, missingval::Pair;
    scaled::Bool, coerce
)
    map(eltypes, metadata) do T, md
        scale, offset = get_scale(md, scaled)
        _mod(T, missingval, scale, offset, coerce)
    end
end

function _mod(::Type{T}, metadata, missingval; scaled::Bool, coerce) where T
    scale, offset = get_scale(metadata, scaled)
    _mod(T, missingval, scale, offset, coerce)
end
function _mod(::Type{T}, missingval, scale, offset, coerce) where T
    if (isnokwornothing(missingval) || !(missingval isa Pair && !(isnothing(last(missingval))))) && 
        isnokwornothing(scale) && isnokwornothing(offset)
        return NoMod{T}(missingval)
    else
        return Mod{T}(missingval, scale, offset, coerce)
    end
end

get_scale(metadata::NoKW, scaled::Bool) = nothing, nothing
function get_scale(metadata, scaled::Bool)
    scale = scaled ? get(metadata, "scale", nothing) : nothing
    offset = scaled ? get(metadata, "offset", nothing) : nothing
    return scale, offset
end

_mod_eltype(::AbstractArray{T}, ::NoMod) where T = T
_mod_eltype(::AbstractArray, m::Mod{T}) where T = T

_mod_inverse_eltype(::AbstractArray{T}, ::NoMod) where T = T
_mod_inverse_eltype(::AbstractArray{T}, m::Mod) where T =
    Base.promote_op(_invertmod, typeof(m.coerce), T, typeof(m))

_write_missingval_pair(A, missingval::Pair; kw...) = missingval
function _write_missingval_pair(A, missingval; 
    verbose=true, eltype, metadata=metadata(A), required=false
)::Pair
    source_mv = Rasters.missingval(A)
    disk_mv = if isnothing(source_mv) || isnothing(missingval)
        if required
            source_mv = isnothing(source_mv) ? missing : source_mv
            _writeable_missing(eltype; verbose)
        else
            nothing
        end
    elseif isnokw(missingval) || ismissing(missingval)
        # See if there is a missing value in metadata
        md_mv = Rasters.missingval(metadata)
        if isnothing(md_mv)
            _writeable_missing(eltype; verbose)
        else
            md_mv
        end
    else
        missingval
    end

    return disk_mv => source_mv  
end

function _read_missingval_pair(var, metadata, missingval)
    if isnokw(missingval)
        mv = Rasters.missingval(var, metadata) 
        isnothing(mv) ? nothing => nothing : mv => missing
    elseif isnothing(missingval)
        nothing => nothing
    elseif missingval isa Pair
        # Pair: inner and outer missing values are manually defined
        missingval
    elseif missingval === Rasters.missingval
        # `missingval` func: detect missing value and keep it as-is
        mv = Rasters.missingval(var, metadata)
        mv => mv
    else
        # Otherwise: detect missing value and convert it to `missingval`
        Rasters.missingval(var, metadata) => missingval
    end
end