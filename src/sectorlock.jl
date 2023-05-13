##############################
# SectorLocks

# These lets us lock part of the array so we can work in parallel
# without always blocking other threads.

# Usually most polygons cover a small subset of a raster.
mutable struct SectorLock
    lock::Threads.SpinLock
    sector::CartesianIndices{2,Tuple{UnitRange{Int},UnitRange{Int}}}
end
SectorLock() = SectorLock(Threads.SpinLock(), CartesianIndices((1:0, 1:0)))

Base.lock(sl::SectorLock) = lock(sl.lock)
Base.unlock(sl::SectorLock) = unlock(sl.lock)
Base.islocked(sl::SectorLock) = islocked(sl.lock)

# One lock for each thread - this is only for threads
struct SectorLocks
    seclocks::Vector{SectorLock}
    spinlock::Threads.SpinLock
end
SectorLocks(n::Int) = SectorLocks([SectorLock() for _ in 1:n], Threads.SpinLock())
SectorLocks() = SectorLocks(_nthreads())

Base.length(sl::SectorLocks) = length(sl.seclocks)
Base.getindex(sl::SectorLocks, i::Int) = sl.seclocks[i]
Base.iterate(sl::SectorLocks, args...) = iterate(sl.seclocks, args...)
Base.eachindex(sl::SectorLocks) = 1:length(sl)

Base.lock(sl::SectorLocks, sector::RasterStack) = Base.lock(sl, first(sector))
function Base.lock(sl::SectorLocks, A::Raster)
    o = otherdims(A, DEFAULT_POINT_ORDER)
    if length(o) > 0
        slice = first(eachslice(A; dims=o))
        return Base.lock(sl, parent(slice))
    else
        return Base.lock(sl, parent(A))
    end
end
# Base.lock(sl::SectorLocks, sector::SubArray) = Base.lock(sl, CartesianIndices(_unitranges(sector.indices...))) 
# Base.lock(sl::SectorLocks, sector::DiskArrays.SubDiskArray) = Base.lock(sl, sector.v)
Base.lock(sl::SectorLocks, sector::Tuple) = Base.lock(sl, CartesianIndices(sector)) 
Base.lock(sl::SectorLocks, sector::AbstractArray) = Base.lock(sl, CartesianIndices(sector))
function Base.lock(seclocks::SectorLocks, sector::CartesianIndices)
    idx = Threads.threadid()
    thread_lock = seclocks.seclocks[idx]
    thread_lock.sector = sector
    # Lock the main lock so no other sectors can be locked
    lock(seclocks.spinlock)
    # Check for any other lock that intersects our sector
    for i in eachindex(seclocks)
        # Ignore a lock from this thread
        i == idx && continue
        seclock = seclocks[i]
        # Most of the time the lock will not intersect
        if islocked(seclock) && _intersects(seclock, sector)
            # If it does, lock the sector and wait 
            # All other threads are also waiting, 
            # but if this doens't happen often its fine.
            lock(seclock)
            unlock(seclock)
        end
    end
    # Lock this sector
    lock(thread_lock)
    # Unlock the main spinlock so other sectors can be locked
    unlock(seclocks.spinlock)
    # Work is now done in the locked sector, knowing no other thread
    # can write to it. `unlock` is run when it is finished.
end
Base.unlock(sl::SectorLocks) = begin
    # println("unlocking thread", idx)
    unlock(sl[Threads.threadid()])
end

_intersects(seclock::SectorLock, sector) = _intersects(seclock.sector, sector) 
_intersects(sector1, sector2) = length(intersect(sector1, sector2)) > 0

_unitranges(x1, xs...) = (_unitranges(x1)..., _unitranges(xs...)...)
_unitranges(x1::UnitRange) = (x1,)
_unitranges(x1::AbstractRange{Int}) = (UnitRange(min(first(x1), last(x1)), max(first(x1), last(x1))),)
_unitranges(x1) = ()
