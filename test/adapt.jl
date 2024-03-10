using Rasters, Adapt, Test
using Rasters.GeoFormatTypes, Rasters.Lookups

struct CustomArray{T,N} <: AbstractArray{T,N}
    arr::Array
end

CustomArray(x::Array{T,N}) where {T,N} = CustomArray{T,N}(x)
Adapt.adapt_storage(::Type{<:CustomArray}, xs::Array) = CustomArray(xs)

Base.size(x::CustomArray, y...) = size(x.arr, y...)
Base.getindex(x::CustomArray, y...) = getindex(x.arr, y...)

@testset "adapt AbstractProjected" begin
    l = Projected([1:10...]; 
        crs=ProjString("+proj="), mappedcrs=EPSG(4326),
        span=Regular(1), metadata=Metadata(:a=>"1", :b=>"2")
    )
    l1 = Adapt.adapt(CustomArray, l)
    @test parent(parent(l1)) isa CustomArray
    @test parent(parent(l1)).arr == [1:10...]
    @test metadata(l1) == NoMetadata()
    l = Mapped([1:10...]; 
        crs=ProjString("+proj="), mappedcrs=EPSG(4326),
        span=Regular(1), metadata=Metadata(:a=>"1", :b=>"2")
    )
    l1 = Adapt.adapt(CustomArray, l)
    @test parent(parent(l1)) isa CustomArray
    @test parent(parent(l1)).arr == [1:10...]
    @test metadata(l1) == NoMetadata()
end

