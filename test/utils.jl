using Rasters, Test
using Rasters: cleankeys, _chunks_to_tuple, StableRange
using Rasters.DiskArrays

@testset "cleankeys" begin
    @test cleankeys(["a", "b", "c"]) == (:a, :b, :c)
    @test cleankeys(("a", "b", "c")) == (:a, :b, :c)
end

@testset "_chunks_to_tuple" begin
    A = zeros(100, 100)
    template = DiskArrays.RechunkedDiskArray(A, DiskArrays.GridChunks(A, (32, 16, 1)))
    # Bool chunk selection is generated from template
    _chunks_to_tuple(template, (X(), Y(), Ti()), true) == (32, 16, 1)
    @test _chunks_to_tuple(template, (X(), Y()), true) == (32, 16)
    @test _chunks_to_tuple(template, (Y(), X()), true) == (32, 16)
    @test _chunks_to_tuple(template, (Y(), X()), false) == nothing

    # Extra dims error
    @test_throws ArgumentError _chunks_to_tuple(template, (Y(), X()), (16, 32, 1))
    # Missing dims filled with 1s
    @test _chunks_to_tuple(template, (Y(), X()), (16,)) == (16, 1)
    # Correct number works as is
    @test _chunks_to_tuple(template, (Y(), X()), (16, 32)) == (16, 32)

    # Named chunks
    @test _chunks_to_tuple(template, (X(), Y()), (X=16, Y=32)) == (16, 32)
    @test _chunks_to_tuple(template, (X(), Y()), (Y=16, X=32)) == (32, 16)
    
    @test _chunks_to_tuple(template, (X(), Y()), (Y(8), X(12))) == (12, 8)
    @test _chunks_to_tuple(template, (X(), Y()), (Y(8), X(12), Ti(5))) == (12, 8)
    @test _chunks_to_tuple(template, (X(), Y(), Ti()), (Y(8), X(12))) == (12, 8, 1)
    @test _chunks_to_tuple(template, (Ti(), X(), Y()), (Y(8), X(12))) == (1, 12, 8)
end

@testset "StableRange" begin
    # Real GMTED 7.5" GeoTIFF parameters — the case PR #1070 was reported against.
    start, st, len = -28493.166784412522, 60.02213698319374, 514
    r = StableRange(; start, step=st, length=len)

    # Step is exactly the input, first value is exactly the input, last
    # value matches the default range to within 1 ULP.
    @test step(r) === st
    @test r[1] === start
    @test length(r) === len
    @test r[end] ≈ start + (len - 1) * st

    # Slicing is bit-stable for every slice — the property that fails on
    # `range(; start, step, length)` (and on `anchored_range` for slices
    # that didn't cover the anchor).
    for slice in (1:100, 1:300, 57:514, 200:400, 100:514, 1:len, 50:200, len÷2:len)
        s = r[slice]
        @test all(s[k] === r[k + first(slice) - 1] for k in eachindex(s))
    end

    # Nested slicing keeps bit-stability — offsets just accumulate.
    s1 = r[57:400]
    s2 = s1[10:200]
    @test all(s2[k] === r[k + 65] for k in eachindex(s2))

    # `Contains`-style lookups agree across slices: searchsortedlast on a
    # slice returns the parent's answer minus the slice's offset.
    for slice in (57:400, 200:514, 100:300)
        s = r[slice]
        shift = first(slice) - 1
        for k in 1:5:length(s)
            v = s[k]                                # mid-cell
            @test searchsortedlast(s, v) == searchsortedlast(r, v) - shift
            @test searchsortedfirst(s, v) == searchsortedfirst(r, v) - shift
            # Slightly off-cell — exercise the rounding path.
            v2 = v + st * 0.3
            @test searchsortedlast(s, v2) == searchsortedlast(r, v2) - shift
        end
    end

    # Negative step (the GeoTIFF Y case).
    r_neg = StableRange(; start=4255944.5659391750, step=-60.02213698319374, length=515)
    @test step(r_neg) === -60.02213698319374
    @test r_neg[1] === 4255944.5659391750
    @test all(r_neg[57:515][k] === r_neg[k + 56] for k in 1:459)

    # Range that never crosses zero — used to fall back to `Base.range`,
    # which made slicing drift. Now bit-stable like every other slice.
    r_pos = StableRange(; start=100.0, step=0.1, length=100)
    @test step(r_pos) === 0.1
    @test r_pos[1] === 100.0
    for slice in (1:50, 25:75, 50:100, 1:100)
        s = r_pos[slice]
        @test all(s[k] === r_pos[k + first(slice) - 1] for k in eachindex(s))
    end

    # Degenerate inputs match the default range's values.
    @test collect(StableRange(; start=0.0, step=1.0, length=10)) == collect(range(; start=0.0, step=1.0, length=10))
    @test collect(StableRange(; start=1.0, step=0.0, length=10)) == collect(range(; start=1.0, step=0.0, length=10))
    @test collect(StableRange(; start=1.0, step=1.0, length=1)) == collect(range(; start=1.0, step=1.0, length=1))
end
