using Rasters, Test
using Rasters: cleankeys, _chunks_to_tuple, anchored_range
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

@testset "anchored_range" begin
    # Real GMTED 7.5" GeoTIFF parameters — the case PR #1070 was reported against.
    start, st, len = -28493.166784412522, 60.02213698319374, 514
    r = anchored_range(start, st, len)

    # Step is exactly the input, first value is exactly the input, last
    # value matches the default range to within 1 ULP (interior values can
    # be more accurate than the default range because of the mid-anchor).
    @test step(r) === st
    @test r[1] === start
    @test length(r) === len
    @test r[end] ≈ start + (len - 1) * st

    # Anchor at the zero-crossing, not at index 1.
    @test r.offset == 476

    # Slicing identity — the property that breaks with `range(; start, step, length)`.
    for slice in (1:100, 1:300, 57:514, 200:400, 100:514, 1:len, 50:200, len÷2:len)
        s = r[slice]
        @test all(s[k] === r[k + first(slice) - 1] for k in eachindex(s))
    end

    # Negative step (the GeoTIFF Y case).
    r_neg = anchored_range(4255944.5659391750, -60.02213698319374, 515)
    @test step(r_neg) === -60.02213698319374
    @test r_neg[1] === 4255944.5659391750
    @test all(r_neg[57:515][k] === r_neg[k + 56] for k in 1:459)

    # Range that never crosses zero — anchor clamps to an endpoint and slicing
    # identity only holds for slices touching that endpoint.
    r_pos = anchored_range(100.0, 0.1, 100)
    @test step(r_pos) === 0.1
    @test r_pos[1] === 100.0
    @test r_pos.offset == 1                  # clamped, no interior zero-crossing
    # Same as the default range in this degenerate case
    @test r_pos === range(; start=100.0, step=0.1, length=100)

    # Degenerate inputs fall back to the default range.
    @test anchored_range(0.0, 1.0, 10) === range(; start=0.0, step=1.0, length=10)
    @test anchored_range(1.0, 0.0, 10) === range(; start=1.0, step=0.0, length=10)
    @test anchored_range(1.0, 1.0, 1) === range(; start=1.0, step=1.0, length=1)
end
