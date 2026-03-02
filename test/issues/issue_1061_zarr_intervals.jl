# Reproduction of GitHub Issue #1061
# "Saving and reopening in Zarr shifts interval from Start to Center"
#
# When saving a Raster with Intervals{Start} to Zarr and reopening,
# the intervals are shifted from Start to Center, and dimension values
# shift accordingly (e.g., 0.0 becomes 0.005).

using Test
using Rasters
using DimensionalData
using DimensionalData.Lookups
using GeoInterface
using ZarrDatasets  # triggers the RastersZarrDatasetsExt extension

@testset "Issue #1061: Zarr roundtrip preserves Intervals{Start}" begin
    mktempdir() do tmpdir
        # --- MWE from the issue: rasterize produces Intervals{Start} ---
        points = GeoInterface.Point.(rand(100), rand(100))
        ras = rasterize(count, points; res=0.01, to=Extent(X=(0, 1), Y=(0, 1)))

        @test locus(dims(ras, X)) == Lookups.Start()
        @test locus(dims(ras, Y)) == Lookups.Start()

        filename = joinpath(tmpdir, "mwe.zarr")
        write(filename, ras)
        ras_disk = Raster(filename)

        @test locus(dims(ras_disk, X)) == Lookups.Start()
        @test locus(dims(ras_disk, Y)) == Lookups.Start()
        @test collect(lookup(ras_disk, X)) ≈ collect(lookup(ras, X))
        @test collect(lookup(ras_disk, Y)) ≈ collect(lookup(ras, Y))
        @test size(ras_disk) == size(ras)
    end
end

@testset "Zarr roundtrip: combinatorial locus × dimension" begin
    loci = (Lookups.Start(), Lookups.Center(), Lookups.End())

    for x_locus in loci, y_locus in loci
        x_name = nameof(typeof(x_locus))
        y_name = nameof(typeof(y_locus))
        @testset "X=$x_name, Y=$y_name" begin
            mktempdir() do tmpdir
                # Build expected index values for each locus
                x_index = if x_locus isa Lookups.Start
                    0.0:0.1:0.9
                elseif x_locus isa Lookups.Center
                    0.05:0.1:0.95
                else  # End
                    0.1:0.1:1.0
                end
                y_index = if y_locus isa Lookups.Start
                    0.0:0.1:0.9
                elseif y_locus isa Lookups.Center
                    0.05:0.1:0.95
                else  # End
                    0.1:0.1:1.0
                end

                ras = Raster(
                    rand(10, 10);
                    dims=(
                        X(Sampled(x_index;
                            order=Lookups.ForwardOrdered(),
                            span=Lookups.Regular(0.1),
                            sampling=Lookups.Intervals(x_locus))),
                        Y(Sampled(y_index;
                            order=Lookups.ForwardOrdered(),
                            span=Lookups.Regular(0.1),
                            sampling=Lookups.Intervals(y_locus))),
                    ),
                )

                # Verify pre-conditions
                @test locus(dims(ras, X)) == x_locus
                @test locus(dims(ras, Y)) == y_locus

                # Roundtrip through Zarr
                filename = joinpath(tmpdir, "test.zarr")
                write(filename, ras)
                ras_disk = Raster(filename)

                # Locus must survive the roundtrip
                @test locus(dims(ras_disk, X)) == x_locus
                @test locus(dims(ras_disk, Y)) == y_locus

                # Dimension values must survive the roundtrip
                @test collect(lookup(ras_disk, X)) ≈ collect(lookup(ras, X))
                @test collect(lookup(ras_disk, Y)) ≈ collect(lookup(ras, Y))

                # Sampling must remain Intervals
                @test sampling(dims(ras_disk, X)) isa Lookups.Intervals
                @test sampling(dims(ras_disk, Y)) isa Lookups.Intervals

                # Data must survive
                @test size(ras_disk) == size(ras)
            end
        end
    end
end
