module GeometryOpsDimensionalDataExt

import DimensionalData as DD
import GeometryOps as GO
import GeoInterface as GI

function GO.polygonize(A::DD.AbstractDimArray; dims=(DD.X(), DD.Y()), crs=GI.crs(A), kw...)
    lookups = DD.lookup(A, dims)
    bounds_vecs = if DD.isintervals(lookups)
        map(DD.intervalbounds, lookups)
    else
        @warn "`polygonsize` is not possible for `Points` sampling, as polygons cover space by definition. Treating as `Intervals`, but this may not be appropriate"
        map(lookups) do l
            Dd.intervalbounds(DD.set(l, DD.Intervals()))
        end
    end
    GO.polygonize(bounds_vecs..., DD.AbstractDimArray; crs, kw...)
end

end
