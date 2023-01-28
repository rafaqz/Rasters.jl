# Load occurrences for the Mountain Pygmy Possum using GBIF.jl
# ## Load
using Rasters, GBIF2
const RS = Rasters

records = GBIF2.occurrence_search("Burramys parvus"; limit=300)

# ## Extract coordinates

# Extract the longitude/latitude value to a `Vector` of points
# (a `Tuple` counts as a `(x, y)` point in GeoInterface.jl):

coords = [(r.decimalLongitude, r.decimalLatitude) for r in records if !ismissing(r.decimalLatitude)]

# ## Get layer / Band
# Get BioClim layers and subset to south-east Australia

A = RasterStack(WorldClim{BioClim}, (1, 3, 7, 12))
e_aus = A[X(138 .. 155), Y(-40 .. -25), RS.Band(1)]

# Plot BioClim predictors and scatter occurrence points on all subplots

#kw = (legend=:none, opacity=0.5, markershape=:cross, markercolor=:black)
#p = plot(se_aus)
#foreach(i -> scatter!(p, coords; subplot=i, kw...), 1:4)
#p

# Then extract predictor variables and write to CSV.
