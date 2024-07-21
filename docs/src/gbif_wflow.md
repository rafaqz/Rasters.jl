Load occurrences for the Mountain Pygmy Possum using GBIF.jl

## Load GBIF

````@example gbif
using Rasters, GBIF2
using RasterDataSources
const RS = Rasters
````

````@example gbif
records = GBIF2.occurrence_search("Burramys parvus"; limit=300)
````

## Extract coordinates

Extract the longitude/latitude value to a `Vector` of points
(a `Tuple` counts as a `(x, y)` point in GeoInterface.jl):

````@example gbif
coords = [(r.decimalLongitude, r.decimalLatitude) for r in records if !ismissing(r.decimalLatitude)]
````

## Get layer / Band
Get BioClim layers and subset to south-east Australia

````@example gbif
A = RasterStack(WorldClim{BioClim}, (1, 3, 7, 12))
se_aus = A[X(138 .. 155), Y(-40 .. -25), RS.Band(1)]
````
Plot BioClim predictors and scatter occurrence points on all subplots

````@example gbif
using Plots
p = plot(se_aus);
kw = (legend=:none, opacity=0.5, markershape=:cross, markercolor=:black)
foreach(i -> scatter!(p, coords; subplot=i, kw...), 1:4)
p
````

Then extract predictor variables and write to CSV.
````@example gbif
using CSV
predictors = collect(extract(se_aus, coords))
CSV.write("burramys_parvus_predictors.csv", predictors)
````

Or convert them to a `DataFrame`:

````@example gbif
using DataFrames
df = DataFrame(predictors)
df[1:5,:]
````
