# Species distribution modelling workflow

This example shows a full Species distribution modelling workflow, from loading data, to cleaning it, to fitting an ensemble and generating predictions.

It uses GBIF and WorldClim data, which are common datasets in ecology. We'll load occurrences for the Mountain Pygmy Possum species using [GBIF2.jl](https://github.com/rafaqz/GBIF2.jl), an interface to the [Global Biodiversity Information Facility](https://www.gbif.org/), and extract environmental variables using BioClim data from [RasterDataSources.jl](https://github.com/EcoJulia/RasterDataSources.jl).

## Load Rasters, ArchGDAL, RasterDataSources and GBIF
The GBIF2 library is used to download occurrence data, RasterDataSources to conveniently access Bioclim data. ArchGDAL is necessary to load in the Bioclim data.

````@example gbif
using Rasters, GBIF2
using RasterDataSources, ArchGDAL
````

Load occurrences for the Mountain Pygmy Possum using GBIF.jl

````@example gbif
records = GBIF2.occurrence_search("Burramys parvus"; limit=300)
````

## Get Bioclimatic variables
Get BioClim layers and subset to south-east Australia.
The first time this is run, this will automatically download and save the files.

````@example gbif
A = RasterStack(WorldClim{BioClim}, (1, 3, 7, 12))
se_aus = A[X(138 .. 155), Y(-40 .. -25), Band(1)]
````
Plot BioClim predictors and scatter occurrence points on all subplots

````@example gbif
# The coordinates from the gbif table
coords = collect(skipmissing(records.geometry))

using CairoMakie
p = Rasters.rplot(se_aus);
for ax in p.content
    if ax isa Axis
        scatter!(ax, coords; alpha=0.5, marker='+', color=:black, markersize = 20)
    end
end
p
````

## Extract bioclim variables at occurrence points
Then extract predictor variables and write to CSV. Use the skipmissing keyword to exclude both missing coordinates and coordinates with missing values in the RasterStack.

````@example gbif
using CSV
presences = extract(se_aus, coords, skipmissing = true)
CSV.write("burramys_parvus_predictors.csv", presences)
````

Or convert them to a `DataFrame`:

````@example gbif
using DataFrames
df = DataFrame(presences)
df[1:5,:]
````

## Sample background points
Next, sample random background points in the Raster. Rasters has a StatsBase extension to make this very straightforward. The syntax and output of `Rasters.sample` is very similar to that of `extract`.

````@example gbif
using StatsBase
background = Rasters.sample(se_aus, 500, skipmissing = true)
````

## Fit a statistical ensemble
In this example, we will [SpeciesDistributionModels.jl](https://github.com/tiemvanderdeure/SpeciesDistributionModels.jl) to fit a statistical ensemble to the occurrence and background data.

First we need to load the models. SDM.jl integrates with MLJ - see the [model browser](https://juliaai.github.io/MLJ.jl/dev/model_browser/#Classification) for what models are available.

````@example gbif
import Maxnet: MaxnetBinaryClassifier
import MLJGLMInterface: LinearBinaryClassifier
# define the models in the ensemble
models = (
    maxnet = MaxnetBinaryClassifier(), 
    maxnet2 = MaxnetBinaryClassifier(features = "lq"),
    glm = LinearBinaryClassifier()
)
````

Next, format the data using `sdmdata`. To test how rigurous our models are, we will use 3-fold cross-validation.

````@example gbif
using SpeciesDistributionModels
const SDM = SpeciesDistributionModels
data = sdmdata(presences, background; resampler = CV(; nfolds = 3))
````

Now, fit the ensemble, passing the data object and the `NamedTuple` of models!

````@example gbif
ensemble = sdm(data, models)
````

Use SDM.jl's evaluate function to see how this ensemble performs.

````@example gbif
SDM.evaluate(ensemble)
````

Not too bad!

## Make predictions of climatic suitability
Use the ensemble to 

````@example gbif
suitability = SDM.predict(ensemble, se_aus, reducer = mean)
````

And let's see what that looks like

````@example gbif
plot(suitability, colorrange = (0,1))
````
