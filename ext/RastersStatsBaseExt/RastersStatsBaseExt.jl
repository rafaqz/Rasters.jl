module RastersStatsBaseExt

using Rasters, StatsBase
using StatsBase.Random

const RA = Rasters

import Rasters: _True, _False, _booltype
import Rasters.DimensionalData as DD

include("sample.jl")

end # Module
