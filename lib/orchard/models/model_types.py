########################################################################
# model_types.py
#
# Contains basic information the models that can be run using the
# python binary (/path/to/orchard/bin/orchard)
########################################################################
from BeamSearch import BeamSearch
from StochasticBeamSearch import StochasticBeamSearch
from ModelData import ModelData

# model type string constants
bs = "bs"
sbs = "sbs"

# model types
MODEL_TYPES = {
    bs: BeamSearch,
    sbs: StochasticBeamSearch
}

# model data types
MODEL_DATA = {
    bs: ModelData,
    sbs: ModelData
}