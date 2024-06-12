#Copyright © 2024, UChicago Argonne, LLC

from .aldchem import ALDideal
from .dose.ideal import ZeroD, WellStirred, ParticlePlugFlow

_ideal_models = {
    'zeroD' : ZeroD,
    'wellstirred' : WellStirred,
    'fluidizedbed' : ParticlePlugFlow
}

def aldmodel(process, model_name, **kwargs):
    if isinstance(process, ALDideal):
        return _ideal_models[model_name](process, **kwargs)

