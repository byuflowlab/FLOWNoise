"""
Tools for computational aeroacoustics (preprocessing, calculations, and
postprocessing) developed at BYU's FLOW lab.

    # AUTHORSHIP
      * Author    : Eduardo J. Alvarez and Tyler Critchfield
      * Email     : Edo.AlvarezR@gmail.com
      * Created   : Nov 2019
      * License   : MIT
"""
module FLOWNoise

# ------------ FLOW CODES ------------------------------------------------------
# GeometricTools https://github.com/byuflowlab/GeometricTools.jl
import GeometricTools
gt = GeometricTools

# # FLOWMath https://github.com/byuflowlab/FLOWMath.jl
# #! To make compatible with v0.6, comment out lines 6 and 7 in the module that include the roots.jl file and export the brent method
# import FLOWMath
# math = FLOWMath

# ------------ GENERIC MODULES -------------------------------------------------
import FFTW
import WAV
import PyPlot

const plt = PyPlot

# ------------ GLOBAL VARIABLES ------------------------------------------------
module_path = splitdir(@__FILE__)[1]                # Path to this module

# ------------ HEADERS ---------------------------------------------------------
# Load headers
for module_name in ["observer", "wopwop", "wav", "postprocessing"]
    include("FLOWNoise_"*module_name*".jl")
end

end # END OF MODULE
