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

# ------------ GENERIC MODULES -------------------------------------------------
import FFTW
import WAV

# ------------ GLOBAL VARIABLES ------------------------------------------------
module_path = splitdir(@__FILE__)[1]                # Path to this module

# ------------ HEADERS ---------------------------------------------------------
# Load headers
for module_name in ["observer", "wopwop", "wav"]#"postprocessing"]
    include("FLOWNoise_"*module_name*".jl")
end

end # END OF MODULE
