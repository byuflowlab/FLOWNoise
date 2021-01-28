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

# FLOWMath https://github.com/byuflowlab/FLOWMath.jl
import FLOWMath
math = FLOWMath


# ------------ GENERIC MODULES -------------------------------------------------
import FFTW
import WAV
import PyPlot.LaTeXStrings: @L_str
import PyPlot
import Dates
import DataFrames
import CSV
import LinearAlgebra: norm, dot, cross, I
import Statistics: mean
import Base: read
using Printf

# Plot Format
const plt = PyPlot
plt.rc("font", family="Times New Roman")            # Text font
plt.rc("mathtext", fontset="stix")                  # Math font
SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 16
plt.rc("font", size=SMALL_SIZE)          # controls default text sizes
plt.rc("axes", titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc("axes", labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc("xtick", labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc("ytick", labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc("legend", fontsize=SMALL_SIZE)    # legend fontsize
plt.rc("figure", titlesize=BIGGER_SIZE)  # fontsize of the figure title

# ------------ GLOBAL VARIABLES ------------------------------------------------
const module_path = splitdir(@__FILE__)[1]                # Path to this module

const def_datasets_psw = Dict()     # Stores PSW data in this dictionary by default
const def_datasets_bpm = Dict()     # Stores BPM data in this dictionary by default

# ------------ HEADERS ---------------------------------------------------------
# Load headers
for module_name in ["observer", "wopwop", "wav", "processing",
                    "timedomain", "freqdomain", "postprocessing"]
    include("FLOWNoise_"*module_name*".jl")
end


function read(io, T, n)
    return [read(io, T) for i in 1:n]
end

end # END OF MODULE
