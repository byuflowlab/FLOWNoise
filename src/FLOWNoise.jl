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



# ---------- 1/3 Octave Bands --------------------
const fc_oto1 = [1, 1.25, 1.6, 2.0, 2.5, 3.15, 4, 5, 6.3, 8]     # (Hz) band center frequencies
const fc_oto = vcat(fc_oto1, fc_oto1*10, fc_oto1*100, fc_oto1*1000, fc_oto1*10000, 1e5)

const flow_oto1 = [0.88, 1.13, 1.414, 1.76, 2.25, 2.825, 3.53, 4.4, 5.65, 7.07] # (Hz) lower bound of 1/3-oct bins
const flow_oto = vcat(flow_oto1, flow_oto1*10, flow_oto1*100, flow_oto1*1000, flow_oto1*10000, 0.88e5)

const fhi_oto1 = [1.13, 1.414, 1.76, 2.25, 2.825, 3.53, 4.4, 5.65, 7.07, 8.8] # (Hz) upper bound of 1/3-oct bins
const fhi_oto = vcat(fhi_oto1, fhi_oto1*10, fhi_oto1*100, fhi_oto1*1000, fhi_oto1*10000, 1.13e5)

end # END OF MODULE
