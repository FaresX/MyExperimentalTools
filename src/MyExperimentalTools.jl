module MyExperimentalTools

using Reexport
# @reexport using DifferentialEquations
# @reexport using Plots
# @reexport using JLD2
# @reexport using Statistics
# @reexport using Smoothers
# @reexport using Peaks
# @reexport using Integrals
# @reexport using DataInterpolations
# @reexport using Unitful
# @reexport using LinearAlgebra
# @reexport using FiniteDifferences
# @reexport using FFTW
# @reexport using NativeFileDialog
# using Plots
using JLD2
using LsqFit
using FFTW
using Statistics
using Smoothers
using Peaks
using Roots
using Integrals
# using DataInterpolations
using Interpolations
using Unitful
using LinearAlgebra
using FiniteDifferences
using NumericalIntegration
using GLMakie
using Optim

export header, fastloaddf
export polyfitting, nonlinearfitting
export prompeaks, findpks, falldownpks!, falldownpks
export findIcs2, findIcs, filterpks!, filterpks, limitpksboundsy!, limitpksboundsy, falldownpks!, falldownpks, filterpksz!
export fastload, resize, xyzforheatmap, loedata, hma, transversal
export intIbias, sortd, normalization, interpVs, alignx, aligny, flip, interpxyz, loadxy, loadxyz, limitdata
export binstep, binstepmap, IV2dIdV, IV2dIdVmap, jBtojx
export symmetrizeRB, SdHdecomposition, analysisf, splitSdH, window
export selectbounds
export solveD, dR, interpmatrix
export lockin!, rcfilter!, lockinhys!

include("Fitting.jl")
include("FittingD.jl")
include("DataAnalysis.jl")
include("SdHAnalysis.jl")
include("CSVFiles.jl")
include("Peaks.jl")
include("Visualization.jl")
include("WalkHysteresis.jl")
# include("RSJ.jl")
# include("NISjunction.jl")

end
