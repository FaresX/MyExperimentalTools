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
using Statistics
using Smoothers
using Peaks
using Integrals
using DataInterpolations
using Unitful
using LinearAlgebra
using FiniteDifferences

export fastload, resize, xyzforheatmap, prompeaks, findpks2, findpks, loedata, hma, transversal
export intIbias, intIbiasmap, sortd, normalization, interpVs, alignx, aligny, flip, interpxyz, loadxy, loadxyz
export binstep, binstepmap, IV2dIdV, IV2dIdVmap, jBtojx

include("DataAnalysis.jl")
# include("RSJ.jl")
# include("NISjunction.jl")

end
