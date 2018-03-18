# BioStructures.jl
# ================
#
# A julia package to read, write and manipulate macromolecular structures.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioStructures.jl/blob/master/LICENSE.md

__precompile__()

module BioStructures

using Libz
using Formatting
using MMTF
import BioCore
import BioCore.distance
import BioSymbols
import BioSequences.AminoAcidSequence

include("model.jl")
include("pdb.jl")
include("mmcif.jl")
include("mmtf.jl")
include("spatial.jl")

end # BioStructures
