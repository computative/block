include("block.jl")

using DelimitedFiles: readdlm
x = readdlm("../resources/data.txt")
print(block(x))