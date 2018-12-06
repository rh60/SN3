using LinearAlgebra
using SparseArrays
using Polynomials
using FastGaussQuadrature
using MATLAB

module dim1
    using Polynomials
    using FastGaussQuadrature
    include("../dim1/lagrange.jl")
end

include("quad.jl")
include("lagrange.jl")
include("mesh.jl")
include("bvp.jl")
