using FastGaussQuadrature
using LinearAlgebra
using SparseArrays
using SuiteSparse
using Polynomials
using MATLAB

struct Quad
    Points::Vector{Float64}
    Weights::Vector{Float64}
    function Quad(nq::Int)
        x,w = gausslegendre(nq)
        new( (x.+1)/2, w/2 )
    end
end

@inline function intSum(values::Vector{Float64},Q::Quad)::Float64
    values⋅Q.Weights
end

@inline function intSum(f,Q::Quad)::Float64
    values = f.(Q.Points)
    values⋅Q.Weights
end


include("lagrange.jl")
include("chebyshev.jl")
include("mesh.jl")
