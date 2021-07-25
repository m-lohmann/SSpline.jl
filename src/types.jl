abstract type Spline end

struct NearestNeighborSpline <: Spline
    a::Vector{Float64}
    x::Vector{Float64}
end
struct LinearSpline <: Spline
    b::Vector{Float64}
    a::Vector{Float64}
    x::Vector{Float64}
end

struct QuadraticSpline <: Spline
    c::Vector{Float64}
    b::Vector{Float64}
    a::Vector{Float64}
    x::Vector{Float64}
end

struct CubicSpline <: Spline
    d::Vector{Float64}
    c::Vector{Float64}
    b::Vector{Float64}
    a::Vector{Float64}
    x::Vector{Float64}
end

struct HermiteSpline <: Spline
    b::Vector{Float64}
    a::Vector{Float64}
    x::Vector{Float64}
end