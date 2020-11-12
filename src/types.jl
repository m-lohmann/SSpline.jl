abstract type Spline end

struct LinearSpline <: Spline
    a::Vector{Float64}
    b::Vector{Float64}
    x::Vector{Float64}
end

struct QuadraticSpline <: Spline
    a::Vector{Float64}
    b::Vector{Float64}
    c::Vector{Float64}
    x::Vector{Float64}
end

struct CubicSpline <: Spline
    a::Vector{Float64}
    b::Vector{Float64}
    c::Vector{Float64}
    d::Vector{Float64}
    x::Vector{Float64}
end