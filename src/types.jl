abstract type Spline end

struct LinearSpline <: Spline
    a
    b
end

struct QuadraticSpline <: Spline
    a
    b
    c
end

struct CubicSpline <: Spline
    a
    b
    c
    d
end