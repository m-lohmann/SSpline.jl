module SSpline

    using LinearAlgebra
    
    import Base: ==, +, -, *, /
    import Base: convert
    
    export  spline0, stepspline, nearestneighbor, nspline,
            spline1,linearspline, lspline,
            spline2, quadraticspline, qspline,
            spline3, cubicspline, cspline, cub,
            naturalspline, clampedspline, periodicspline,
            hermitespline, hspline,
            nakspline,
            interpolate, interp, extrapolate, extrap,
            deriv1, grad, slope,
            deriv2, curvature,
            deriv3, curvrate,
            CubicSpline, HermiteSpline, QuadraticSpline, LinearSpline, NearestNeighborSpline

    include("types.jl")
    include("splinefunctions.jl")
    include("interpolation.jl")
    include("extrapolation.jl")
    include("derivatives.jl")
end #module