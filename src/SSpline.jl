module SSpline

    using LinearAlgebra#, OffsetArrays
     
    import Base: ==, +, -, *, /
    import Base: convert
    
    export linearspline, lspline, spline1,
           quadraticspline, qspline, spline2,
           cubicspline, cspline, spline3, cubicsplineip,
           constrainedspline,
           naturalspline, clampedspline, periodicspline,
           interp, interptest, extrap,
           slope, deriv1,
           curvature, deriv2,
           curvrate, deriv3
    
    include("types.jl")
    include("ssplinefunctions.jl")
end #module