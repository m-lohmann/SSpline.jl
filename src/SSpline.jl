module SSpline

    using LinearAlgebra#, OffsetArrays
     
    import Base: ==, +, -, *, /
    import Base: convert
    
    export stepspline, nearestneighbor, spline0,
           linearspline, lspline, spline1,
           quadraticspline, qspline, spline2,
           cubicspline, cspline, spline3, cubicsplineip,
           constrainedspline,
           naturalspline, clampedspline, periodicspline,
           interp, interptest, extrap,
           slope, deriv1,
           curvature, deriv2,
           curvrate, deriv3,
           extrapolate,
           CubicSpline, QuadraticSpline, LinearSpline
           #testconstrained
           #runtests
    
    #include("./../test/runtests.jl")
    include("types.jl")
    include("splinefunctions.jl")
    include("interpolation.jl")
    include("extrapolation.jl")
    include("derivatives.jl")
    include("constrained_spline.jl")
    #include("testconstrained.jl")
    #include(".\\test\runtests.jl")
end #module