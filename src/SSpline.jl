module SSpline

    using LinearAlgebra
    using GRUtils
    
    import Base: ==, +, -, *, /
    import Base: convert
    
    export spline1, spline2, spline3,
           linearspline, quadraticspline, cubicspline,
           naturalspline, clampedspline, periodicspline,
           cspl, slope, curvature, curvrate
    
    include("types.jl")
    include("ssplinefunctions.jl")
    include("spectralfunctions.jl")
    include("types.jl")
    include("spectralops.jl")
end
