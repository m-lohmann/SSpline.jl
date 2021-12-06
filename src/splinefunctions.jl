"""
    nearestneighbor(x, y)
    
Nearest neighbor interpolation with the spline values in the center of intervals.

- aliases: `spline0`, `stepspline`
"""
function nearestneighbor(x::Vector, y::Vector)
    NearestNeighborSpline(y, x)
end

const spline0 = nearestneighbor
const stepspline = nearestneighbor
const nspline = nearestneighbor


"""
    linearspline(x, y)

Linear spline interpolation.

- aliases: `spline1`, `lspline`
"""
function linearspline(x::Vector, y::Vector)
    #n = length(x)
    b = Vector{Float64}(undef, length(x))

    @inbounds for i in 1:length(x)-1
        b[i] = (y[i+1] - y[i]) / (x[i+1] - x[i])
    end
    #yy = y[begin+1:end].-y[begin:end-1]
    #xx = x[begin+1:end].-x[begin:end-1]
    LinearSpline(b, y, x)
end

const spline1 = linearspline
const lspline = linearspline


"""
    quadraticspline(x::Vector, y::Vector, b1 = 0.0)

Quadratic spline interpolation.

- aliases: `spline2`, `qspline`

Optional argument `b1` sets the 1st derivative at the first knot. Default value = `0.0`.
"""
function quadraticspline(x::Vector, y::Vector, b1 = 0.0)
    n = length(x)
    b = Vector{Float64}(undef, n-1)
    c = Vector{Float64}(undef, n-1)
    b[1] = b1

    Δx = Vector{Float64}(undef, n-1)
    Δy = Vector{Float64}(undef, n-1)

    @inbounds for i in 1:n-1
        Δx[i] = x[i+1] - x[i]
        Δy[i] = y[i+1] - y[i]
    end

    @inbounds for i in 1:n-1
        b[i+1] = -b[i] + 2 * Δy[i] / Δx[i]
    end

    @inbounds for i in 1:n-1
         c[i] = (b[i+1] - b[i]) / (2 * Δx[i])
    end
    QuadraticSpline(c, b, y, x)
end

const spline2 = quadraticspline
const qspline = quadraticspline


"""
    cubicspline(x::Vector, y::Vector, style...)

- alias: `spline3`, `cspline`

Creates a cubic spline. Boundary conditions are defined by `style`.

`style...`: settings for boundary conditions.

- `:natural`: zero curvature at the start and end points of the spline.
    - aliases: `cubicspline(x, y)`, `naturalspline(x, y)`
- `:clamped, ds, de`: predefined 1st derivative of start (ds) and end points (de).
    - alias: `clampedspline(x, y, ds, de)`
- `:periodic`: periodic boundary conditions. `y` values at start and end points have to be equal,
conditions at start and end points are identical.
    - alias: `periodicspline(x, y)`
- `:parabolic`, `:notaknot`: parabolic runout spline, also known as not-a-knot spline
- `:cubic`: cubic runout spline
- `:constrained`: constrained cubic spline, see Kruger’s paper
- `:hermite`: hermite spline
    - `:kruger`: constrained cubic spline, see Kruger’s paper. No overshooting.
    - `:free`: hermite spline with free derivative conditions. No overshooting.
    - `:predef`: predefined derivatives at knots

# Examples:
```jldoctest
julia> x=[0,1,2,3,4]; y=[-1,1,5,2,-3]; spline3(x,y,:natural)
CubicSpline([1.0, -3.0, 2.0, 0.0], [0.0, 3.0, -6.0, 0.0], [1.0, 4.0, 1.0, -5.0], [-1.0, 1.0, 5.0, 2.0])
```
"""
function cubicspline(x::Vector, y::Vector, style...)
    if style[1] == :periodic
        y[1] != y[end] ? error("Periodic spline error: different y values at endpoints!") : nothing
    elseif style[1] == :constrained
        return hermitespline(x, y, :kruger)
    elseif style[1] == :hermite
        if length(style) == 1
            return hermitespline(x, y, :kruger)
        elseif length(style) == 2
            return hermitespline(x, y, style[2])
        else
            throw(DomainError(style[1], "Spline style does not exist."))
        end
    else
        n = length(x)
        #Δx= Vector{Float64}(undef, n-1)
        Δx = zeros(Float64, n-1)
        #Δy= Vector{Float64}(undef, n-1)
        Δy = zeros(Float64, n-1)
        #λ = Vector{Float64}(undef, n)    # upper diagonal
        λ = zeros(Float64, n)
        δ = ones(n) * 2 # main diagonal
        #μ = Vector{Float64}(undef, n)    # lower diagonal
        μ = zeros(Float64, n)
        #r = Vector{Float64}(undef, n)    # derivatives
        r = zeros(Float64, n)
        #M = Vector{Float64}(undef, n)
        #M = zeros(Float64, n)

        @inbounds for i in 1:n-1
            Δx[i] = x[i+1] - x[i]
            Δy[i] = y[i+1] - y[i]
        end
        if style == :notaknot
            A = zeros(n, n)
            #=@inbounds=# for i in 2:n-1
                A[i, i-1] = Δx[i-1]
                A[i, i] = 2.0 * (Δx[i-1] + Δx[i])
                A[i, i+1] = Δx[i+1]
                r[i] = 6.0 * (Δy[i] / Δx[i] - Δy[i-1] / Δx[i-1])
            end
            #A[begin, begin] = -Δx[begin + 2]
            #A[begin, begin + 1] = Δx[begin + 1] + Δx[begin + 2]
            #A[begin, begin + 2] = -Δx[begin + 1]
            A[1,1] = -(x[3] - x[2])
            A[1,2] = x[3] - x[1]
            A[1,3] = -(x[2] - x[1])
            #A[end, end - 2] = -Δx[end - 1]
            #A[end, end - 1] = Δx[end - 2] + Δx[end - 1]
            #A[end, end]   = Δx[end - 2]
            A[n-1, n-3] = -(x[n-1] - x[n-2])
            A[n-1, n-2] = x[n-1] - x[n-3]
            A[n-1, n-1] = -(x[n-2] - x[n-3])
            #A[n-1, n-1] = -(x[n-3] - x[n-2])
            #A[n-1, n-2] = x[n-3] - x[n-1]
            #A[n-1, n-3] = -(x[n-2] - x[n-1])
            #for i in 2:n-1
            #    r[i] = 6.0 / (Δx[i-1] + Δx[i]) * (Δy[i] / Δx[i] - Δy[i-1] / Δx[i-1])
            #    #r[i] = 6.0 * (Δy[i] / Δx[i] - Δy[i-1] / Δx[i-1])
            #end
        else
            #=@inbounds=# for i in 2:n-1
                λ[i] = Δx[i] / (Δx[i-1] + Δx[i])
                μ[i] = 1.0 - λ[i]
                r[i] = (6.0 / (Δx[i-1] + Δx[i])) * (Δy[i] / Δx[i] - Δy[i-1] / Δx[i-1])
            end
        end

        e = length(Δx)
        r = rbound(r, e, Δx, Δy, style...)
        if style ≠ :notaknot
            λ, μ = λµbound(λ, μ, e, Δx, style...)
            A = Tridiagonal(μ[2:n], δ, λ[1:n-1])
            M = A \ r
            M = Mbound(M,style...)
        else
            M = A \ r
        end
        b = zeros(Float64, n-1) # Vector{Float64}(undef, n-1)
        c = zeros(Float64, n-1) # Vector{Float64}(undef, n-1)
        d = zeros(Float64, n-1) # Vector{Float64}(undef, n-1)

        #=@inbounds=# for i in 1:n-1
            d[i] = (M[i+1] - M[i]) / (6 * Δx[i])
            c[i] = 0.5 * M[i] #
            b[i] = Δy[i] / Δx[i] - Δx[i] / 6.0 * (2.0 * M[i] + M[i+1])
        end
    end
    if style[1] in (:hermite, :constrained)
        nothing
    else
        return CubicSpline(d, c, b, y, x)
    end
end


const spline3 = cubicspline
const cspline = cubicspline


"""
    rbound(r,Δx,Δy,style...)

Returns vector `r` according to boundary conditions defined by `style...`
Helper function for `cubicspline`.

Returns: `r`
"""
function rbound(r, e, Δx, Δy, style...)
    if style[1] in (:natural, :notaknot)
        r[1]   = 0.0
        r[end] = 0.0
    elseif style[1] == :parabolic
        r[1] = r[2]
        r[end] = r[end-1]
    elseif style[1] == :clamped
        ds     = style[2]
        de     = style[3]
        r[1]   = 6.0 / Δx[1] * (Δy[1] / Δx[1] - ds)    #ds = start derivative
        r[end] = -6.0 / Δx[e] * (Δy[e] / Δx[e] - de)   #de = end derivative
    elseif style[1] == :periodic
        r[end] = 6.0 / (Δx[1] + Δx[e]) * (Δy[1] / Δx[1] - Δy[e] / Δx[e])
        r[1]   = r[end]
    end
    return r
end

"""
    Mbound(M,style...)

Returns vector `M` (moments) according to boundary conditions defined by `style...`

Helper function for `cubicspline`.
"""
function Mbound(M, style...)
    if style[1] == :natural
        M[1]   = 0.0
        M[end] = 0.0
    elseif style[1] == :periodic
        M[1] = M[end]
    elseif style[1] == :parabolic
        M[1] = M[2]
        M[end] = M[end-1]
    elseif style[1] == :notaknot
        #M[1] = M[2]
        #M[end] = M[end-1]
    end
    return M
end

function λµbound(λ,μ,e,Δx,style...)
    if style[1] == :natural
        λ[1]   = 0.0
        μ[end] = 0.0
    elseif style[1] == :clamped
        λ[1]   = 1.0
        μ[end] = 1.0
    elseif style[1] == :parabolic
        λ[1] = λ[2]
        μ[end] = λ[end-1]
    elseif style[1] == :periodic
        λ[end] = Δx[1] / (Δx[1] + Δx[e-1])
        μ[end] = 1 - Δx[1] / (Δx[1] + Δx[e])
    end
    return λ, µ
end

"""
    cubicspline(x, y)
Fall back function, setting the default boundary condition to `:natural`

- alias: `naturalspline`
"""
function cubicspline(x, y)
    cubicspline(x, y, :natural)
end


"""
    naturalspline(x::Vector, y::Vector)

Equivalent to `cubicspline(x::Vector, y::Vector, :natural)`.
"""
naturalspline(x, y) = cubicspline(x, y, :natural)


"""
    clampedspline(x::Vector, y::Vector, ds = 0.0, de = 0.0)

Equivalent to `cubicspline(x::Vector, y::Vector, ds, de)`.

Returns a clamped cubic spline.

Default arguments: 1st and last derivatives = `0.0`
"""
function clampedspline(x, y, ds = 0.0, de = 0.0)
    cubicspline(x, y, :clamped, ds, de)
end


"""
    periodicspline(x::Vector, y::Vector)

Equivalent to `cubicspline(x::Vector, y::Vector, :periodic)`.

`y` values of start and end points have to be equal.
"""
function periodicspline(x, y)
    cubicspline(x, y, :periodic)
end

notaknotspline(x, y) = cubicspline(x, y, :notaknot)


"""
    hermitespline(x::Vector, y::Vector, style::Symbol = :kruger, args...)

Creates cubic hermite interpolation spline with defined 1st derivatives at all knots.

Allowed `style`s:
- `:kruger`: Default value. Defines 1st derivatives using harmonic means to prevent overshoot at extremal knots.
- `:free`: Simple mean values for 1st derivatives.
- `:predef`: user defined 1st derivatives at all knots.
"""
function hermitespline(x::Vector, y::Vector, style::Symbol = :kruger, args...)
    n = length(x)
    # 1st derivatives
    b = Vector{Float64}(undef, n)
    # finite differences
    Δx = Vector{Float64}(undef, n-1)
    Δy = Vector{Float64}(undef, n-1)

    @inbounds for i in 1:n-1
        Δx[i] = x[i+1] - x[i]
        Δy[i] = y[i+1] - y[i]
    end

    #generate 1st derivatives at kntos

    # calculate derivatives based on C.J.C. Kruger,
    # [`Constrained Cubic Spline Interpolation for Chemical Engineering Applications`](https://pages.uoregon.edu/dgavin/software/spline.pdf), 2002
    # which utilizes the harmonic mean of the neighboring slopes of each knot.
    if style == :kruger
        @inbounds for i in 2:n-1
            if Δy[i-1] * Δy[i] > 0.0
                b[i] = 2.0 / (Δx[i] / Δy[i] + Δx[i-1] / Δy[i-1]) #harmonic mean
            end
        end
        # boundary knots
        b[1] = 1.5 * Δy[1]/Δx[1] - b[2] / 2
        b[n] = 1.5 * Δy[n-1]/Δx[n-1] - b[n-1] / 2

    # 1st derivatives at knots are arithmetic means of neighboring slopes
    elseif style == :free
        @inbounds for i in 2:n-1
            b[i] = Δy[i] / Δx[i]
        end
        # boundary knots
        b[1] = Δy[1] / Δx[1]
        b[n] = Δy[n-1] / Δx[n-1]
    elseif style == :predef
        if length(args) == 0
            throw(ArgumentError("Derivative vector missing."))
        end
        if length(args[1]) ≠ n-1
            throw(DimensionMismatch("x, y, and derivative vectors must have the same length"))
        end
        @inbounds for i in 1:length(args[1])
            b[i] = args[1][i]
        end
    else
        throw(DomainError(style, "Derivative constraint style does not exist!"))
    end

    HermiteSpline(b, y, x)
end

const hspline = hermitespline

    
#    for i in 1:n-1
#        Φ0 = 2 * Δx[i]^3 - 3 * Δx[i]^2 + 1
#        Φ1 = Δx[i]^3 - 2 * Δx[i]^2 + Δx[i]
#        Ψ0 = -2 * Δx[i]^3 + 3 * Δx[i]^3
#        Ψ1 = Δx[i]^3 - Δx[i]^2
#
#        Φ0_ = 6 * Δx[i]^2 - 6 * Δx[i]
#        Φ1_ = 3 * Δx[i]^2 - 4 * Δx[i] + 1
#        Ψ0_ = -6 * Δx[i]^2 + 6 * Δx[i]
#        Ψ1_ = 3 * Δx[i]^2 - 2 * Δx[i]
#    
#        Φ0__ = 12 * Δx[i] - 6
#        Φ1__ = 6 * Δx[i] - 4
#        Ψ0__ = -12 * Δx[i] + 6
#        Ψ1__ = 6 * Δx[i] - 2
#    
#        Φ0___ = 12
#        Φ1___ = 6
#        Ψ0___ = -12
#        Ψ1___ = 6
#        #xb = (xs[i] - spl.x[j]) / (spl.x[j+1] - spl.x[j])
#        c[i] = (y[i] * Φ0_ + Δx[i] * b[i] * Φ1_ + y[i+1] * Ψ0_ + Δx[i] * b[i+1] * Ψ1_) / Δx[i]
#        #d[i] = (y[i] * Φ0__ + Δx[i] * b[i] * Φ1__ + y[i+1] * Ψ0__ + Δx[i] * b[i+1] * Ψ1__) / Δx[i]
#    end

    #CubicSpline(d, c, b, y, x)



"""
    nakspline(x::Vector, y::Vector)

Creates not-a-knot spline.
"""
function nakspline(x::Vector, y::Vector)
    n = length(x)
    Δx = zeros(Float64, n)
    Δy = zeros(Float64, n)
    r = zeros(Float64, n-1)
    for i in 1:n-1
        Δx[i] = x[i+1] - x[i]
        Δy[i] = y[i+1] - y[i]
    end
    A = zeros(n-1, n-1)

    for i in 2:n-2
        A[i,i-1] = Δx[i]
        A[i,i]   = 2 * (Δx[i] + Δx[i+1])
        A[i,i+1] = Δx[i+1]
        r[i] = 3 * Δy[i+1] / Δx[i+1] - 3 * Δy[i] / Δx[i]
    end

    A[1,1] = Δx[2]
    A[1,2] = -(Δx[1] + Δx[2])
    A[1,3] = Δx[1]

    A[n-1, n-3] = Δx[n-1]
    A[n-1, n-2] = -(Δx[n-2] + Δx[n-1])
    A[n-1, n-1] = Δx[n-2]

    r[1] = r[2]
    r[n-1] = r[n-2]

    #return A,r

    M = A \ r

    b = zeros(Float64, n-1) # Vector{Float64}(undef, n-1)
    c = zeros(Float64, n-1) # Vector{Float64}(undef, n-1)
    d = zeros(Float64, n-1) # Vector{Float64}(undef, n-1)
    for i in 1:n-2
        #d[i] = (M[i+1] - M[i]) / (6 * Δx[i])
        d[i] = (M[i+1] - M[i]) / (3 * Δx[i])
        c[i] = 0.5 * M[i] #
        #b[i] = Δy[i] / Δx[i] - Δx[i] / 6.0 * (2.0 * M[i] + M[i+1])
        b[i] = Δy[i] / Δx[i] - (M[i+1] + 2 * M[i]) * Δx[i] / 3
    end
    return CubicSpline(d, c, b, y, x)
end
