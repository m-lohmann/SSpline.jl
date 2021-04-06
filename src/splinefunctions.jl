"""
    `stepspline(x::Vector,y::Vector)`
    
    Nearest neighbor interpolation with the spline values in the center of intervals.

    result: `NearestNeighborSpline(a,x)`
"""
function nearestneighbor(x::Vector,y::Vector)
    NearestNeighborSpline(y,x)
end

stepspline(x,y) = nearestneighbor(x,y)
spline0(x,y) = nearestneighbor(x,y)


"""
    `linearspline(x::Vector,y::Vector`

    Linear spline interpolation.

    result: `LinearSpline(a,b,x)`
"""
function linearspline(x::Vector,y::Vector)
    n = length(x)
    Δx= zeros(n-1)
    Δy= zeros(n-1)
    a=zeros(n-1)
    b=zeros(n-1)
    xn=copy(x)

    @inbounds for i = 1:n-1
        Δx[i] = x[i+1] - x[i]
        Δy[i] = y[i+1] - y[i]
    end

    @inbounds for i in 1:n-1
        # reusing Δx[i] as storage for coefficients a[i]
        # to minimize memory usage
        #Δx[i] = Δy[i]/Δx[i] #a[i]
        a[i] = Δy[i]/Δx[i]
        #reuse Δy[i] as b[i]
        b[i] = y[i]
        #Δy[i] = y[i] #b[i]
    end
    #push!(Δx,Δx[end])
    #push!(Δy,Δy[end]+Δx[end]*(xn[end]-xn[end-1]))
    LinearSpline(a,b,xn)
end

const spline1 = linearspline
const lspline = linearspline


function quadraticspline(x::Vector,y::Vector)
    n = length(x)
    Δx= zeros(n-1)
    Δy= zeros(n-1)
    d = zeros(n-1)
    f = zeros(n-1)
    #M = zeros(n)
    @inbounds for i in 1:n-1
        Δx[i] = x[i+1] - x[i] #h_i = (x_i+1 - x_i)
        Δy[i] = y[i+1] - y[i]
        d[i]  = Δy[i]/Δx[i]
        f[i] = 3*(Δy[i+1]-Δy[i])/(Δx[i+1]+Δx[i])
    end

    a = zeros(n-1)
    b = zeros(n-1)
    c = zeros(n-1)
    @inbounds for i in 1:n-1
        #a[i] = (f[i+1]-f[i]) / (6*Δx[i])
        b[i] = Δy[i]/Δx[i]
        c[i]=y[i]
    end
    QuadraticSpline(a,b,c,x)
end

const spline2 = quadraticspline
const qspline = quadraticspline


"""
    cubicspline(x::Vector,y::Vector,style...)

  * alias: `spline3(x::Vector,y::Vector,style...)`, `cspline(x,y,style...)`

Creates `CubicSpline(a,b,c,d)` object, containing the 4 spline coefficient vectors for a cubic spline.

style... settings for boundary conditions:

`:natural`: zero curvature at the start and end points of the spline.

  * alias: `cubicspline(x,y)`, `naturalspline(x,y)`

`:clamped,ds,de`: predefined 1st derivative of start (ds) and end points (de).

  * alias: `clampedspline(x,y,ds,de)`

`:periodic`: periodic boundary conditions. `y` values at start and end points have to be equal,
conditions at start and end points are identical.

  * alias: `periodicspline(x,y)`

`:parabolic`: parabolic runout spline

`:cubic`: cubic runout spline

Example:

    julia> x=[0,1,2,3,4]; y=[-1,1,5,2,-3]; spline3(x,y,:natural)
    CubicSpline([1.0, -3.0, 2.0, 0.0], [0.0, 3.0, -6.0, 0.0], [1.0, 4.0, 1.0, -5.0], [-1.0, 1.0, 5.0, 2.0])
"""
function cubicspline(x::Vector,y::Vector,style...)
    if style[1] == :periodic
        y[1] != y[end] ? error("Periodic spline error: different y values at endpoints!") : nothing
    elseif style[1] == :constrained
        constrainedspline(x,y)
    else
        n = length(x)
        Δx= zeros(n-1)
        Δy= zeros(n-1)
        λ = zeros(n)    # upper diagonal    
        δ = ones(n)*2  # main diagonal
        μ = zeros(n)    # lower diagonal
        r = zeros(n)
        M = zeros(n)

        @inbounds for i in 1:n-1
            Δx[i] = x[i+1] - x[i] #h_i = (x_i+1 - x_i)
            Δy[i] = y[i+1] - y[i]
        end

        @inbounds for i in 2:n
            λ[i] = Δx[i] / (Δx[i-1] + Δx[i])
            µ[i] = 1.0 - λ[i]
            r[i] = (6.0 / (Δx[i-1] + Δx[i])) * (Δy[i] / Δx[i] - Δy[i-1] / Δx[i-1])
        end

        e=length(Δx)
        r = rbound(r,e,Δx,Δy,style...)
        λ,µ = λµbound(λ,µ,e,Δx,style...)
        A=Tridiagonal(µ[2:n],δ,λ[1:n-1])
        M=A\r
        M = Mbound(M,style...)
        a = zeros(n-1) #
        b = zeros(n-1) #
        c = zeros(n-1)
        d = zeros(n-1) #

        @inbounds for i in 1:n-1
            # a[i] reuses λ[i] to optimize memory and
            # number of allocations
            #λ[i]= (M[i+1]-M[i])/(6*Δx[i]) # a[i]
            a[i]= (M[i+1]-M[i])/(6*Δx[i]) #
            b[i]= 0.5*M[i] #

            # c[i] reuses µ[i] to optimize memory and
            # number of allocations
            #µ[i]= Δy[i]/Δx[i] - Δx[i]/6.0 * (2.0* M[i]+M[i+1]) #c[i]
            c[i]= Δy[i]/Δx[i] - Δx[i]/6.0 * (2.0* M[i]+M[i+1])
            d[i]= y[i] #
        end
        #CubicSpline(λ,M.*0.5,µ,y,x)
    end
    spl=CubicSpline(a,b,c,d,x)

    push!(d,y[end])
    push!(c,slope(spl,x[end])[2][1])
    push!(b,curvature(spl,x[end])[2][1])
    push!(a,curvrate(spl,x[end])[2][1])
    CubicSpline(a,b,c,d,x)
end

const spline3 = cubicspline
const cspline = cubicspline


"""
`cubicspline(x,y)` falls back to the default boundary condition `:natural`
"""
function cubicspline(x,y)
    cubicspline(x,y,:natural)
end


"""
`naturalspline(x::Vector,y::Vector)`

Equivalent to `cubicspline(x::Vector,y::Vector,:natural)`.

Returns `CubicSpline(a,b,c,d)` containing all spline coefficients for a natural spline.
"""
function naturalspline(x::Vector,y::Vector)
    cubicspline(x::Vector,y::Vector,:natural)
end


"""
`clampedspline(x::Vector,y::Vector,ds,de)`

Equivalent to `cubicspline(x::Vector,y::Vector,ds,de)`.

Returns `CubicSpline(a,b,c,d)` containing all spline coefficients for a clamped spline.
"""
function clampedspline(x::Vector,y::Vector,ds,de)
    cubicspline(x::Vector,y::Vector,:clamped,ds,de)
end


"""
`periodicspline(x::Vector,y::Vector)`

Equivalent to `cubicspline(x::Vector,y::Vector,:periodic)`.
`y` values of start and end points have to be equal.

Returns `CubicSpline(a,b,c,d)` containing all spline coefficients for a periodic spline.
"""
function periodicspline(x::Vector,y::Vector)
    cubicspline(x::Vector,y::Vector,:periodic)
end


"""
`rbound(r,Δx,Δy,style...)`

Returns vector `r` according to boundary conditions defined by `style...`

Returns: `r`
"""
function rbound(r,e,Δx,Δy,style...)
    if style[1] == :natural
        r[1]   = 0.0
        r[end] = 0.0
    elseif style[1] == :clamped
        ds     = style[2]
        de     = style[3]
        r[1]   = 6.0/Δx[1]*(Δy[1]/Δx[1]-ds) #ds = start derivative
        r[end]   = -6.0/Δx[e]*(Δy[e]/Δx[e]-de)   #de = end derivative
    elseif style[1] == :periodic
        r[end] = 6.0/(Δx[1]+Δx[e])*(Δy[1]/Δx[1] - Δy[e]/Δx[e])
        r[1]   = r[end]
    end
    return r
end

"""
`Mbound(M,style...)`

Returns: vector `M` (moments) according to boundary conditions defined by `style...`
"""
function Mbound(M,style...)
    if style[1] == :natural
        M[1]   = 0.0
        M[end] = 0.0
#   elseif style[1] == :clamped
#       ds=style[2]
#       de=style[3]
    elseif style[1] == :periodic
        M[1] = M[end]
        #λ[end] = Δx[1]/(Δx[1]+Δx[e-1])
        #µ[end] = 1-Δx[1]/(Δx[1]+Δx[e])
    elseif style[1] == :parabolic
        M[1] = M[2]
        M[end] = M[end-1]
    elseif style[1] == :cubic
        M[1] = 2*M[2]-M[3]
        M[end] = 2*M[end-1] - M[end-2]
    end
    return M
end

function λµbound(λ,µ,e,Δx,style...)
    if style[1] == :natural
        λ[1]   = 0.0
        µ[end] = 0.0
    elseif style[1] == :clamped
        λ[1]   = 1.0
        µ[end] = 1.0
    elseif style[1] == :periodic
        #M[1]   = M[end]
        λ[end] = Δx[1]/(Δx[1]+Δx[e-1])
        µ[end] = 1-Δx[1]/(Δx[1]+Δx[e])
    end
    return λ,µ
end