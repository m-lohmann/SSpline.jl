"""
`naturalspline(x::Vector,y::Vector)`

Equivalent to `cubicspline(x::Vector,y::Vector,"natural")`.

Returns `CubicSpline(a,b,c,d)` containing all spline coefficients for a natural spline.
"""
function naturalspline(x::Vector,y::Vector)
    cubicspline(x::Vector,y::Vector,"natural")
end

"""
`clampedspline(x::Vector,y::Vector,ds,de)`

Equivalent to `cubicspline(x::Vector,y::Vector,ds,de)`.

Returns `CubicSpline(a,b,c,d)` containing all spline coefficients for a clamped spline.
"""
function clampedspline(x::Vector,y::Vector,ds,de)
    cubicspline(x::Vector,y::Vector,"clamped",ds,de)
end

"""
`periodicspline(x::Vector,y::Vector)`

Equivalent to `cubicspline(x::Vector,y::Vector,"periodic")`.
`y` values of start and end points have to be equal.

Returns `CubicSpline(a,b,c,d)` containing all spline coefficients for a periodic spline.
"""
function periodicspline(x::Vector,y::Vector)
    cubicspline(x::Vector,y::Vector,"periodic")
end

function linearspline(x::Vector,y::Vector)
end
const spline1 = linearspline

function quadraticspline(x::Vector,y::Vector)
end
const spline2 = quadraticspline


"""
    cubicspline(x::Vector,y::Vector,style...)

alias:
    
    spline3(x::Vector,y::Vector,style...)

Creates `CubicSpline(a,b,c,d)` object, containing the 4 spline coefficient vectors for a cubic spline.

style... settings for boundary conditions:

`"natural"`: zero curvature at the start and end points of the spline.

`"clamped",ds,de`: defined 1st derivative of start (ds) and end points (de). Hermite spline.

`"periodic"`: periodic boundary conditions. `y` values at start and end points have to be equal,
conditions at start and end points are identical.


Example:

    julia> x=[0,1,2,3,4]; y=[-1,1,5,2,-3]; spline3(x,y,"natural")
    CubicSpline([1.0, -3.0, 2.0, 0.0], [0.0, 3.0, -6.0, 0.0], [1.0, 4.0, 1.0, -5.0], [-1.0, 1.0, 5.0, 2.0])
"""
function cubicspline(x::Vector,y::Vector,style...)
    if style[1] == "periodic"
        y[1] != y[end] ? error("Periodic spline error: different y values at endpoints!") : nothing
    end
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
    @inbounds for i in 2:n-1
        λ[i] = Δx[i] / (Δx[i-1] + Δx[i])
        µ[i] = 1.0 - λ[i]
        r[i] = (6.0 / (Δx[i-1] + Δx[i])) * (Δy[i] / Δx[i] - Δy[i-1] / Δx[i-1])
    end

    A=Tridiagonal(µ[2:n],δ,λ[1:n-1])
    r = rbound(r,Δx,Δy,style...)
    M=A\r
    M,λ,µ = bound(M,λ,µ,Δx,style...)
    a = zeros(n-1)
    b = zeros(n-1)
    c = zeros(n-1)
    d = zeros(n-1)

    @inbounds for i in 1:n
        a[i]= (M[i+1]-M[i])/(6*Δx[i])
        b[i]= 0.5*M[i]
        c[i]= Δy[i]/Δx[i] - Δx[i]/6.0 * (2.0* M[i]+M[i+1])
        d[i]= y[i]
    end
    CubicSpline(a,b,c,d)
end
const spline3 = cubicspline


function cspl(x,y,style...)

end

"""
`cspl(x,y,dλ,style...)`

Returns tuple of vectors of spline values `s` at locations `xs`,
as defined by the grid points at `x+n*dλ` (smoothed curve).

Returns: `xs, s`
"""
function cspl(x,y,dλ,style...)
    Spl=cubicspline(x,y,style...)
    xs=collect(x[1]:dλ:x[end])#,(x[end]-x[1])/dλ)  # spline x vector
    s=zeros(length(xs))         # spline spectral vector
    j=1
    @inbounds for i in 1 : length(xs)
        if j < length(x)
            while xs[i] > x[j+1]
                j+=1
            end
            delta = xs[i]-x[j]
            s[i]=Spl.a[j]*delta^3.0 + Spl.b[j]*delta^2.0 + Spl.c[j]*delta + Spl.d[j]
        end
    end
    return xs,s
end

function plotspline(x::Vector,y::Vector)
    plot(x,y)
end

"""
`rbound(r,Δx,Δy,style...)`

Returns vector `r` according to boundary conditions defined by `style...`

Returns: `r`
"""
function rbound(r,Δx,Δy,style...)
    e=length(r)-1
    if style[1] == "natural"
        r[1]   = 0.0
        r[end] = 0.0
    elseif style[1] == "clamped"
        ds     = style[2]
        de     = style[3]
        r[1]   = 6.0/Δx[1]*(Δy[1]/Δx[1]-ds) #ds = start derivative
        r[end]   = -6.0/Δx[e]*(Δy[e]/Δx[e]-de)   #de = end derivative
    elseif style[1] == "periodic"
        r[end] = 6.0/(Δx[1]+Δx[e])*(Δy[1]/Δx[1] - Δy[e]/Δx[e])
        r[1]   = r[end]
    end
    return r
end

"""
`bound(M,λ,μ,Δx,style...)`

Returns tuple of vectors `M` (momenta), `λ` (upper diagonal), `μ` (lower diagonal)
according to boundary conditions defined by `style...`.

Returns: `M, λ, μ`
"""
function bound(M,λ,μ,Δx,style...)
    e=length(λ)-1
    if style[1] == "natural"
        M[1]   = 0.0
        M[end] = 0.0
        λ[1]   = 0.0
        µ[end] = 0.0
    elseif style[1] == "clamped"
        ds=style[2]
        de=style[3]
        λ[1]   = 1.0
        µ[1]   = 1.0
    elseif style[1] == "periodic"
        #M[1]   = M[end]
        λ[end] = Δx[1]/(Δx[1]+Δx[e-1])
        µ[end] = 1-Δx[1]/(Δx[1]+Δx[e])
    end
    return M,λ,µ
end

"""
`slope(x,y,dλ,style...)`

Returns the 1st derivative of the spline function.
"""
function slope(x,y,dλ,style...)
    Spl=cubicspline(x,y,style...)
    xs=collect(x[1]:dλ:x[end])#,(x[end]-x[1])/dλ)  # spline x vector
    s=zeros(length(xs))         # spline spectral vector
    j=1
    @inbounds for i in 1 : length(xs)
        if j < length(x)
            while xs[i] > x[j+1]
                j+=1
            end
            delta = xs[i]-x[j]
            s[i]=3.0*Spl.a[j]*delta^2.0 + 2.0*Spl.b[j]*delta + Spl.c[j]
        end
    end
    return xs,s
end

"""
`curvature(x,y,dλ,style...)`

Returns the 2nd derivative of the spline function (s'')
"""
function curvature(x,y,dλ,style...)
    Spl=cubicspline(x,y,style...)
    xs=collect(x[1]:dλ:x[end])#,(x[end]-x[1])/dλ)  # spline x vector
    s=zeros(length(xs))         # spline spectral vector
    j=1
    @inbounds for i in 1 : length(xs)
        if j < length(x)
            while xs[i] > x[j+1]
                j+=1
            end
            delta = xs[i]-x[j]
            s[i]=6.0*Spl.a[j]*delta + 2.0*Spl.b[j]
        end
    end
    return xs,s
end

"""
`curvrate(x,y,dλ,style...`

Returns the 3rd derivative of the spline function (s''').
"""
function curverate(x,y,dλ,style...)
    Spl=cubicspline(x,y,style...)
    xs=collect(x[1]:dλ:x[end])#,(x[end]-x[1])/dλ)  # spline x vector
    s=zeros(length(xs))         # spline spectral vector
    j=1
    @inbounds for i in 1 : length(xs)
        if j < length(x)
            while xs[i] > x[j+1]
                j+=1
            end
            delta = xs[i]-x[j]
            s[i]=6.0*Spl.a[j]
        end
    end
    return xs,s
end