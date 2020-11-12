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
    n = length(x)
    Δx= zeros(n-1)
    Δy= zeros(n-1)

    @inbounds for i = 1:n-1
        Δx[i] = x[i+1] - x[i]
        Δy[i] = y[i+1] - y[i]
    end

    @inbounds for i in 1:n-1
        # reusing Δx[i] as storage for coefficients a[i]
        # to minimize memory usage
        Δx[i] = Δy[i]/Δx[i]
    end
    # y[i] = b[i]
    LinearSpline(Δx,y,x)
end

const spline1 = linearspline
const lspline = linearspline


function quadraticspline(x::Vector,y::Vector)
    n = length(x)
    Δx= zeros(n-1)
    Δy= zeros(n-1)
    δ = ones(n)  # main diagonal
    μ = zeros(n-1)    # lower diagonal
    r = zeros(n)
    #M = zeros(n)
    @inbounds for i in 1:n-1
        Δx[i] = x[i+1] - x[i] #h_i = (x_i+1 - x_i)
        Δy[i] = y[i+1] - y[i]
    end
    @inbounds for i in 2:n-1
        µ[i]=1.0
        r[i]=2*Δy[i]/Δx[i]
    end
    A=Bidiagonal(δ, µ, :L)
    M=A\r
    a = zeros(n-1)
    b = zeros(n-1)
    c = zeros(n-1)
    @inbounds for i in 1:n
        c[i]=y[i]
    end
end

const spline2 = quadraticspline
const qspline = quadraticspline


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
    e=length(Δx)
    A=Tridiagonal(µ[2:n],δ,λ[1:n-1])
    r = rbound(r,e,Δx,Δy,style...)
    M=A\r
    λ,µ = λµbound(λ,µ,e,Δx,style...)
    M = Mbound(M,e,style...)
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
    CubicSpline(a,b,c,d,x)
end
const spline3 = cubicspline
const cspline = cubicspline


function cubicsplineip(x::Vector,y::Vector,style...)
    if style[1] == "periodic"
        y[1] != y[end] ? error("Periodic spline error: different y values at endpoints!") : nothing
    elseif style[1] == "constrained"
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
        @inbounds for i in 2:n-1
            λ[i] = Δx[i] / (Δx[i-1] + Δx[i])
            µ[i] = 1.0 - λ[i]
            r[i] = (6.0 / (Δx[i-1] + Δx[i])) * (Δy[i] / Δx[i] - Δy[i-1] / Δx[i-1])
        end

        e=length(Δx)
        r = rbound(r,e,Δx,Δy,style...)
        λ,µ = λµbound(λ,µ,e,Δx,style...)
        A=Tridiagonal(µ[2:n],δ,λ[1:n-1])
        M=A\r
        M = Mbound(M,e,style...)
        #a = zeros(n-1)
        #b = zeros(n-1)
        c = zeros(n-1)
        #d = zeros(n-1)

        @inbounds for i in 1:n
            # a[i] reuses λ[i] to optimize memory and
            # number of allocations
            λ[i]= (M[i+1]-M[i])/(6*Δx[i]) # a[i]
            #b[i]= 0.5*M[i]

            # c[i] reuses µ[i] to optimize memory and
            # number of allocations
            µ[i]= Δy[i]/Δx[i] - Δx[i]/6.0 * (2.0* M[i]+M[i+1]) #c[i]
            #d[i]= y[i]
        end
        CubicSpline(λ,M.*0.5,µ,y,x)
    end
end

"""
`cubicspline(x,y)` falls back to the default boundary condition `"natural"`
"""
function cubicspline(x,y)
    cubicspline(x,y,"natural")
end

"""
`cubicspline(x,y,"constrained")`

Cubic spline interpolation based on C.J.C. Kruger, [`Constrained Cubic Spline Interpolation for Chemical Engineering Applications`](http://www.korf.co.uk/spline.pdf), 2002
"""
function constrainedspline(x,y)
    length(x)!=length(y) ? error("Vector length mismatch. Length(x): $(length(x)) != length(y): $(length(y))") : nothing
    length(x) <3 ? error("Too few points: $(length(x)). Spline needs at least 3 points") : nothing

    n = length(x)-1
    Δx= zeros(n)
    Δy= zeros(n)
    σ= zeros(n+1) #slope
    a = zeros(n+1)
    b = zeros(n+1)
    c = zeros(n+1)
    d = zeros(n+1)

    @inbounds for i in 1:n
        Δx[i] = x[i+1] - x[i] #h_i = (x_i+1 - x_i)
        Δy[i] = y[i+1] - y[i]
    end

    # f[i]      left of  xi,yi with value y=fi(xi)
    # f[i+1]    right of xi,yi with value y=fi+1(xi)
    # fi'(xi) = fi+1'(xi) = f'(xi)
    @inbounds for i in 2:n-1
        sc = Δy[i-1]*Δy[i]
        if sc > 0.0 #slope does not change sign
            σ[i] = 2.0/(Δx[i]/Δy[i] + Δx[i-1]/Δy[i-1])
            #σ[i] = (Δy[i]/Δx[i] + Δy[i-1]/Δx[i-1]) / 2.0
        elseif sc <= 0.0
            σ[i] = 0.0 #slope changes sign
        end
        println("σ[$i]=$(σ[i])")
        #σ[i] = 2.0 / ((x[i+1]-x[i]) / (y[i+1]-y[i]) + (x[i]-x[i-1]) / (y[i]-y[i-1]))
    end

    σ[1] = 3.0*Δy[1] / (2.0*Δx[1]) - σ[2]/2.0
    println("σ[1]=$(σ[1])")
    σ[n] = 3.0*Δy[n-1] / (2.0*Δx[n-1]) - σ[n-1]/2.0
    println("σ[n]=$(σ[n])")
    #σ[1] = 1.5*(y[2]-y[1]) / (x[2]-x[1]) - σ[2]/2.0
    #σ[end]=1.5*(y[end]-y[end-1]) / (x[end]-x[end-1]) - σ[end-1]/2.0
    
    @inbounds for i in 2:n
        κa =  -2.0*(σ[i]+2*σ[i-1])/Δx[i-1] + 6*Δy[i-1]/(Δx[i-1]*Δx[i-1])
        println("κa[$i]=$(κa[i])")
        κb =  2.0*(2.0*σ[i]+σ[i-1])/Δx[i-1] - 6*Δy[i-1]/(Δx[i-1]*Δx[i-1])
        println("κb[$i]=$(κb[i])")
        
        d[i]=(κb-κa)/(6*Δx[i-1])
        println("d[$i]=$(d[i])")
        c[i]=(x[i]*κa - x[i-1]*κb)/(2*Δx[i-1])
        println("c[$i]=$(c[i])")
        b[i]=(Δy[i-1]- c[i]*(x[i]*x[i] - x[i-1]*x[i-1]) - d[i]*(x[i]*x[i]*x[i]-x[i-1]*x[i-1]*x[i-1]))/Δx[i-1]
        println("b[$i]=$(b[i])")
        a[i] = y[i-1] - b[i]*x[i-1] - c[i]*x[i-1]*x[i-1] - d[i]*x[i-1]*x[i-1]*x[i-1]
        println("a[$i]=$(a[i])")
    end
    e=zeros(n+1)
    CubicSpline(d,c,b,a,x)
end

"""
`interp(Spl::Spline,dλ)`

Returns tuple of vectors of spline values `s` at locations `xs`,
as defined by the grid points at `x+n*dλ` (smoothed curve).

Returns: `xs, s`
"""
function interp(Spl::CubicSpline,dλ)
    #xs=collect(Spl.x[1]:dλ:Spl.x[end]) #,(x[end]-x[1])/dλ)  # spline x vector
    xs=collect(dλ)
    s=zeros(length(xs))         # spline spectral vector
    j=1
    @inbounds for i in 1 : length(xs)
        if j < length(Spl.x)
            while xs[i] > Spl.x[j+1]
                j+=1
            end
            delta = xs[i]-Spl.x[j]
            s[i]=Spl.a[j]*delta^3.0 + Spl.b[j]*delta^2.0 + Spl.c[j]*delta + Spl.d[j]
        end
    end
    return xs,s
end

function interp(Spl::QuadraticSpline,dλ)
    #xs=collect(Spl.x[1]:dλ:Spl.x[end]) #,(x[end]-x[1])/dλ)  # spline x vector
    xs=collect(dλ)
    s=zeros(length(xs))         # spline spectral vector
    j=1
    @inbounds for i in 1 : length(xs)
        if j < length(Spl.x)
            while xs[i] > Spl.x[j+1]
                j+=1
            end
            delta = xs[i]-Spl.x[j]
            s[i]=Spl.a[j]*delta^2.0 + Spl.b[j]*delta + Spl.c[j]
        elseif j >= length(Spl.x)
            delta = xs[i]-Spl.x[j]
            s[i]=Spl.a[j]*delta^2.0 + Spl.b[j]*delta + Spl.c[j]
        end
    end
    return xs,s
end

function interp(Spl::LinearSpline,λ)
    #xs=collect(Spl.x[1]:dλ:Spl.x[end]) #,(x[end]-x[1])/dλ)  # spline x vector
    xs=collect(λ)
    s=zeros(length(xs))         # spline spectral vector
    j=1
    @inbounds for i in 1 : length(xs)
        if j < length(Spl.x)
            while xs[i] > Spl.x[j+1]
                j+=1
            end
            delta = xs[i]-Spl.x[j]
            s[i]=Spl.a[j]*delta + Spl.b[j]
        end
    end
    return xs,s
end

function interptest(Spl::CubicSpline,λ)
    xs=collect(λ) #,(x[end]-x[1])/dλ)  # spline x vector
    s=zeros(length(xs))         # spline spectral vector
    j=1
    @inbounds for i in 1 : length(xs)
        if j < length(Spl.x)
            while xs[i] > Spl.x[j+1]
                j+=1
            end
            delta = xs[i]-Spl.x[j]
            s[i]=Spl.a[j]*delta + Spl.b[j]
        end
    end
    return xs,s
end


"""
`rbound(r,Δx,Δy,style...)`

Returns vector `r` according to boundary conditions defined by `style...`

Returns: `r`
"""
function rbound(r,e,Δx,Δy,style...)
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
function Mbound(M,e,style...)
    if style[1] == "natural"
        M[1]   = 0.0
        M[end] = 0.0
    elseif style[1] == "clamped"
        ds=style[2]
        de=style[3]
    elseif style[1] == "periodic"
        M[1] = M[end]
        #λ[end] = Δx[1]/(Δx[1]+Δx[e-1])
        #µ[end] = 1-Δx[1]/(Δx[1]+Δx[e])
    elseif style[1] == "parabolic"
        M[1] = M[2]
        M[end] = M[end-1]
    elseif style[1] == "cubic"
        M[1] = 2*M[2]-M[3]
        M[end] = 2*M[end-1] - M[end-2]
    end
    return M
end

function λµbound(λ,µ,e,Δx,style...)
    if style[1] == "natural"
        λ[1]   = 0.0
        µ[end] = 0.0
    elseif style[1] == "clamped"
        ds=style[2]
        de=style[3]
        λ[1]   = 1.0
        µ[end] = 1.0
    elseif style[1] == "periodic"
        #M[1]   = M[end]
        λ[end] = Δx[1]/(Δx[1]+Δx[e-1])
        µ[end] = 1-Δx[1]/(Δx[1]+Δx[e])
    end
    return λ,µ
end

"""
`slope(x,y,dλ,style...)`

Returns the pair of `λ, s'` (1st derivative of the spline function).
"""
function slope(Spl::Spline,dλ)
    xs=collect(Spl.x[1]:dλ:Spl.x[end])#,(x[end]-x[1])/dλ)  # spline x vector
    s=zeros(length(xs))         # spline spectral vector
    j=1
    @inbounds for i in 1 : length(xs)
        if j < length(Spl.x)
            while xs[i] > Spl.x[j+1]
                j+=1
            end
            delta = xs[i]-Spl.x[j]
            typeof(Spl) == CubicSpline ? s[i]= 3.0*Spl.a[j]*delta^2.0 + 2.0*Spl.b[j]*delta + Spl.c[j] :
            typeof(Spl) == QuadraticSpline ? s[i]= 2.0*Spl.a[j]*delta + Spl.b[j] :
            typeof(Spl) == LinearSpline ? s[i] = Spl.a[j] : nothing
        end
    end
    return xs,s
end

const deriv1 = slope

"""
`curvature(x,y,dλ,style...)`

Returns the 2nd derivative of the spline function (s'')
"""
function curvature(Spl::Spline,dλ)
    xs=collect(Spl.x[1]:dλ:Spl.x[end])#,(x[end]-x[1])/dλ)  # spline x vector
    s=zeros(length(xs))         # spline spectral vector
    j=1
    @inbounds for i in 1 : length(xs)
        if j < length(Spl.x)
            while xs[i] > Spl.x[j+1]
                j+=1
            end
            delta = xs[i]-Spl.x[j]
            typeof(Spl) == CubicSpline ? s[i]=6.0*Spl.a[j]*delta + 2.0*Spl.b[j] :
            typeof(Spl) == QuadraticSpline ? s[i] = 2.0*Spl.a[j] :
            typeof(Spl) == LinearSpline ? s[i] = 0.0 : nothing
        end
    end
    return xs,s
end

const deriv2 = curvature

"""
`curvrate(x,y,dλ,style...`

Returns the 3rd derivative of the spline function (s''').
"""
function curvrate(Spl::Spline,dλ)
    xs=collect(Spl.x[1]:dλ:Spl.x[end])#,(x[end]-x[1])/dλ)  # spline x vector
    s=zeros(length(xs))         # spline spectral vector
    j=1
    @inbounds for i in 1 : length(xs)
        if j < length(Spl.x)
            while xs[i] > Spl.x[j+1]
                j+=1
            end
            typeof(Spl) == CubicSpline ? s[i]=6.0*Spl.a[j] :
            typeof(Spl) == QuadraticSpline ? s[i] = 0.0 :
            typeof(Spl) == LinearSpline ? s[i] = 0.0 : nothing
        end
    end
    return xs,s
end

const deriv3 = curvrate


function extrap(Spl::Spline,x::Real,type::AbstractString)
    if type == "zero"
        e=0.0
    elseif type == "boundary"
        if typeof(Spl) == CubicSpline
            if x < x[1]
                e= interp(spl3,)[2][end]
            elseif x > x[end]
                e= Spl.d[end]
            end
        elseif typeof(Spl) == QuadraticSpline
            if x < x[1]
                e=Spl.c[1]
            elseif x > x[end]
                e=Spl.c[end]
            end
        elseif typeof(Spl) == LinearSpline
            if x < x[1]
                e=Spl.b[1]
            elseif x > x[end]
                e=Spl.b[end]
            end
            e
        end
    elseif type == "linear"
        if typeof(Spl) == CubicSpline
            e
        elseif typeof(Spl) == QuadraticSpline
            e
        elseif typeof(Spl) == LinearSpline
            e
        end
    elseif type == "parabolic"
        if typeof(Spl) == CubicSpline
            e
        elseif typeof(Spl) == QuadraticSpline
            e
        elseif typeof(Spl) == LinearSpline
            e
        end
    end
    return e
end