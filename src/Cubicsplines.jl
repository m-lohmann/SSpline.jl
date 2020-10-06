using LinearAlgebra, OffsetArrays, GRUtils#, GR

abstract type Spline end

struct CubicSpline <: Spline
    a
    b
    c
    d
end

function cubicspline(x::Vector,y::Vector,style...)
    n1=length(x)
    n  =n1-1

    x = OffsetArray(x,0:n)
    y = OffsetArray(y,0:n)
    Δx = OffsetArray(zeros(n1),0:n)
    Δy = OffsetArray(zeros(n1),0:n)
    λ  = OffsetArray(zeros(n1),0:n)    #upper diagonal
    d  = OffsetArray(ones(n1),0:n)*2   #main diagonal
    µ  = OffsetArray(zeros(n1),0:n)    #lower diagonal
    r  = OffsetArray(zeros(n1),0:n)
    M  = OffsetArray(zeros(n1),0:n)

    for i in 0:n-1
        Δx[i] = x[i+1] - x[i] #h_i = (x_i+1 - x_i)
        Δy[i] = y[i+1] - y[i]
    end


    for i in 1:n-1
        λ[i] = Δx[i] / (Δx[i-1] + Δx[i])
        µ[i] = 1 - λ[i]
        r[i] = (6 / (Δx[i-1] + Δx[i])) * (Δy[i] / Δx[i] - Δy[i-1] / Δx[i-1])
    end

    A=Tridiagonal(µ[1:end],collect(d),λ[0:end-1])
    r = boundary_r(r,Δx,Δy,style...)
    M=OffsetArray(A\collect(r),0:n)
    M,λ,µ = boundarycondition(M,λ,µ,Δx,style...)

    println("A=$A")
    println("M=$M")

    a = OffsetArray(zeros(n1),0:n)
    b = OffsetArray(zeros(n1),0:n)
    c = OffsetArray(zeros(n1),0:n)
    d = OffsetArray(zeros(n1),0:n)

    for i in 0:n-1
        a[i]= (M[i+1]-M[i])/(6*Δx[i])
        b[i]= 0.5*M[i]
        c[i]= Δy[i]/Δx[i] - Δx[i]/6.0 * (2.0* M[i]+M[i+1])
        d[i]= y[i]
    end

    println("a=$a")
    println("b=$b")
    println("c=$c")
    println("d=$d")
    result = CubicSpline(a,b,c,d)
    return result
end

function cspl(x,y,style::AbstractString...)
    dx = 0.1 * (x[end]-x[1])
    println("dx=$dx")
    cspl(x,y,dx,style...)
end

function cspl(x,y,dλ::Float64,style::AbstractString...)
    println("Computing Spl...")
    Spl=cubicspline(x,y,style...)
    println("cubicspline computed...")
    xs=collect(x[1]:dλ:x[end])#,(x[end]-x[1])/dλ)  # spline x vector
    s=zeros(length(xs))         # spline spectral vector
    j=1
    for i in 0 : length(xs)
        if j < length(x)
            if x[j] <= xs[i] && x[j+1] >= xs[i]
                #println("xs[$i]= $(xs[i]), x[$j] = $(x[j])")
                delta = xs[i]-x[j]
                s[i]=Spl.a[j]*delta^3 + Spl.b[j]*delta^2 + Spl.c[j]*delta + Spl.d[j]
            else
                j+=1
                delta = xs[i]-x[j]
                s[i]=Spl.a[j]*delta^3 + Spl.b[j]*delta^2 + Spl.c[j]*delta + Spl.d[j]
            end
        end
    end

    plot(x,y)
    hold(true)
    scatter(x,y)
    plot(xs,s)
    return xs,s
end

function boundary_r(r,Δx,Δy,style...)
    if style[1] == "natural"
        println("style=$(style[1])")
        r[0]   = 0.0
        r[end] = 0.0
    elseif style[1] == "clamp"
        println("style=$(style[1])")
        ds     = style[2]
        de     = style[3]
        println("derivatives:\nds=$ds\nde=$de")
        r[0]   = 6.0/Δx[1]*(Δy[1]/Δx[1]-ds) #ds = start derivative
        r[end] = 6.0/Δx[end]*(Δy[end]/Δx[end]-de)   #de = end derivative
    elseif style[1] == "periodic"
        println("style=$(style[1])")
        r[end] = 6.0/(Δx[1]+Δx[end])*(Δy[1]/Δx[1] - Δy[end]/Δx[end])
    end
    println("Boundary conditions for r computed...")
    return r
end

function boundarycondition(M,λ,µ,Δx,style...)
    if style[1] == "natural"
        M[0]   = 0
        M[end] = 0
        λ[0]   = 0
        µ[end] = 0
    elseif style[1] == "clamped"
        ds=style[2]
        de=style[3]
        λ[0]   = 1
        µ[0]   = 1
    elseif style[1] == "periodic"
        M[0]   = M[end]
        λ[end] = Δx[0]/(Δx[0]+Δx[end])
        µ[end] = 1-λ[end]
    end
    return M,λ,µ
end