"""
extrapolate(Spl::Spline,lambdas,env::SpecEnvironment)
"""
function extrapolate(Spl,xrange,extr::Symbol)
    envx=collect(xrange)
    xs=0.0  # start of extrapolation range
    xe=0.0  # end of exrapolation range
    envx[1] < Spl.x[1] ? xs=envx[1] : xs=Spl.x[1]
    envx[end] > Spl.x[end] ? xe=envx[end] : xe=Spl.x[end]
    a=copy(Spl.a)
    b=copy(Spl.b)
    x=copy(Spl.x)

    if typeof(Spl) == QuadraticSpline
        c=copy(Spl.c)
    elseif typeof(Spl) == CubicSpline
        c=copy(Spl.c)
        d=copy(Spl.d)
    end

    if extr in (:zero, :boundary)
        if typeof(Spl) == LinearSpline
            if xs < x[1]
                pushfirst!(a,0.0)
                pushfirst!(b,extr == :zero ? 0.0 : Spl.b[1])
                pushfirst!(x,xs)
            end
            if xe > x[end]
                #push!(a,0.0)
                a[end]= 0.0
                push!(a,0.0)
                b[end]= extr == :zero ? 0.0 : Spl.b[end]
                push!(b,b[end])
                #push!(x,nextfloat(Spl.x[end]))
                x[end] = nextfloat(Spl.x[end])
                push!(x,xe)
            end
            LinearSpline(a,b,x)

        elseif typeof(Spl) == QuadraticSpline
            if xs < x[1]
                pushfirst!(a,0.0)
                pushfirst!(b,0.0)
                pushfirst!(c,extr == :zero ? 0.0 : Spl.c[1])
                pushfirst!(x,xs)
            end
            if xe > x[end]
                a[end] = 0.0
                b[end] = 0.0
                c[end] = (extr == :zero) ? 0.0 : Spl.c[end]
                push!(a,0.0)
                push!(b,0.0)
                push!(c,(extr == :zero) ? 0.0 : Spl.c[end])
                x[end] = nextfloat(Spl.x[end])
                push!(x,xe)
            end
            QuadraticSpline(a,b,c,x)

        elseif typeof(Spl) == CubicSpline
            if xs < x[1]
                pushfirst!(a,0.0)
                pushfirst!(b,0.0)
                pushfirst!(c,0.0)
                pushfirst!(d,(extr == :zero) ? 0.0 : Spl.d[1])
                pushfirst!(x,xs)
            end
            if xe ≥ x[end]
                #a[end] = 0.0
                #b[end] = 0.0
                #c[end] = 0.0
                push!(a,0.0)
                push!(b,0.0)
                push!(c,0.0)
                #d[end] = (extr == :zero) ? 0.0 : Spl.d[end]
                push!(d,(extr == :zero) ? 0.0 : Spl.d[end])
                x[end] = nextfloat(x[end])
                push!(x,xe)
            end
            CubicSpline(a,b,c,d,x)
        end

    elseif extr == :linear
        if typeof(Spl) == LinearSpline
            if xs < x[1]
                pushfirst!(a,a[1])
                pushfirst!(b,b[1]-a[1]*(x[1]-xs))
                pushfirst!(x,xs)
            end
            if xe ≥ x[end]
                push!(a,a[end])
                push!(b,b[end]+a[end]*(xe-x[end]))
                push!(x,xe)
            end
            LinearSpline(a,b,x)
        elseif typeof(Spl) == QuadraticSpline
            if envx[1] < x[1]
                pushfirst!(a,0.0)
                pushfirst!(b,b[1])
                pushfirst!(c,c[1]-b[1]*(x[1]-xs))
                pushfirst!(x,xs)
            end
            if xe ≥ x[end]
                push!(a,0.0)
                push!(b,b[end])
                push!(c,c[end]+b[end]*(xe-x[end]))
                push!(x,xs)
            end
            QuadraticSpline(a,b,c,x)
        elseif typeof(Spl) == CubicSpline
            if xs < x[1]
                pushfirst!(a,0.0)
                pushfirst!(b,0.0)
                pushfirst!(c,c[1])
                pushfirst!(d,d[1]-c[1]*(x[1]-xs))
                pushfirst!(x,xs)
            end
            if xe ≥ x[end]
                #a[end] = 0.0
                #b[end] = 0.0
                push!(a,0.0)
                push!(b,0.0)
                dx=x[end]-x[end-1]
                push!(c,3*a[end]*dx^2+2*b[end]*dx+c[end])
                x[end] = nextfloat(x[end])

                push!(d,a[end-1]*dx^3+b[end-1]*dx^2+c[end-1]*dx+d[end-1])
                push!(x,xe)
            end
            CubicSpline(a,b,c,d,x)
        end
    elseif extr == :quadratic
        if typeof(Spl) == LinearSpline
            if xs < x[1]
                pushfirst!(a,x[1]*((b[3]-b[2])+x[2]*(b[1]-b[3])+x[3]*(b[2]-b[1]))/((x[1]-x[2])*(x[1]-x[3])*(x[2]-x[3])))
                pushfirst!(b,b[1]/(a[1]*(x[1]-xs)))
                pushfirst!(x,xs)
            end
            if xe ≥ x[end]
                push!(a,0.0)
                push!(b,b[end]/(a[end])*(x[end]-xe))
                push!(x,xe)
            else nothing
            end
            LinearSpline(a,b,x)
        elseif typeof(Spl) == QuadraticSpline
            if xs < x[1]
                pushfirst!(a,0.0)
                pushfirst!(b,0.0)
                pushfirst!(c,c[1]/(b[1]*(x[1]-xs)))
                pushfirst!(x,xs)
            end
            if xe ≥ x[end]
                push!(a,0.0)
                push!(b,0.0)
                push!(c,c[end]/(b[end])*(x[end]-xe))
                push!(x,xe)
            end
            QuadraticSpline(a,b,c,x)
        elseif typeof(Spl) == CubicSpline
            if xs < x[1]
                pushfirst!(a,0.0)
                pushfirst!(b,0.0)
                pushfirst!(c,0.0)
                pushfirst!(d,d[1]/(c[1]*(x[1]-xs)))
                pushfirst!(x,xs)
            end
            if xe ≥ x[end]
                push!(a,0.0)
                push!(b,0.0)
                push!(c,0.0)
                push!(d,d[end]/(c[end])*(x[end]-envx[end]))
                push!(x,xe)
            end
            CubicSpline(a,b,c,d,x)
        end
    end
end