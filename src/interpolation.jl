"""
`interp(Spl::Spline,lambdas)`

Returns tuple of vectors of spline values `s` at locations `xs`,
as defined by the grid points at `x+n*dλ` (smoothed curve).

Returns: `xs, s`
"""
function interp(Spl::CubicSpline,lambdas)
    #xs=collect(Spl.x[1]:dλ:Spl.x[end]) #,(x[end]-x[1])/dλ)  # spline x vector
    xs=collect(lambdas)
    s=zeros(length(xs))         # spline spectral vector
    j=1                         # count through knots
    @inbounds for i in 1 : length(xs)
        if j < length(Spl.x)
            while xs[i] ≥ Spl.x[j+1]
                j+=1
            end
            if xs[i] == Spl.x[end]
                s[i] = Spl.d[end]
            else
            delta = xs[i]-Spl.x[j]
            s[i] = Spl.a[j]*delta^3.0 + Spl.b[j]*delta^2.0 + Spl.c[j]*delta + Spl.d[j]
            end
        end
    end
    return xs,s
end


"""
`interp(Spl::QuadraticSpline,lambdas)`

Create interpolated vectors based on a quadratic spline

"""
function interp(Spl::QuadraticSpline,lambdas)
    #xs=collect(Spl.x[1]:dλ:Spl.x[end]) #,(x[end]-x[1])/dλ)  # spline x vector
    xs=collect(lambdas)
    s=zeros(length(xs))         # spline spectral vector
    j=1
    @inbounds for i in 1 : length(xs)
        if j < length(Spl.x)
            while xs[i] ≥ Spl.x[j+1]
                j+=1
            end
            delta = xs[i]-Spl.x[j]
            s[i]=Spl.a[j]*delta^2.0 + Spl.b[j]*delta + Spl.c[j]
        elseif j ≥ length(Spl.x)
            delta = xs[i]-Spl.x[j]
            s[i]=Spl.a[j]*delta^2.0 + Spl.b[j]*delta + Spl.c[j]
        end
    end
    return xs,s
end

function interp(Spl::LinearSpline,lambdas)
    #xs=collect(Spl.x[1]:dλ:Spl.x[end]) #,(x[end]-x[1])/dλ)  # spline x vector
    xs=collect(lambdas)
    s=zeros(length(xs))         # spline spectral vector
    j=1
    @inbounds for i in 1 : length(xs)
        if j < length(Spl.x)
            while xs[i] > Spl.x[j+1]
                j+=1
            end
            if xs[i] == Spl.x[end]
                s[i] = Spl.a[end]*(Spl.x[end]-Spl.x[end-1])+Spl.b[end]
            else
                delta = xs[i]-Spl.x[j]
                s[i]=Spl.a[j]*delta + Spl.b[j]
            end
        end
    end
    return xs,s
end

function interp(Spl::NearestNeighborSpline,lambdas)
    xs=collect(lambdas)
    s=zeros(length(xs))
    j=1
    @inbounds for i in 1:length(xs)
        if j < length(Spl.x)
            while xs[i] > Spl.x[j+1]
                j+=1
            end
            if xs[i] == Spl.x[end]
                s[i] = Spl.a[end]
            elseif 2 * xs[i] < Spl.x[j]+Spl.x[j+1]
                s[i]=Spl.a[j]
            else
                s[i]=Spl.a[j+1]
            end
        end
    end
    return xs,s
end