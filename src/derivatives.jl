"""
`slope(Spl::Spline,dλ)`

Returns the pair of `λ, s'` (1st derivative of the spline function).
"""
function slope(Spl::Spline,dλ)
    xs=collect(dλ)#,(x[end]-x[1])/dλ)  # spline x vector
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
            typeof(Spl) == LinearSpline ? s[i] = Spl.a[j] :
            typeof(Spl) == NearestNeighborSpline ? s[i] = 0.0 : nothing
        end
    end
    return xs,s
end

const deriv1 = slope

"""
`curvature(Spl::Spline,dλ)`

Returns the 2nd derivative of the spline function (s'')
"""
function curvature(Spl::Spline,dλ)
    xs=collect(dλ)#,(x[end]-x[1])/dλ)  # spline x vector
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
            typeof(Spl) == LinearSpline ? s[i] = 0.0 : 
            typeof(Spl) == NearestNeighborSpline ? s[i] = 0.0 : nothing
        end
    end
    return xs,s
end

const deriv2 = curvature

"""
`curvrate(Spl::Spline,dλ)`

Returns the 3rd derivative of the spline function (s''').
"""
function curvrate(Spl::Spline,dλ)
    xs=collect(dλ)#,(x[end]-x[1])/dλ)  # spline x vector
    s=zeros(length(xs))         # spline spectral vector
    j=1
    @inbounds for i in 1 : length(xs)
        if j < length(Spl.x)
            while xs[i] > Spl.x[j+1]
                j+=1
            end
            typeof(Spl) == CubicSpline ? s[i]=6.0*Spl.a[j] :
            typeof(Spl) == QuadraticSpline ? s[i] = 0.0 :
            typeof(Spl) == LinearSpline ? s[i] = 0.0 :
            typeof(Spl) == NearestNeighborSpline ? s[i] = 0.0 : nothing
        end
    end
    return xs,s
end

const deriv3 = curvrate