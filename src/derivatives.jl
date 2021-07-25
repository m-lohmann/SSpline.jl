"""
    slope(spl::Spline,dλ)

Returns the pair of `λ, s'` (1st derivative of the spline function).
"""
function slope(spl::Spline,dλ)
    xs=collect(dλ)#,(x[end]-x[1])/dλ)  # spline x vector
    s=zeros(length(xs))         # spline spectral vector
    j=1
    @inbounds for i in 1 : length(xs)
        if j < length(spl.x)
            while xs[i] > spl.x[j+1]
                j+=1
            end
            delta = xs[i]-spl.x[j]
            if typeof(spl) == CubicSpline
                s[i]= 3.0 * spl.d[j] * delta * delta + 2.0 * spl.c[j] * delta + spl.b[j]
            elseif typeof(spl) == HermiteSpline
                xb = (xs[i] - spl.x[j]) / (spl.x[j+1] - spl.x[j])
                xb2 = xb * xb
                #Φ0_ = 6 * xb^2 - 6 * xb
                #Φ1_ = 3 * xb^2 - 4 * xb + 1
                #Ψ0_ = -6 * xb^2 + 6 * xb
                #Ψ1_ = 3 * xb^2 - 2 * xb
                h = spl.x[j+1] - spl.x[j]
                #s[i] = (spl.a[j] * Φ0_ + h * spl.b[j] * Φ1_ + spl.a[j+1] * Ψ0_ + h * spl.b[j+1] * Ψ1_) / h
                s[i] = (spl.a[j] * (6 * xb2 - 6 * xb) + h * spl.b[j] * (3 * xb2 - 4 * xb + 1) + spl.a[j+1] * (-6 * xb2 + 6 * xb) + h * spl.b[j+1] * (3 * xb2 - 2 * xb)) / h
            elseif typeof(spl) == QuadraticSpline
                s[i]= 2.0 * spl.c[j] * delta + spl.b[j]
            elseif typeof(spl) == LinearSpline
                s[i] = spl.b[j]
            elseif typeof(spl) == NearestNeighborSpline
                s[i] = 0.0
            end
        end
    end
    return xs,s
end

# alias
const deriv1 = slope

"""
    curvature(spl::Spline,dλ)

Returns the 2nd derivative of the spline function (s'')
"""
function curvature(spl::Spline,dλ)
    xs=collect(dλ)#,(x[end]-x[1])/dλ)  # spline x vector
    s=zeros(length(xs))         # spline spectral vector
    j=1
    @inbounds for i in 1 : length(xs)
        if j < length(spl.x)
            while xs[i] > spl.x[j+1]
                j+=1
            end
            delta = xs[i]-spl.x[j]
            if typeof(spl) == CubicSpline
                s[i] = 6.0 * spl.d[j] * delta + 2.0 * spl.c[j]
            elseif typeof(spl) == HermiteSpline
                xb = (xs[i] - spl.x[j]) / (spl.x[j+1] - spl.x[j])
                #Φ0__ = 12 * xb - 6
                #Φ1__ = 6 * xb - 4
                #Ψ0__ = -12 * xb + 6
                #Ψ1__ = 6 * xb - 2
                h = spl.x[j+1] - spl.x[j]
                #s[i] = (spl.a[j] * Φ0__ + h * spl.b[j] * Φ1__ + spl.a[j+1] * Ψ0__ + h * spl.b[j+1] * Ψ1__) / h
                s[i] = (spl.a[j] * (12 * xb - 6) + h * spl.b[j] * (6 * xb - 4) + spl.a[j+1] * (-12 * xb + 6) + h * spl.b[j+1] * (6 * xb - 2)) / h
                #s[i] = 2.0 * spl.b[j]
            elseif typeof(spl) == QuadraticSpline
                s[i] = 2.0 * spl.c[j]
            elseif typeof(spl) == LinearSpline
                s[i] = 0.0
            elseif typeof(spl) == NearestNeighborSpline
                s[i] = 0.0
            else
                nothing
            end
        end
    end
    return xs,s
end

# alias
const deriv2 = curvature

"""
    curvrate(spl::Spline,dλ)

Returns the 3rd derivative of the spline function (s''').
"""
function curvrate(spl::Spline,dλ)
    xs=collect(dλ)#,(x[end]-x[1])/dλ)  # spline x vector
    s=zeros(length(xs))         # spline spectral vector
    j=1
    @inbounds for i in 1 : length(xs)
        if j < length(spl.x)
            while xs[i] > spl.x[j+1]
                j+=1
            end
            typeof(spl) == CubicSpline ? s[i]= 6.0 * spl.d[j] :
            typeof(spl) == QuadraticSpline ? s[i] = 0.0 :
            typeof(spl) == LinearSpline ? s[i] = 0.0 :
            typeof(spl) == NearestNeighborSpline ? s[i] = 0.0 : nothing
        end
    end
    return xs,s
end

# alias
const deriv3 = curvrate