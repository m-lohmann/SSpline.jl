"""
    interp(spl::CubicSpline, xrange)

Returns tuple of vectors of spline values `s` at locations `xs`,
as defined by the grid points at `x+n*dλ` (smoothed curve).

# Examples:

for x=[0,10,30,50,70,90,100] and y=[30,130,150,150,170,220,320]

```jldoctest
julia> cs = spline3(x,y,:natural)
CubicSpline([-0.015813397129186602, 0.009126794258373205, -0.0006937799043062205, -0.0013516746411483252, 0.007350478468899521, -0.01305023923444976], [0.0, -0.4744019138755981, 0.07320574162679427, 0.03157894736842105, -0.049521531100478466, 0.3915071770334928], [11.58133971291866, 6.83732057416268, -1.1866028708133973, 0.9090909090909092, 0.5502392344497609, 7.389952153110048], [30.0, 130.0, 150.0, 150.0, 170.0, 220.0], [0.0, 10.0, 30.0, 50.0, 70.0, 90.0, 100.0])

julia> interp(cs, 30.0:0.2:320.0)
([30.0, 30.2, 30.4, 30.6, 30.8, 31.0, 31.2, 31.4, 31.6, 31.8  …  318.2, 318.4, 318.6, 318.8, 319.0, 319.2, 319.4, 319.6, 319.8, 320.0], [150.00000000000003, 149.76560210526316, 149.53702736842106, 149.31424248803827, 149.09721416267942, 148.88590909090908, 148.68029397129186, 148.48033550239234, 148.28600038277511, 148.09725531100477  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
```
"""
function interpolate(spl::CubicSpline, xrange)
    xs=collect(xrange)
    s=zeros(length(xs))         # spline spectral vector
    j=1                         # count through knots
    @inbounds for i in 1 : length(xs)
        if j < length(spl.x)
            while (xs[i] > spl.x[j+1])
                j+=1
            end
            delta = xs[i]-spl.x[j]
            delta2 = delta * delta
            s[i] = spl.d[j] * delta2 * delta + spl.c[j] * delta2 + spl.b[j] * delta + spl.a[j]
        end
    end
    return xs,s
end

"""
    interpolate(spl::HermiteSpline, xrange)

Create interpolated vectors based on a Hermite spline

# Examples:

for x=[0,10,30,50,70,90,100] and y=[30,130,150,150,170,220,320]

```jldoctest
julia> hs = hermitespline(x,y,:kruger)
HermiteSpline([14.909090909090908, 0.18181818181818182, 0.0, 0.0, 0.5714285714285714, 0.16, 14.92], [30.0, 130.0, 150.0, 150.0, 170.0, 220.0, 320.0], [0.0, 10.0, 30.0, 50.0, 70.0, 90.0, 100.0])

julia> interp(hs, 30.0:0.2:320.0)
([30.0, 30.2, 30.4, 30.6, 30.8, 31.0, 31.2, 31.4, 31.6, 31.8  …  318.2, 318.4, 318.6, 318.8, 319.0, 319.2, 319.4, 319.6, 319.8, 320.0], [150.0, 150.0, 150.00000000000003, 149.99999999999997, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
```
"""
function interpolate(spl::HermiteSpline, xrange)
    xs = collect(xrange)
    s = Vector{Float64}(undef, length(xs))
    j = 1
    @inbounds for i in 1 : length(xs)
        if j < length(spl.x)
            while (xs[i] > spl.x[j+1])
                j+=1
            end
            #scaling and translation factor
            xb = (xs[i] - spl.x[j]) / (spl.x[j+1] - spl.x[j])
            xb2 = xb * xb
            #basis functions for function values
            Φ0 = 2.0 * xb2 * xb - 3.0 * xb2 + 1.0
            Φ1 = xb2 * xb - 2.0 * xb2 + xb
            #basis functions for derivatives
            Ψ0 = -2.0 * xb2 * xb + 3.0 * xb2
            Ψ1 = xb2 * xb - xb2
            Δx = spl.x[j+1] - spl.x[j]

            s[i] = spl.a[j] * Φ0 + Δx * spl.b[j] * Φ1 + spl.a[j+1] * Ψ0 + Δx * spl.b[j+1] * Ψ1
        end
    end
    return xs, s
end


"""
    interpolate(spl::QuadraticSpline, xrange)

Create interpolated vectors based on a quadratic spline

# Examples:

for x=[0,10,30,50,70,90,100] and y=[30,130,150,150,170,220,320]

```jldoctest
julia> qs = spline2(x,y)
QuadraticSpline([1.0, -0.95, 0.9, -0.85, 0.925, -1.1, 0.0], [0.0, 20.0, -18.0, 18.0, -16.0, 21.0, -1.0], [30.0, 130.0, 150.0, 150.0, 170.0, 220.0, 0.0], [0.0, 10.0, 30.0, 50.0, 70.0, 90.0, 100.0])

julia> interp(qs, 30.0:0.2:320.0)
([30.0, 30.2, 30.4, 30.6, 30.8, 31.0, 31.2, 31.4, 31.6, 31.8  …  318.2, 318.4, 318.6, 318.8, 319.0, 319.2, 319.4, 319.6, 319.8, 320.0], [150.0, 146.436, 142.94400000000002, 139.52399999999997, 136.176, 132.9, 129.69600000000003, 126.56400000000002, 123.50399999999998, 120.51599999999999  …  2.2006465129835e-310, 2.2034094153224e-310, 2.2061740509752e-310, 2.20894041994186e-310, 2.2117085222224e-310, 2.21447835781683e-310, 2.2172499267252e-310, 2.2200232289474e-310, 2.22279826448355e-310, 2.22557503333353e-310])
```
"""
function interpolate(spl::QuadraticSpline, xrange)
    xs = collect(xrange)
    s = zeros(length(xs))         # spline spectral vector
    j = 1
    @inbounds for i in 1 : length(xs)
        if j < length(spl.x)
            while xs[i] > spl.x[j+1]
                j+=1
            end
            delta = xs[i] - spl.x[j]
            s[i] = spl.c[j] * delta * delta + spl.b[j] * delta + spl.a[j]
        end
    end
    return xs,s
end


"""
    interpolate(spl::LinearSpline, xrange)

Create interpolated vectors based on a linear spline

# Examples:

for x=[0,10,30,50,70,90,100] and y=[30,130,150,150,170,220,320]

```jldoctest
julia> ls = spline1(x,y)
LinearSpline([10.0, 1.0, 0.0, 1.0, 2.5, 10.0], [30.0, 130.0, 150.0, 150.0, 170.0, 220.0], [0.0, 10.0, 30.0, 50.0, 70.0, 90.0, 100.0])

julia> interp(ls, 30.0:0.2:320.0)
([30.0, 30.2, 30.4, 30.6, 30.8, 31.0, 31.2, 31.4, 31.6, 31.8  …  318.2, 318.4, 318.6, 318.8, 319.0, 319.2, 319.4, 319.6, 319.8, 320.0], [150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
```
"""
function interpolate(spl::LinearSpline, xrange)
    xs = collect(xrange)
    s = zeros(length(xs))         # spline spectral vector
    j = 1
    @inbounds for i in 1 : length(xs)
        if j < length(spl.x)
            while xs[i] > spl.x[j+1]
                j+=1
            end
            s[i] = spl.b[j] * (xs[i] - spl.x[j]) + spl.a[j]
        end
    end
    return xs,s
end


"""
    interpolate(spl::NearestNeighborSpline, xrange)

Create interpolated vectors based on a nearest neighbor spline.

# Examples:

for x=[0,10,30,50,70,90,100] and y=[30,130,150,150,170,220,320]

```jldoctest
julia> nns = spline0(x,y)
NearestNeighborSpline([30.0, 130.0, 150.0, 150.0, 170.0, 220.0, 320.0], [0.0, 10.0, 30.0, 50.0, 70.0, 90.0, 100.0])

julia> interp(nns, 30.0:0.2:320.0)
([30.0, 30.2, 30.4, 30.6, 30.8, 31.0, 31.2, 31.4, 31.6, 31.8  …  318.2, 318.4, 318.6, 318.8, 319.0, 319.2, 319.4, 319.6, 319.8, 320.0], [150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
```
"""
function interpolate(spl::NearestNeighborSpline, xrange)
    xs = collect(xrange)
    s = zeros(length(xs))
    j = 1
    @inbounds for i in 1:length(xs)
        if j < length(spl.x)
            while xs[i] > spl.x[j+1]
                j+=1
            end
            if xs[i] == spl.x[end]
                s[i] = spl.a[end]
            elseif 2 * xs[i] < spl.x[j] + spl.x[j+1]
                s[i]=spl.a[j]
            else
                s[i]=spl.a[j+1]
            end
        end
    end
    return xs, s
end

# alias
const interp = interpolate