"""
	extrapolate(spl, xrange, extr::Symbol)

Extrapolate a spline to the limits given by `xrange`, using the extrapolation type `extr`.

- alias: `extrap`

Allowed modes for `extr`:
- `:zero`
- `:boundary`
- `:linear`
- `:quadratic`
- `:smooth` (smooth blending to zero for Hermite splines)
"""
function extrapolate(spl, xrange, extr::Symbol)
    envx = Float64.(collect(xrange))
    xs = 0.0 # start of extrapolatin range
    xe = 0.0 # end of extrapolation range
    envx[1] < spl.x[1] ? xs = envx[1] : xs = spl.x[1]
    envx[end] > spl.x[end] ? xe = envx[end] : xe = spl.x[end]
    
    x = copy(spl.x)
    a = copy(spl.a)
	nf = nextfloat(Float64(x[end]))
    
    if typeof(spl) in (LinearSpline, QuadraticSpline, CubicSpline, HermiteSpline)
        b = copy(spl.b)
    end
    if typeof(spl) in (QuadraticSpline, CubicSpline)
        c = copy(spl.c)
    end
    if typeof(spl) == CubicSpline
        d = copy(spl.d)
    end
	if extr == :zero
		if xs < x[1]
			if typeof(spl) == HermiteSpline
				pushfirst!(a, 0.0)
				pushfirst!(b, 0.0)
				pushfirst!(x, prevfloat(x[1]))
			end
			x[1] = prevfloat(x[1])
			pushfirst!(x, xs)
			pushfirst!(a, 0.0)
			if typeof(spl) in (LinearSpline, QuadraticSpline, CubicSpline, HermiteSpline)
				pushfirst!(b, 0.0)
			end
			if typeof(spl) in (QuadraticSpline,CubicSpline)
				pushfirst!(c, 0.0)
			end
			if typeof(spl) == CubicSpline
				pushfirst!(d, 0.0)
			end
		end
		if xe ≥ x[end]
		if typeof(spl) == HermiteSpline
			push!(x, nf)
			push!(x, xe)
			push!(a, 0.0)
			push!(a, 0.0)
			push!(b, 0.0)
			push!(b, 0.0)
		end
			if typeof(spl) in (LinearSpline, QuadraticSpline, CubicSpline)
				x[end] = nf
				a[end] = 0.0
				push!(x, xe)
				push!(a, 0.0)
				push!(b, 0.0)
			end
			if typeof(spl) in (QuadraticSpline, CubicSpline)
				push!(c, 0.0)
			end
			if typeof(spl) == CubicSpline
				push!(d, 0.0)
			end
		end
	elseif extr == :boundary
		de = x[end] - x[end-1]
		if xs < x[1]
			if typeof(spl) == HermiteSpline
				pushfirst!(a, a[1])
				pushfirst!(b, 0.0)
				pushfirst!(x, prevfloat(x[1]))
			end
			pushfirst!(x, xs)
			pushfirst!(a, a[1])
			if typeof(spl) in (LinearSpline, QuadraticSpline, CubicSpline, HermiteSpline)
				pushfirst!(b, 0.0)
			end
			if typeof(spl) in (QuadraticSpline, CubicSpline)
				pushfirst!(c, 0.0)
			end
			if typeof(spl) == CubicSpline
				pushfirst!(d, 0.0)
			end
		end
		if xe ≥ x[end]
			if typeof(spl) == LinearSpline
				push!(a, a[end])
				push!(b, 0.0)
			end
			if typeof(spl) == QuadraticSpline
				push!(a, a[end])
				push!(b, 0.0)
				push!(c, 0.0)
				x[end] = nf
			end
			if typeof(spl) == CubicSpline
				de2 = de * de
				push!(a,a[end])
				push!(b,0.0)
				push!(c,0.0)
				push!(d,0.0)
				x[end] = nf
			end
			if typeof(spl) == HermiteSpline
				push!(a, a[end])
				push!(a, a[end])
				push!(x, nf)
				push!(b, 0.0)
				push!(b, 0.0)
			end
		end
		push!(x,xe)
	elseif extr == :linear
		if xs < x[1]
			pushfirst!(x, xs)
			pushfirst!(a, a[1] + b[1] * (x[1] - x[2]))
			if typeof(spl) in (LinearSpline, QuadraticSpline, CubicSpline, HermiteSpline)
				pushfirst!(b, b[1])
			end
			if typeof(spl) in (QuadraticSpline,CubicSpline)
				pushfirst!(c, 0.0)
			end
			if typeof(spl) == CubicSpline
				pushfirst!(d, 0.0)
			end
		end
		if xe ≥ x[end]
			dn = xe - x[end]
			dn2 = dn * dn
			de = x[end] - x[end-1]
			de2 = de * de
			if typeof(spl) == LinearSpline
				push!(a, a[end] + b[end] * dn)
				push!(b, b[end])
			end
			if typeof(spl) == QuadraticSpline
				push!(b, b[end] + c[end] * dn)
				push!(a, a[end] + b[end] * dn)
				push!(c, 0.0)
			end
			if typeof(spl) == CubicSpline
				push!(b, 3 * d[end] * de2 + 2 * c[end] * de + b[end])
				push!(a, a[end] + b[end] * dn)
				push!(c, 0.0)
				push!(d, 0.0)
				x[end] = nf
			end
			if typeof(spl) == HermiteSpline
				push!(a, a[end] + b[end] * dn)
				push!(b, b[end])
			end
		end
		push!(x,xe)
    elseif extr == :quadratic
        if xs < x[1]
            pushfirst!(x, xs)
            if typeof(spl) in (LinearSpline, QuadraticSpline, CubicSpline)
                pushfirst!(b, x[1] * ((a[3] - a[2]) + x[2] * (a[1] - b[3]) + x[3] * (a[2] - a[1])) /
							((x[1] - x[2]) * (x[1] - x[3]) * (x[2] - x[3])))
                pushfirst!(a, a[1] / (b[1] * (x[1] - xs)))
            end
            if typeof(spl) in (QuadraticSpline, CubicSpline)
                pushfirst!(c, 0.0)
            end
            if typeof(spl) == CubicSpline
                pushfirst!(d, 0.0)
            end
        end
        if xe ≥ x[end]
            push!(x, xe)
            if typeof(spl) in (LinearSpline, QuadraticSpline, CubicSpline)
                push!(b, x[end-2] * ((a[end] - a[end-1]) + x[end-1] * (a[end-2] - b[end]) + x[end] * (a[end-1] - a[end-2])) /
						((x[end-2] - x[end-1]) * (x[end-2] - x[end]) * (x[end-1] - x[end])))
                push!(a, a[end]/(b[end]*(x[end]-xs)))
            end
            if typeof(spl) in (QuadraticSpline, CubicSpline)
                push!(c, 0.0)
            end
            if typeof(spl) == CubicSpline
                push!(d, 0.0)
            end
        end
	elseif extr == :smooth
		pushfirst!(x, xs)
		pushfirst!(x, prevfloat(xs))
		pushfirst!(a, 0.0)
		pushfirst!(a, 0.0)
		push!(x, xe)
		push!(x, nextfloat(xe))
		push!(a, 0.0)
		push!(a, 0.0)
		smooth = hermitespline(x, a, :kruger)
		a = smooth.a
		b = smooth.b
    else
		throw(DomainError(extr, "Extrapolation type does not exist."))
    end
    if typeof(spl) == NearestNeighborSpline
        NearestNeighborSpline(a, x) 
    elseif typeof(spl) == LinearSpline
        LinearSpline(b, a, x)
    elseif typeof(spl) == QuadraticSpline
        QuadraticSpline(c, b, a, x)
	elseif typeof(spl) == CubicSpline
        CubicSpline(d, c, b, a, x)
	else
		HermiteSpline(b, a, x)
    end
end

"""
	extrapolate(spl, env = SpectralVis.SPECENV)

Extrapolate spline according to spectral environment.
"""
function extrapolate(spl, env = SpectralVis.SPECENV)
    extrapolate(spl, env.λmin:Δλ:λmax, env.ex)
end

# alias
extrap = extrapolate