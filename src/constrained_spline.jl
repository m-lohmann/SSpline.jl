"""
`cubicspline(x,y,:constrained)`

Cubic spline interpolation based on C.J.C. Kruger, [`Constrained Cubic Spline Interpolation for Chemical Engineering Applications`](https://pages.uoregon.edu/dgavin/software/spline.pdf), 2002
"""
function constrainedspline(x,y)
    length(x)!=length(y) ? error("Vector length mismatch. Length(x): $(length(x)) != length(y): $(length(y))") : nothing
    length(x) < 3 ? error("Too few points: $(length(x)). Cubic splines need at least 3 points") : nothing

    n = length(x)-1
    Δx= zeros(n+1)
    Δy= zeros(n+1)
    σ = zeros(n+1) #slope

    @inbounds for i in 1:n-1 # 1:n-1
        Δx[i] = x[i+1] - x[i] #h_i = (x_i+1 - x_i)
        Δy[i] = y[i+1] - y[i]
    end

    @inbounds for i in 2:n-1
        slope = Δy[i-1]*Δy[i]
        if slope > 0.0
			σ[i] = 2.0/(Δx[i+1]/Δy[i] + Δx[i-1]/Δy[i-1])
		elseif slope <= 0.0
			σ[i] = 0.0
		end
    end

    σ[1] = 3.0* Δy[1]  /(2.0* Δx[1])   - σ[2]/2.0
    σ[n] = 3.0* Δy[n-1]/(2.0* Δx[n-1]) - σ[n-1]/2.0

    a = zeros(n+1)
    b = zeros(n+1)
    c = zeros(n+1)
    d = zeros(n+1)

    @inbounds for i in 2:n
        #κa =  -2.0*(σ[i+1]+2*σ[i])/Δx[i] + 6*Δy[i]/Δx[i]^2 # f''[i-1]
        #κb =  2.0*(2.0*σ[i+1]+σ[i])/Δx[i] - 6*Δy[i]/Δx[i]^2 # f''[i]
        κa = -2.0 * (σ[i] + 2* σ[i-1]) / Δx[i-1] + 6 * Δy[i-1] / (Δx[i-1]*Δx[i-1])
        κb = 2.0* (2.0 * σ[i] + σ[i-1]) / Δx[i-1] - 6 * Δy[i-1] / (Δx[i-1] * Δx[i-1])
        a[i]=(κb-κa)/(6*Δx[i-1])
        b[i]=(x[i]*κa - x[i-1]*κb)/(2*Δx[i-1])
        c[i]=(Δy[i-1]- b[i]*(x[i]^2 - x[i-1]^2) - a[i]*(x[i]^3 - x[i-1]^3))/Δx[i-1]
        d[i] = y[i-1] - c[i]*x[i-1] - b[i]*x[i-1]^2 - a[i]*x[i-1]^3
    end
    CubicSpline(a,b,c,d,x)
end
