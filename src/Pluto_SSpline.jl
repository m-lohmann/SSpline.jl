### A Pluto.jl notebook ###
# v0.14.0

using Markdown
using InteractiveUtils

# ╔═╡ 90f87210-59c6-11eb-2d1c-8790aa130324
import Pkg

# ╔═╡ 9af1fa20-07fa-11eb-2d8d-239d4cc3f3ae
using SSpline

# ╔═╡ 2de193c2-080c-11eb-0451-5d5c46e11324
using Colors, Plots, OffsetArrays

# ╔═╡ 7b664047-b676-4a20-9c8c-5e610497f9fb
begin
	#import DarkMode
    #DarkMode.enable(theme="material-palenight", cm_config=Dict("tabSize" => 4))
	#DarkMode.enable()
    #DarkMode.Toolbox(theme="default")
end

# ╔═╡ 0c4701a0-dd1c-4f4b-91eb-cb94bc9d4440
md"cell width"

# ╔═╡ 679c9e8a-ef4f-4b2b-a1ce-245b0ffa773d
html"""<style>
main {
    max-width: 1000px;
}
"""

# ╔═╡ 3eb09f10-63f2-11eb-38c5-db7e4547b561
plotly()

# ╔═╡ 8890e470-63f0-11eb-3df9-19f6bddcf0fd
md"""
Create x and y vectors of original points:
"""

# ╔═╡ 26bca8a0-080c-11eb-33b2-7bf3e34cf10a
begin
	x=[0,10,30,50,70,90,100]
	y=[30,130,150,150,170,220,320]
	#x=[0,10,30,45,70,90,100]
	#y=[20,50,70,-70,60,50,20]
	#x=[0,10,30,50,70,80,100]
	#y=[150,200,200,200,180,100,0]
end

# ╔═╡ 9240d97e-63eb-11eb-0a53-278a141ee9a6
md"""Create nearest neighbor spline:

`nearestneighbor(x, y)`, `stepspline(x, y)`, `spline0(x,y)`

Result:

`NearestNeighborSpline(a, x)`
"""

# ╔═╡ 39f7e0b0-63ec-11eb-0a32-eb2d79a226d5
nns=nearestneighbor(x,y)

# ╔═╡ 2f3b8ef0-5920-11eb-3e40-39ee8b03512a
md"""Create linear spline:

`linearspline(x, y)`, `lspline(x, y)`, `spline1(x,y)`

Result:

`LinearSpline(a, b, x)`
"""

# ╔═╡ 6b7f17a0-59c6-11eb-2bd7-81340fa06b7c
ls=linearspline(x,y)

# ╔═╡ 028b3dc0-63f1-11eb-2b34-ab517656c8d3
md"""Create cubic splines:

* `cubicspline(x::Vector,y::Vector,style...)`
   * shorthands for `cubicspline(...)`:
     * `cspline(...)`
     * `spline3(...)`

`style...` =

* `:natural`
  * natural spline, zero curvature at both boundaries
    * shorthand: `cubicspline(x::Vector,y::Vector)`

* `:clamped, ds, de`
  * clamped spline, arbitrary slopes at start (`ds`) and end (`de`) boundaries
* `:periodic`
   * periodic spline, identical values, slope and curvature at both boundaries.
* `:constrained` (in the works)
   * constrained spline, no overshooting of the interpolation function. zero slope at spline values if at least one slope is zero or the slope has a a sign change.

Result:

object of type `CubicSpline` with the fields `a, b, c, d, x`
"""

# ╔═╡ a520a3e0-63f1-11eb-3d5c-31ad2e72d89c
md"spline interpolations:"

# ╔═╡ 55e3217a-e35e-496e-9f8d-2e042392bb68
function ocons(x,y)
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
        a[i-1]=(κb-κa)/(6*Δx[i-1])
        b[i-1]=(x[i]*κa - x[i-1]*κb)/(2*Δx[i-1])
        c[i-1]=(Δy[i-1]- b[i]*(x[i]^2 - x[i-1]^2) - a[i]*(x[i]^3 - x[i-1]^3))/Δx[i-1]
        d[i-1] = y[i-1] - c[i]*x[i-1] - b[i]*x[i-1]^2 - a[i]*x[i-1]^3
    end
    CubicSpline(a,b,c,d,x)
end

# ╔═╡ 017de170-5d17-11eb-2964-3b6996cb7eb8
begin
	csnat=cubicspline(x,y,:natural)
	csclamp=cubicspline(x,y,:clamped,-10,10)
	cstrain=ocons(x,y)
	#csconstr=conspline(x,y)
end

# ╔═╡ 6b1bc1ee-59c6-11eb-3a5c-2d0c5ab35d0c
begin
	nni=interp(nns,collect(0:1:100))
	lint=interp(ls,collect(0:100))
	cintn=interp(csnat,collect(0:100))
	cintc=interp(csclamp,collect(0:100))
	cints=interp(cstrain,collect(0:100))
end

# ╔═╡ 51f3c956-9462-490d-8375-a5001509a9d5
begin
	plot(nni[1],nni[2],color=:red, style=:dot, label="nearest neighbor", legend=:topleft)
	plot!(lint[1],lint[2], color=:orange1, label="linear")
	plot!(cintn[1],cintn[2],color=:green4,style=:dash, label="cubic natural")
	#plot!(cintc[1],cintc[2],color=:blue2,style=:dashdot, label="cubic clamped, slope at endpoints = 0, 0")
	plot!(cints[1],cints[2],color=:red, label="cubic constrained")
	#plot!(cintcon[1],cintcon[2],color=:black,style=:dash,label="cubic constrained")
	scatter!(x,y,label="knots",color=:black, markersize=5)
end

# ╔═╡ 50b92ab0-5920-11eb-3f85-01c690d30a78
lspline(x,y)

# ╔═╡ c17b386b-d17f-4930-bd4f-bf5b4c0b302a
function constrained(x,y)
	
end

# ╔═╡ 013f0c19-c116-48b3-99af-e9ca7a7974bc
x

# ╔═╡ 6f8fadd2-9157-452f-b66e-7d798077356f
y

# ╔═╡ 7893260e-b68e-470a-b42d-6f9563bd3407
constrainedspline(x,y)

# ╔═╡ 1433e532-2f95-11eb-238b-8d3969489ce0
function conspline(x,y)
    length(x)!=length(y) ? error("Vector length mismatch. Length(x): $(length(x)) != length(y): $(length(y))") : nothing
    length(x) < 3 ? error("Too few points: $(length(x)). Spline needs at least 3 points") : nothing

    n = length(x)#-1
    global Δx= zeros(n)
    Δy= zeros(n)
    σ= zeros(n+1) #slope
    a = zeros(n+1)
    b = zeros(n+1)
    c = zeros(n+1)
    d = zeros(n+1)
    for i in 1:n-1
        Δx[i] = x[i+1] - x[i] #h_i = (x_i+1 - x_i)
        Δy[i] = y[i+1] - y[i]
    end

    @inbounds for i in 1:n-1
        slopechange = y[i]*y[i-1]
        if slopechange > 0.0
			σ[i] = 2.0/(Δx[i]/Δy[i] + Δx[i-1]/Δy[i-1])
		elseif slopechange <= 0.0
			σ[i] = 0.0
		end
    end
	σ[1] = 1.5 * (y[2]-y[1]) / (x[2]-x[1]) - σ[2] / 2.0
    σ[n] = 1.5 * (y[n]-y[n-1]) / (x[n]-x[n-1]) - σ[n-1] / 2.0
	
	κi_1=0.0
	κi=0.0
	
	global κineg=zeros(n)
	global κipos=zeros(n)
	κineg[1]=0.0
	κipos[1]=0.0
    for i in 2:n
        κi_1 = -2*(σ[i]+2*σ[i-1]) / (x[i]-x[i-1]) + 6*(y[i]-y[i-1])/(x[i]-x[i-1])^2 # f''[i-1]
		#κineg[i]=κi_1
        κi =  2*(2*σ[i]+σ[i-1])/ (x[i]-x[i-1]) - 6*(y[i]-y[i-1])/(x[i]-x[i-1])^2 # f''[i]
		#κipos[i]=κi
        #d[i]=(κi-κi_1)/(6*(x[i]-x[i-1]))
		#c[i]=(x[i]*κi_1 - x[i-1]*κi)/(2*(x[i]-x[i-1]))
		b[i]=((y[i]-y[i-1])-c[i]*(x[i]^2-x[i-1]^2)-d[i]*(x[i]^3-x[i-1]^3))/(x[i]-x[i-1])
		a[i] = y[i-1] - b[i]*x[i-1] - c[i]*x[i-1]^2 - d[i]*x[i-1]^3
    end
	CubicSpline(d[2:end],c[2:end],b[2:end],vcat(a[2:end-1],y[end]),x)
end

# ╔═╡ 73c8b780-38a7-11eb-247b-f3e0fe47db66
function cons(x,y)
    length(x)!=length(y) ? error("Vector length mismatch. Length(x): $(length(x)) != length(y): $(length(y))") : nothing
    length(x) < 3 ? error("Too few points: $(length(x)). Spline needs at least 3 points") : nothing

    n = length(x)#-1
    Δx= zeros(n)
    Δy= zeros(n)
    σ= zeros(n) #slope
    a = zeros(n)
    b = zeros(n)
    c = zeros(n)
    d = zeros(n)

    @inbounds for i in 1:n-1
        Δx[i] = x[i+1] - x[i] #h_i = (x_i+1 - x_i)
        Δy[i] = y[i+1] - y[i]
    end

    σ[1] = 3.0*Δy[1] / (2.0*Δx[1]) - σ[2]/2.0
    σ[n] = 3.0*Δy[n-1] / (2.0*Δx[n-1]) - σ[n-1]/2.0
    @inbounds for i in 2:n-1
        slopechange = Δy[i-1]*Δy[i]
        if slopechange > 0.0
			σ[i] = 2.0/(Δx[i]/Δy[i] + Δx[i-1]/Δy[i-1])
		elseif slopechange <= 0.0
			σ[i] = 0.0
		end
    end
    σ[1] = 3.0*(y[2]-y[1])/(2.0*(x[2]-x[1])) - σ[2]/2.0
    σ[n] = 3.0*(y[n]-y[n-1])/(2.0*(x[n]-x[n-1])) - σ[n-1]/2.0

    @inbounds for i in 1:n
        κa =  -2.0*(σ[i+1]+2*σ[i])/Δx[i] + 6*Δy[i]/Δx[i]^2 # f''[i-1]
        κb =  2.0*(2.0*σ[i+1]+σ[i])/Δx[i] - 6*Δy[i]/Δx[i]^2 # f''[i]
        #a[i]=(κb-κa)/(6*Δx[i])
        #b[i]=(x[i+1]*κa - x[i]*κb)/(2*Δx[i])
        c[i]=(Δy[i]- b[i]*(x[i+1]^2 - x[i]^2) - a[i]*(x[i+1]^3 - x[i]^3))/Δx[i]
        d[i] = y[i] - c[i]*x[i] - b[i]*x[i]^2 - a[i]*x[i]^3
    end
    CubicSpline(a,b,c,d,x)
end

# ╔═╡ 5793d204-3531-43c3-aab6-67f9d68b0a1f
ocons(x,y)

# ╔═╡ bdf29d20-49c9-11eb-3993-393f8ddf9ef6


# ╔═╡ c509f3b0-38b2-11eb-0388-a30708ab7cc8
md"Splines"

# ╔═╡ ca6d58ae-38b2-11eb-3fc1-bf442497c9df
md"Interpolations"

# ╔═╡ 8322c440-38b2-11eb-0d19-27627b9a99c0
md"Graphs"

# ╔═╡ 4ce80fc2-2ff2-11eb-3ef4-8315ae577413
md"""

Spline degrees:

* linearspline, lspline, spline1

* quadraticspline, qspline, spline2

* cubicspline, cspline, spline3

* cubicsplineip

Spline styles:

* constrainedspline, naturalspline, clampedspline, periodicspline

Interpolations and extrapolations:

* interp, interptest, extrap, extrapolate

Derivatives:

* slope, deriv1

* curvature, deriv2

* curvrate, deriv3

Types:

* CubicSpline, QuadraticSpline, LinearSpline
"""

# ╔═╡ 5ceaf4c0-87e3-11eb-0821-733526b5a838
x

# ╔═╡ c2cf3461-052d-4265-a7ca-a354c921a160
y

# ╔═╡ c6a21b54-182f-4fc7-af2d-948c9b1d26f6
begin
	slope=zeros(length(x))
	n=length(x)
	for i in 2:n-1
	if (y[i]-y[i-1])*(y[i+1]-y[i]) ≤ 0.0
		slope[i]=0.0
	else
		slope[i]=2.0/((x[i+1]-x[i])/(y[i+1]-y[i]) + (x[i]-x[i-1])/(y[i]-y[i-1]))
	end
	slope[1] = 1.5*(y[2]-y[1])/(x[2]-x[1]) - slope[2]/2.0
    slope[n] = 1.5*(y[n]-y[n-1])/(x[n]-x[n-1]) - slope[n-1]/2.0
	end
	slope
end

# ╔═╡ cd86e220-87ed-11eb-273f-f184755a385c
slope(nns,collect(30.0:1.0:320.0))

# ╔═╡ ca80905a-5845-4ac3-bcef-068e47f366b5
begin
	curveminus=zeros(n)
	curveplus=zeros(n)
	for i in 2:n-1
		cm = 2.0*(slope[i]+2*slope[i-1])/(x[i]-x[i-1]) + 6*(y[i])-y[i-1]/(x[i]-x[i-1])^2
		cp = 2.0*(2.0*slope[i]+slope[i-1])/(x[i]-x[i-1]) - 6*(y[i])-y[i-1]/(x[i]-x[i-1])^2
		curveminus[i]=cm
		curveplus[i]=cp
	end
	curveminus,curveplus
end

# ╔═╡ 04111612-7c8f-47c8-a9a3-c20dcd693ae1
begin
	d=zeros(n)
	c=zeros(n)
	b=zeros(n)
	a=zeros(n)
	for i in 2:n
	cm = -2*(slope[i]+2*slope[i-1])/(x[i]-x[i-1]) + 6*(y[i]-y[i-1])/(x[i]-x[i-1])^2
	cp =2*(2*slope[i]+slope[i-1])/(x[i]-x[i-1]) - 6*(y[i]-y[i-1])/(x[i]-x[i-1])^2
		d[i-1]=(cp-cm)/(6*(x[i]-x[i-1]))
		c[i-1]=(x[i]*cm-x[i-1]*cp)/(2*(x[i]-x[i-1]))
		b[i-1]=((y[i]-y[i-1])-c[i]*(x[i]^2-x[i-1]^2)-d[i]*(x[i]^3-x[i-1]^3))/(x[i]-x[i-1])
		a[i-1]=y[i-1]-b[i]*x[i-1]-c[i]*x[i-1]^2-d[i]*x[i-1]^3
	end
	test=CubicSpline(d[1:n-1],c[1:n-1],b[1:n-1],a[1:n-1],x)
end

# ╔═╡ a455dd2f-bb06-44f4-9988-d0ee0b4df1f3
ctest=interp(test,0:1:100)

# ╔═╡ b25f5b8f-e460-4ab8-80e2-9f49791d6b90
begin
	plot(ctest[1],ctest[2])
	plot!(x,y)
end

# ╔═╡ Cell order:
# ╠═9af1fa20-07fa-11eb-2d8d-239d4cc3f3ae
# ╠═90f87210-59c6-11eb-2d1c-8790aa130324
# ╠═2de193c2-080c-11eb-0451-5d5c46e11324
# ╠═7b664047-b676-4a20-9c8c-5e610497f9fb
# ╟─0c4701a0-dd1c-4f4b-91eb-cb94bc9d4440
# ╟─679c9e8a-ef4f-4b2b-a1ce-245b0ffa773d
# ╠═3eb09f10-63f2-11eb-38c5-db7e4547b561
# ╟─8890e470-63f0-11eb-3df9-19f6bddcf0fd
# ╠═26bca8a0-080c-11eb-33b2-7bf3e34cf10a
# ╟─9240d97e-63eb-11eb-0a53-278a141ee9a6
# ╠═39f7e0b0-63ec-11eb-0a32-eb2d79a226d5
# ╠═cd86e220-87ed-11eb-273f-f184755a385c
# ╟─2f3b8ef0-5920-11eb-3e40-39ee8b03512a
# ╠═6b7f17a0-59c6-11eb-2bd7-81340fa06b7c
# ╟─028b3dc0-63f1-11eb-2b34-ab517656c8d3
# ╠═017de170-5d17-11eb-2964-3b6996cb7eb8
# ╠═a520a3e0-63f1-11eb-3d5c-31ad2e72d89c
# ╠═6b1bc1ee-59c6-11eb-3a5c-2d0c5ab35d0c
# ╠═51f3c956-9462-490d-8375-a5001509a9d5
# ╠═55e3217a-e35e-496e-9f8d-2e042392bb68
# ╠═50b92ab0-5920-11eb-3f85-01c690d30a78
# ╠═c17b386b-d17f-4930-bd4f-bf5b4c0b302a
# ╠═013f0c19-c116-48b3-99af-e9ca7a7974bc
# ╠═6f8fadd2-9157-452f-b66e-7d798077356f
# ╠═7893260e-b68e-470a-b42d-6f9563bd3407
# ╟─1433e532-2f95-11eb-238b-8d3969489ce0
# ╟─73c8b780-38a7-11eb-247b-f3e0fe47db66
# ╠═5793d204-3531-43c3-aab6-67f9d68b0a1f
# ╠═bdf29d20-49c9-11eb-3993-393f8ddf9ef6
# ╟─c509f3b0-38b2-11eb-0388-a30708ab7cc8
# ╟─ca6d58ae-38b2-11eb-3fc1-bf442497c9df
# ╟─8322c440-38b2-11eb-0d19-27627b9a99c0
# ╟─4ce80fc2-2ff2-11eb-3ef4-8315ae577413
# ╠═5ceaf4c0-87e3-11eb-0821-733526b5a838
# ╠═c2cf3461-052d-4265-a7ca-a354c921a160
# ╠═c6a21b54-182f-4fc7-af2d-948c9b1d26f6
# ╠═ca80905a-5845-4ac3-bcef-068e47f366b5
# ╠═04111612-7c8f-47c8-a9a3-c20dcd693ae1
# ╠═a455dd2f-bb06-44f4-9988-d0ee0b4df1f3
# ╠═b25f5b8f-e460-4ab8-80e2-9f49791d6b90
