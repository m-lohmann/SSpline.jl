### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 1ff97a4d-92bf-4f4f-ab60-c5f311b6b2b4
using Plots, SSpline

# ╔═╡ 3a8ef900-d5d9-11eb-1e22-edc33d208d0e
#import Pkg; Pkg.add(path="c:/Users/4vektor/.julia/dev/SSpline")

# ╔═╡ 1e7150a7-40ce-44e7-a719-fa995c10ac86
let
    #import Pkg
	#Pkg.activate(".")
    #Pkg.add("SSpline")
end

# ╔═╡ fcc81fdb-a8c0-489c-b36e-4065a91dca6e
pyplot()

# ╔═╡ c56bd9d8-6829-47c8-b610-e0bb0c1f88ed
html"""<style> main { max-width: 900px;} """

# ╔═╡ 101c094a-b0eb-46fb-a06d-1415de4018da
begin
	x=[10,20,35,50,55,57,60];      y=[5,30,20,22,18,30,12];
	#x=[0,10,30,50,70,90,100];     y=[30,130,150,150,170,220,320];
	#x=[0.0,25.0,50.0,75.0,100.0]; y=[20.0,50.0,10.0,60.0,35.0];
	#x=[0,10,30,45,70,90,100];     y=[20,50,70,-70,60,50,20];
	#x=[0,10,30,50,70,80,82];      y=[150,200,200,200,180,100,0];
end

# ╔═╡ ee8e5634-b8c3-4a99-8dad-2ee228221609
begin
	ra = x[1]:0.2:x[end] # interpolation range
end

# ╔═╡ a8eb6119-8ba8-4a60-9be4-7787731da630
begin
	nn = spline0(x,y)
	ls = spline1(x,y)
	ni = interp(nn, ra)
	li = interp(ls, ra)
end

# ╔═╡ 2cf7565e-6d87-4219-bd8b-64a0a3e7ec41
begin
	plot(size = (850,500), legend = :bottom, title = "nearest neighbor and linear interpolations")
	scatter!(x, y, label = "knots")
	plot!(ni[1], ni[2], color = :green, label = "nearest neighbor")
	plot!(li[1], li[2], color = :blue, style = :dashdot, label = "linear")
end

# ╔═╡ 8081af65-98bc-4160-8b73-915f671e21a9
begin
	qs = spline2(x,y)
	qs2 = spline2(x,y, 6.1)
	qi = interp(qs, ra)
	qi2 = interp(qs2, ra)
	qsl = slope(qs2, ra)
end

# ╔═╡ 4f01a157-9a1a-4bfc-b367-df142313e7e6
begin
	plot(size = (850,500), legend = :topleft, title = "linear and quadratic interpolations")
	scatter!(x, y, label = "knots")
	plot!(li[1], li[2], color = :grey, label = "linear")
	plot!(qi[1], qi[2], color = :orange, style = :dashdot, label = "quadratic, f'(x₁) = 0.0 (default)")
	plot!(qi2[1], qi2[2], color = :blue, thickness = 2.0, style = :dash, label = "quadratic, f'(x₁) = 6.1")
end

# ╔═╡ 0bcea8ae-002b-4fb6-be6a-89a4b1abadb1
begin
	csn = spline3(x,y, :natural)
	csc = spline3(x,y, :clamped, -5.0, 5.0)
	csp = spline3(x,y,:notaknot)
	ci = interp(csn,ra)
	cic = interp(csc, ra)
	csl = slope(csn, ra)
	cci = interp(csp,ra)
end

# ╔═╡ f72de351-b2ae-4a96-b6ef-d8627a45c70a
begin
	plot(size = (850,500), legend = :topleft, title = "quadratic vs. cubic natural")
	scatter!(x, y, label = "knots", color = :white, markersize = 10, markershape = :hexagon)
	plot!(qi2[1], qi2[2], color = :green, style = :dashdot, thickness = 2.0, label = "quadratic, slope 6.1 at start")
	plot!(ci[1], ci[2], color = :purple, style = :dash, thickness = 2.0, label = "cubic natural")
	plot!(cic[1], cic[2], color = :orange2, thickness = 2.0, label = "cubic clamped (-5.0, 5.0)")
	plot!(cci[1], cci[2], color = :skyblue, label = "notaknot f")
end

# ╔═╡ 03e94833-ab83-46c8-9eb9-4612727d072d
begin
	xh=[0,10,30,50,70,90,100]
	yh=[30,130,150,150,170,220,320]
	
	rah = xh[1]:0.2:xh[end] 
	
	hk = hspline(xh, yh, :kruger)
	hki = interp(hk, rah)
	
	hf = hspline(xh, yh, :free)
	hfi = interp(hf, rah)
	
	#hg = hspline(xh, yh, :geom)
	#hgi = interp(hg, rah)
end

# ╔═╡ 16be7b7b-c784-4519-8ec8-13a8c1ac2d26
begin
	plot(size = (850,500), legend = :bottom, title = "hermite splines")
	scatter!(xh, yh, label = "knots", color = :white, markersize = 10, markershape = :hexagon)
	plot!(hki[1], hki[2], color = :orange2, label = "hermite kruger (harmonic mean)")
	#plot!(hgi[1], hgi[2], color = :blue, style = :dashdot, label = "hermite geometric mean")
	plot!(hfi[1], hfi[2], color = :green, style = :dash, label = "hermite free")
end

# ╔═╡ 5cefbebb-fe49-4ef1-b3fb-5e6577341383
begin
	cel = extrapolate(hk, -10:110, :linear)
	ces = extrapolate(hk, -10:110, :smooth)
	ceil = interp(cel, -10:0.2:110)
	ceis = interp(ces, -10:0.2:110)
end

# ╔═╡ df66e2da-d6a0-4afe-8fc6-dc36398a5bc2
begin
	plot(size = (850,400), legend = :bottom, title = "hermite extrapolations")
	scatter!(xh, yh, label = "knots", color = :white, markersize = 10, markershape = :hexagon)
	plot!(hki[1], hki[2], color = :orange, width = 2.0, label = "hermite kruger")
	plot!(ceil[1], ceil[2], color = :skyblue, style = :dash, label = "herm. kruger linear extr.")
	plot!(ceis[1], ceis[2], color = :green, style = :dashdot, label = "herm. kruger smooth extr.")
end

# ╔═╡ 7fc4737e-9b75-4d53-bae1-f977a1466d8e
begin
	cnn = spline2(xh,yh)
	ceq = extrapolate(cnn, -10:110, :zero)
	ceiq = interp(ceq, -10:0.2:110)
	plot(size = (850,400), legend = :bottom, title = "cubic extrapolations")
	scatter!(xh, yh, label = "knots", color = :white, markersize = 10, markershape = :hexagon)
	plot!(ceiq[1], ceiq[2], color = :red, style = :dot, label = "cub. nat. lin.")
	scatter!(ceq.x, ceq.a)
end

# ╔═╡ a5fb28d3-4b1b-4bda-8d20-1b5796f1374f
ceiq[2][end]

# ╔═╡ 680e2580-f4bf-47d1-bca8-19086f017b04
ceq

# ╔═╡ d7189e79-63d5-46bd-bf8c-d37f5ed8d50c
ceq.a[end-1]

# ╔═╡ 014045fc-e2b2-48cf-aeb6-0e71cbe18426
begin
	#xx=[0,10,30,50,70,90,100];	yy=[30,130,150,150,170,220,320]
	#xx=[-30,-15, -10,-7, -5, 0, 5, 7, 10, 15, 30]; yy=[0, 20, 15, 20, 8, 8, 8, 20, 15, 20, 0];
	#xx=[10,18,20,35,50,55,57,60]; yy=[0,20,20,20,22,18,30,12];
	xx=[10,18,20,35,50,55,57,60]; yy=[12,30,18,22,20,20,20,0];
	x2=[10,18,20,35,50,55,57,60]; y2=[0,20,20,20,22,18,30,12];
	cn = cspline(xx,yy,:natural)
	cp = cspline(xx,yy,:parabolic)
	cc = cspline(xx,yy,:notaknot)
	
	nx,ny = interp(cn, xx[1]:1:xx[end])
	n1x,n1y = deriv1(cn, xx[1]:1:xx[end])
	n2x,n2y = deriv2(cn, xx[1]:1:xx[end])
	n3x,n3y = deriv3(cn, xx[1]:1:xx[end])
	nkx,nky = deriv2(cn, xx)
	
	px,py = interp(cp, xx[1]:1:xx[end])
	p1x,p1y = deriv1(cp, xx[1]:1:xx[end])
	p2x,p2y = deriv2(cp, xx[1]:1:xx[end])
	p3x,p3y = deriv3(cp, xx[1]:1:xx[end])
	pkx,pky = deriv2(cp, xx)
	
	cx,cy = interp(cc, xx[1]:1:xx[end])
	c1x,c1y = deriv1(cc, xx[1]:1:xx[end])
	c2x,c2y = deriv2(cc, xx[1]:1:xx[end])
	c3x,c3y = deriv3(cc, xx[1]:1:xx[end])
	ckx,cky = deriv2(cc, xx)
	
	plot(size = (850,500), xlimits = (xx[1], xx[end]), legend = :outertopright, title = "extrapolations", fg_color_grid = :black, gridalpha =  0.5)
	scatter!(xx, yy, label = "knots", color = :white, markersize = 10, markershape = :hexagon)
	#plot!(e,f, color = :blue)
	#plot!(hx,hy, color = :red, label = "natural")
	plot!(nx,ny, color = :blue, style = :solid, label = "natural f")
	#plot!(n1x,10*n1y, color = :blue, style = :dash, label = "natural f'")
	#plot!(n2x,10*n2y, color = :blue, style = :dashdot, label = "natural f''")
	#plot!(n3x,10*n3y, color = :blue, style = :dot, label = "natural f'''")
	#scatter!(nkx, nky, color = :blue)
	plot!(px,py, color = :red, style = :dash, label = "parabolic f")
	#plot!(p1x,10*p1y, color = :red, style = :solid, label = "parabolic f'")
	#plot!(p2x,10*p2y, color = :red, style = :dashdot, label = "parabolic f''")
	#plot!(p3x,400*p3y, color = :red, style = :dot, label = "parabolic f'''")
	#scatter!(pkx, pky, color = :red)
	plot!(cx, cy,color = :orange, style = :dash, label = "notaknot f")
	#plot!(c1x,10*c1y, color = :green, style = :solid, label = "notaknot f'")
	#plot!(c2x,10*c2y, color = :orange, style = :dashdot, label = "notaknot f''")
	#plot!(c3x,10*c3y, color = :orange, style = :dot, label = "notaknot f'''")
	#scatter!(nkx, nky, color = :green)
end

# ╔═╡ 2ef87041-69f6-4e37-8a46-22e480a2b417
begin
        n = length(x)-1
        Δx= Vector{Float64}(undef, n-1)
        Δy= Vector{Float64}(undef, n-1)
        λ = Vector{Float64}(undef, n)    # upper diagonal    
        δ = ones(n) * 2 # main diagonal
        μ = Vector{Float64}(undef, n)    # lower diagonal
        r = Vector{Float64}(undef, n)    # derivatives
        M = Vector{Float64}(undef, n)

        @inbounds for i in 1:n-1
            Δx[i] = x[i+1] - x[i]
            Δy[i] = y[i+1] - y[i]
        end
            A = zeros(n, n)
            @inbounds for i in 2:n-2
                A[i, i-1] = Δx[i-1]
                A[i, i] = 2.0 * (Δx[i-1] + Δx[i])
                #A[i, i] = 2.0 * (Δx[i-1] + Δx[i])
                A[i, i+1] = Δx[i+1]
            end
            #A[1, 1] = -Δx[2]
            #A[1, 2] = Δx[1] + Δx[2]
            #A[1, 3] = -Δx[1]
            A[begin, begin] = -Δx[begin + 2]
            A[begin, begin + 1] = Δx[begin + 1] + Δx[begin + 2]
            A[begin, begin + 2] = -Δx[begin + 1]

            A[end, end-2] = -Δx[end - 1]
            A[end, end-1] = Δx[end - 2] + Δx[end - 1]
            A[end, end]   = -Δx[end - 2]
            for i in 2:n-1
            #    r[i] = 6.0 / (Δx[i-1] + Δx[i]) * (Δy[i] / Δx[i] - Δy[i-1] / Δx[i-1])
                r[i] = 6.0 * (Δy[i] / Δx[i] - Δy[i-1] / Δx[i-1])
            end
end

# ╔═╡ b4a2a3c5-7f10-4330-b7ad-adbb54addf01
A

# ╔═╡ Cell order:
# ╠═3a8ef900-d5d9-11eb-1e22-edc33d208d0e
# ╠═1e7150a7-40ce-44e7-a719-fa995c10ac86
# ╠═1ff97a4d-92bf-4f4f-ab60-c5f311b6b2b4
# ╠═fcc81fdb-a8c0-489c-b36e-4065a91dca6e
# ╠═c56bd9d8-6829-47c8-b610-e0bb0c1f88ed
# ╠═101c094a-b0eb-46fb-a06d-1415de4018da
# ╠═ee8e5634-b8c3-4a99-8dad-2ee228221609
# ╠═a8eb6119-8ba8-4a60-9be4-7787731da630
# ╠═2cf7565e-6d87-4219-bd8b-64a0a3e7ec41
# ╠═8081af65-98bc-4160-8b73-915f671e21a9
# ╠═4f01a157-9a1a-4bfc-b367-df142313e7e6
# ╠═0bcea8ae-002b-4fb6-be6a-89a4b1abadb1
# ╠═f72de351-b2ae-4a96-b6ef-d8627a45c70a
# ╠═03e94833-ab83-46c8-9eb9-4612727d072d
# ╠═16be7b7b-c784-4519-8ec8-13a8c1ac2d26
# ╠═5cefbebb-fe49-4ef1-b3fb-5e6577341383
# ╠═df66e2da-d6a0-4afe-8fc6-dc36398a5bc2
# ╠═7fc4737e-9b75-4d53-bae1-f977a1466d8e
# ╠═a5fb28d3-4b1b-4bda-8d20-1b5796f1374f
# ╠═680e2580-f4bf-47d1-bca8-19086f017b04
# ╠═d7189e79-63d5-46bd-bf8c-d37f5ed8d50c
# ╠═014045fc-e2b2-48cf-aeb6-0e71cbe18426
# ╠═2ef87041-69f6-4e37-8a46-22e480a2b417
# ╠═b4a2a3c5-7f10-4330-b7ad-adbb54addf01
