### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ f041ef20-f968-11eb-1013-ef1b063dc4cc
begin
	using Pkg
	Pkg.activate(".")
	using Plots
	using SSpline
end

# ╔═╡ 0d9037b8-b0ff-4a4a-8d7f-80b29338ac19
html"""<style> main { max-width: 900px;} """

# ╔═╡ 470b7ee9-b8f7-4c8d-abfd-7fa871abc36f
x = collect(1:10);y = rand(10);

# ╔═╡ d0253bf6-4564-4c3f-9c80-3f59bf3f5370
md"""

(1)

    s[i](x)    = a[i] * (x-x[i])^3 + b[i] * (x-x[i])^2 + c[i] * (x-x[i]) + d[i]

----

    s'[i](x)   = 3 * a[i] * (x-x[i])^2 + 2 * b[i] * (x-x[i]) + c[i]

    s''[i](x)  = 6 * a[i] * (x-x[i]) + 2 * b[i]

    s'''[i](x) = 6 * a[i]

(2)

    s_i(x_i)     = y_i, i=0..n-1

																s_n-1 = y_n

(3)

    s_i(x_i+1)   = s_i+1(x_i+1),      i = 0..n-2

(4)

    s'_i(x_i+1)  = s'_i+1(x_i+1),     i = 0..n-2

(5)

    s''_i(x_i+1) = s''_i+1(x_i+1),    i = 0..n-2


Natural

    s''_0(x_0) = 0
    s''_n-1(x_n) = 0

Clamped

    s'_0(x_0) = y'_0
    s'_n-1(x_n) = y'_n

Periodic

    s'_0(x_0) = s'_n-1(x_n)
    s''_0(x_0) = s''_n-1(x_n)

NAK

    s'''_0(x_1) = s'''_1(x_1)
    s'''_n-2(x_n-1) = s'''_n-1(x_n-1)

----

	h = x_i+1 - x_i

from (1) `s[i](x) = a[i] * (x-x[i])^3 + b[i] * (x-x[i])^2 + c[i] * (x-x[i]) + d[i]`

    d_i = y_i

    σ_i = s''(x_i) = 6 * a_i * (x_i - x_i) + 2 * b_i = 2 * b_i,                          i = 0..n

    σ_i = 2 * b_i

    b_i = σ_i / 2

----

    σ_i+1 = 6 * a_i * h + 2 * b_i

    σ_i+1 = 6 * a_i * h + 2 * σ_i / 2

from (5) `s''_i(x_i+1) = s''_i+1(x_i+1),    i = 0..n-2`

	a_i = (σ_i+1 - σ_i) / (6 * h)

from (3) `s_i(x_i+1) = s_i+1(x_i+1),      i = 0..n-2`

    y_i+1 = a_i * h^3 + b_i * h^2 + c_i * h + y_i

    y_i+1 = (σ_i+1 - σ_i) / (6 * h) * h^3 + σ_i / 2 * h^2 + c_i * h + y_i

    y_i+1 - y_1 = (σ_i+1 - σ_1) / (6 * h) * h^3 + σ_i / 2 * h^2) + c_i * h

    (y_i+1 - y_1) / h = (σ_i+1 - σ_1) / 6 * h + σ_i / 2 * h + c_i

    (y_i+1 - y_1) / h - h * ( (σ_i+1 - σ_1) / 6 + σ_i / 2) = c_i

	c_i = (y_i+1 - y_1) / h - h * ( (σ_i+1 - σ_1) / 6 + 3σ_i / 6 )

    c_i = (y_i+1 - y_1) / h - h * (σ_i+1 - σ_1 + 3 * σ_i) / 6

    c_i = (y_i+1 - y_1) / h - h * (2 * σ_i + σ_i+1) / 6

from (4) `s'_i(x_i+1)  = s'_i+1(x_i+1),     i = 0..n-2`

	c_i+1 = 3 * a_i * h^2 + 2 * b_i * h + c_i

Insert a_i, b_i, c_i, and use

from (2) `s_i(x_i) = y_i, i=0..n-1)`

    c_i+1 = 3 * (σ_i+1 - σ_i) / (6 * h) * h^2 + 2 * σ_i / 2 * h + (y_i+1 - y_1) / h - h * (2 * σ_i + σ_i+1) / 6

	
"""

# ╔═╡ c82ec4aa-827e-40f0-ab10-543380e5d5d7
begin
	range = collect(1:0.1:10)
	
	nak = nakspline(x,y)
	kx, ky = interp(nak,range)
	
	nat = spline3(x,y,:natural)
	nx, ny = interp(nat,range)
end

# ╔═╡ 88c79600-cc20-4ad9-9a98-6d05e2f7512e
nakspline(x,y)

# ╔═╡ 10458fe7-350f-481f-99ca-29e5c185d3b9


# ╔═╡ Cell order:
# ╠═f041ef20-f968-11eb-1013-ef1b063dc4cc
# ╠═0d9037b8-b0ff-4a4a-8d7f-80b29338ac19
# ╠═470b7ee9-b8f7-4c8d-abfd-7fa871abc36f
# ╟─d0253bf6-4564-4c3f-9c80-3f59bf3f5370
# ╟─c82ec4aa-827e-40f0-ab10-543380e5d5d7
# ╠═88c79600-cc20-4ad9-9a98-6d05e2f7512e
# ╠═10458fe7-350f-481f-99ca-29e5c185d3b9
