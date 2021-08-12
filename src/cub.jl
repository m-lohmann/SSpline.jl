function cub(x::Vector, y::Vector, style...)
    h[i] = x[i] - x[i-1]

    s[i] = 1/h[i] * ( (x[i] - xs)^3 / 6 * M[i-1] + h[i]^3 / 3 * M[i] + y[i-1] - h[i]^2 / 6 * M[i-1] * x[i] - xs + y[i] - h[i]^2 / 6  * M[i] * h[i] )

    return CubicSpline(d, c, b, y, x)
end