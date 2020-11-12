using SSpline
export SSpline, clampedspline, constrainedspline, cspline, cubicspline, cubicsplineip, curvature, curvrate, deriv1, deriv2, deriv3, extrap, interp, interptest, linearspline, lspline, naturalspline, periodicspline, qspline, quadraticspline, slope, spline1, spline2, spline3

function testconstrained()
    x=[0,10,30,50,70,90,100];
    y=[30,130,150,150,170,220,320];
    atest=[30, 109.9]
    btest=[14.09, 2.3181]
    ctest=[0, -0.01818]
    dtest=[-0.0409, -0.000454]

    csp=constrainedspline(x,y)
    eps=1e-6
    for i in 1:2
        println("csp.d[$i]: $(csp.d[i]) : $(dtest[i])")
        csp.d[i]-dtest[i] < eps
        println("csp.c[$i]: $(csp.c[i]) : $(ctest[i])")
        csp.c[i]-ctest[i] < eps
        println("csp.b[$i]: $(csp.b[i]) : $(btest[i])")
        csp.b[i]-btest[i] < eps
        println("csp.a[$i]: $(csp.a[i]) : $(atest[i])")
        csp.a[i]-atest[i] < eps
    end
end