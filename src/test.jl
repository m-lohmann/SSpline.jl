using SSpline, Plots

function splinetest()

    x=collect(1:10)
    y=[1.0,1.65,3.0,1.95,1.5,2.5,0.5,0.5,-3.75,0.0]
    int=collect(1.0:0.1:10.0)
    #Spline creation
    nn=nearestneighbor(x,y)
    ln=linearspline(x,y)
    cn=cubicspline(x,y,:natural)
    cc=cubicspline(x,y,:clamped,0.0,0.0)
    co=constrainedspline(x,y)

    #Spline interpolation
    nnx,nny=interp(nn,int)
    lnx,lny=interp(ln,int)
    cnx,cny=interp(cn,int)
    ccx,ccy=interp(cc,int)
    cox,coy=interp(co,int)

    #Plots
    scatter(x,y,label="knots")
    plot!(nnx,nny,label="nearest neighbor")
    plot!(lnx,lny,label="linear")
    plot!(cnx,cny,color=:orange,label="cubic natural")
    plot!(ccx,ccy,linestyle= :dashdot,color=:blue,label="cubic clamped (0,0)")
    plot!(cox,coy,linestyle= :dot, color=:red,label="cubic constrained")
end

splinetest()
