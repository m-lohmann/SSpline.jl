using SSpline
using Test

@testset "" begin
    

end

@testset "Nearest Neighbor Spline" begin
    x=[0,1,2,3,4]
    y=[0,1,5,2,2]
    @test typeof(nearestneighbor(x,y))==NearestNeighborSpline
end

@testset "Nearest Neighbor Spline Naming" begin
    x=[0,1,2,3,4]
    y=[0,1,5,2,2]
    nn1=nearestneighbors(x,y)
    nn2=stepspline(x,y)
    nn3=spline0(x,y)
    @test nn1.a == [0.0,1.0,5.0,2.0,2.0]
    @test nn1.x == [0.0,1.0,2.0,3.0,4.0]

    @test nn2.a == nn1.a
    @test nn2.x == nn1.x
    @test nn3.a == nn1.a
    @test nn3.x == nn1.x
end

@testset "Linear Spline Type" begin
    x=[0,1,2,3,4]
    y=[0,1,5,2,2]
    @test typeof(linearspline(x,y))==LinearSpline
end

@testset "Linear Spline Naming" begin
    x=[0,1,2,3,4]
    y=[0,1,5,2,2]

    lin1=linearspline(x,y)
    lin2=lspline(x,y)
    lin3=spline1(x,y)

    # name tests, spline generation

    @test lin1.b == [1.0, 4.0, -3.0, 0.0, 0.0]
    @test lin1.a == [0.0,1.0,5.0,2.0,2.0]
    @test lin1.x == [0.0,1.0,2.0,3.0,4.0]

    @test lin2.b == lin1.b
    @test lin2.a == lin1.a
    @test lin2.x == lin1.x

    @test lin3.b == lin1.b
    @test lin3.a == lin1.a
    @test lin3.x == lin1.x
end

@testset "Linear Spline Interpolation" begin
    # interpolation test
    x=[0,1,2,3,4]
    y=[0,1,5,2,2]

    lin1=linearspline(x,y)

    intlin=interp(lin1,collect(0.0:0.5:4.0))

    @test intlin == ([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0], [0.0, 0.5, 1.0, 3.0, 5.0, 3.5, 2.0, 2.0, 2.0])

end

@testset "Cubic Spline Type" begin
    x=[0,1,2,3,4]
    y=[0,1,5,2,2]

    @test typeof(cubicspline(x,y)) == CubicSpline
end=

