include("newton.jl")
include("alg.jl")

R = 6371000
OR = 20180000
C = 299792458


"""
Returns coordinates (x,y,z) on earth surface
given (x, y)
"""
function surfaceCords(x, y)
    return (x, y, sqrt(R*R - x*x - y*y))
end

"""
Returns coordinates (x,y,z) of a satelite
given (x, y)
"""
function sateliteCords(x, y)
    return (x, y, sqrt(OR*OR - x*x - y*y))
end

"""
Return travel time of light between points
"""
function time(p, s)
    return distance(p, s)/C
end

"""
Given x, y of receiver and satelites,
returns matrix of GPS equations
(Calculates z coordinates and time)
"""
function generateMatrix(p, x1, y1, x2, y2, x3, y3, x4, y4, t)
    s1 = sateliteCords(x1, y1)
    s2 = sateliteCords(x2, y2)
    s3 = sateliteCords(x3, y3)
    s4 = sateliteCords(x4, y4)
    m = Array{Float64, 2}(4, 4)
    m[1, 1:4] = [s1[1], s1[2], s1[3], t+time(p, s1)]
    m[2, 1:4] = [s2[1], s2[2], s2[3], t+time(p, s2)]
    m[3, 1:4] = [s3[1], s3[2], s3[3], t+time(p, s3)]
    m[4, 1:4] = [s4[1], s4[2], s4[3], t+time(p, s4)]

    return m
end

"""
Generates random (x, y), given radius
"""
function randomXY(r)
    x = rand()*10^8%r
    y = rand()*10^8%floor(sqrt(r*r-x*x))
    return (x, y)
end

"""
Returns distance between points
"""
function distance(p, s)
    return sqrt((p[1]-s[1])^2+(p[2]-s[2])^2+(p[3]-s[3])^2)
end

"""
Returns random equations data
(matrix of equations, receiver coordinates)
"""
function randomData()
    (x, y) = randomXY(R)
    p = surfaceCords(x, y)
    t = rand()
    (x1, y1) = randomXY(OR)
    (x2, y2) = randomXY(OR)
    (x3, y3) = randomXY(OR)
    (x4, y4) = randomXY(OR)
    m = generateMatrix(p, x1, y1, x2, y2, x3, y3, x4, y4, t)
    return (m, [p[1], p[2], p[3], t])
end

"""
Given GPS matrix, solution
returns Newtons error
"""
function testNewton(m, v)
    pr = newton(m)
    if distance((v[1], v[2], v[3]), pr) > 1 # second solution
        pr = newton(m, v)
    end
    return distance(v, pr)
end


"""
Given GPS matrix, solution
returns Algebraics method error
"""
function testAlgebraic(m, v)
    (p1, p2) = algebraic(m)
    return min(distance(v, p1), distance(v, p2))
end

"""
n - number of tests
"""
function bigTest(n = 100)
    maxNewtonError = 0
    maxAlgebraicError = 0

    for i in 1:n
        (m, p) = randomData()
        errorNewton = testNewton(m, p)
        if errorNewton > maxNewtonError
            maxNewtonError = errorNewton
        end
        errorAlgebraic = testAlgebraic(m, p)
        if errorAlgebraic > maxAlgebraicError
            maxAlgebraicError = errorAlgebraic
        end
    end
    return (maxNewtonError, maxAlgebraicError)
end
