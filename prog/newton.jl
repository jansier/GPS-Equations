"""
Compute (x-a)^2 + (y-b)^2 + (z-c)^2 - C^2(t-T)^2
Given [a,b,c,t] and [x,y,z,T]
"""
function f(cv, v)
    return (v[1]-cv[1])^2 + (v[2]-cv[2])^2 + (v[3]-cv[3])^2 - C*C*((v[4]-cv[4])^2)
end

"""
Compute vector of f values
Given matrix of constants and vector of variables
"""
function F(mc, v)
    V = zeros(v)
    V[1] = f(mc[1,1:4], v)
    V[2] = f(mc[2,1:4], v)
    V[3] = f(mc[3,1:4], v)
    V[4] = f(mc[4,1:4], v)
    return V
end

"""
Compute derivative of (x-c)^2 in x
"""
function df(c, x)
    return 2*x-2*c
end

"""
Compute gradient of
(x-a)^2 + (y-b)^2 + (z-c)^2 - C^2(t-T)^2 in [x, y, z]
Given [a, b, c] and [x, y, z]
"""
function gradient(cv, v)
    g = zeros(v)
    g[1] = df(cv[1], v[1])
    g[2] = df(cv[2], v[2])
    g[3] = df(cv[3], v[3])
    g[4] = -C^2*df(cv[4], v[4])
    return g
end

"""
Compute jacobian of GPS equations in [x, y, z, T]
Given matrix 4x4 with rows like
[a1, b1, c1, t1] and vector v
"""
function jacobian(cm, v)
    j = zeros(cm)
    j[1, 1:4] = gradient(cm[1, 1:4], v)
    j[2, 1:4] = gradient(cm[2, 1:4], v)
    j[3, 1:4] = gradient(cm[3, 1:4], v)
    j[4, 1:4] = gradient(cm[4, 1:4], v)
    return j
end

"""
Find solution for GPS equations with Newton method
Given matrix 4x4 and starting vector
(default: speed of light, numbers type)
"""
function newton(mc, v=[0,0,0,0], c=299792458, T=Float64)
    global C = c
    v = convert(Array{T,1}, v)
    mc = convert(Array{T,2}, mc)
    # print("v0: ", v, "\n")

    for i in 1:20

        Fv = F(mc, v)
        J = jacobian(mc, v)
        h = -inv(J)*Fv
        v = v+h
        # print("v", i, ": ", v, "\n")
    end
    return v
end
