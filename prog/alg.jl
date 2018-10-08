function algebraic(m, C=299792458)

    q1 = m[1,1:4]
    q2 = m[2,1:4]
    q3 = m[3,1:4]
    q4 = m[4,1:4]

    a₁ = q1[1]; b₁ = q1[2]; c₁ = q1[3]; t₁ = q1[4];
    a₂ = q2[1]; b₂ = q2[2]; c₂ = q2[3]; t₂ = q2[4];
    a₃ = q3[1]; b₃ = q3[2]; c₃ = q3[3]; t₃ = q3[4];
    a₄ = q4[1]; b₄ = q4[2]; c₄ = q4[3]; t₄ = q4[4];

    U⃗x = [2(a₂-a₁), 2(a₃-a₁), 2(a₄-a₁)]
    U⃗y = [2(b₂-b₁), 2(b₃-b₁), 2(b₄-b₁)]
    U⃗z = [2(c₂-c₁), 2(c₃-c₁), 2(c₄-c₁)]
    U⃗T = [2*(C*C)*(t₁-t₂), 2*(C*C)*(t₁-t₃), 2*(C*C)*(t₁-t₄)]

    temp₁ = a₁*a₁ - a₂*a₂ + b₁*b₁ - b₂*b₂ + c₁*c₁ - c₂*c₂ + C*C*t₂*t₂- C*C*t₁*t₁
    temp₂ = a₁*a₁ - a₃*a₃ + b₁*b₁ - b₃*b₃ + c₁*c₁ - c₃*c₃ + C*C*t₃*t₃- C*C*t₁*t₁
    temp₃ = a₁*a₁ - a₄*a₄ + b₁*b₁ - b₄*b₄ + c₁*c₁ - c₄*c₄ + C*C*t₄*t₄- C*C*t₁*t₁

    W⃗ = [temp₁, temp₂, temp₃]

    l₁ = det([U⃗y U⃗z U⃗x])
    l₂ = det([U⃗y U⃗z U⃗T])
    l₃ = det([U⃗y U⃗z W⃗])

    m₁ = det([U⃗x U⃗z U⃗y])
    m₂ = det([U⃗x U⃗z U⃗T])
    m₃ = det([U⃗x U⃗z W⃗])

    n₁ = det([U⃗x U⃗y U⃗z])
    n₂ = det([U⃗x U⃗y U⃗T])
    n₃ = det([U⃗x U⃗y W⃗])

    e₁ = -l₂/l₁
    e₂ = -l₃/l₁
    e₃ = -m₂/m₁
    e₄ = -m₃/m₁
    e₅ = -n₂/n₁
    e₆ = -n₃/n₁

    G = e₁*e₁ + e₃*e₃ + e₅*e₅ - C*C
    H = 2(e₁*e₂+e₃*e₄+e₅*e₆-a₁*e₁-b₁*e₃-c₁*e₅+C*C*t₁)
    I = e₂*e₂+e₄*e₄+e₆*e₆-2(a₁*e₂+b₁*e₄+c₁*e₆)+a₁*a₁+b₁*b₁+c₁*c₁-C*C*t₁*t₁

    Δ = H*H-4*G*I

    T₁ = (-H+√Δ)/(2G)
    T₂ = (-H-√Δ)/(2G)

    x₁ = e₁*T₁+e₂
    y₁ = e₃*T₁+e₄
    z₁ = e₅*T₁+e₆

    x₂ = e₁*T₂+e₂
    y₂ = e₃*T₂+e₄
    z₂ = e₅*T₂+e₆

    return ([x₁,y₁,z₁,T₁], [x₂,y₂,z₂,T₂])
end
