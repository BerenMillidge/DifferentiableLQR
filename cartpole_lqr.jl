using ControlSystems

A,B,Q,R = randn(4,4),randn(4,4),randn(4,4), randn(4,4)

S = ControlSystems.care(A,B,Q,R)
Sdis = ControlSystems.dare(A,B,Q,R)

function LQR(A,B,Q,R)
    S = ControlSystems.dare(A,B,Q,R)
    println("S: $S")
    K = inv(((transpose(B)*S*B)+R)) * (transpose(B)*S*A)
    return K
end

function simulate_modelled_system(A,B,Q,R,x0,T,Ahat, Bhat)
    K = LQR(Ahat,Bhat,Q,R)
    println("K $K")
    x = [x0]
    u = [0.0]
    for t in 1:T
        u_t = -K * x[t]
        xnext = A*x[t] + (B * u_t)
        push!(x, xnext)
        push!(u, u_t)
    end
    return x,u
end

A = [1 1; 0 1]
B = [0;1]
Q = [1.0 0.0; 0.0 0.0]
R = 1.0
x0 = [-1.0; 0.0]
T = 30
Ahat = A +(0.4 * randn(2,2))
Bhat = B + (0.4 * randn(2))
A*x0 + B*1

x,u = simulate_modelled_system(A,B,Q,R,x0,T,Ahat, Bhat)

using Plots
using LinearAlgebra
plot(u)
plot(x)

plot(cat(x...,dims=2)[2,:])

function simulate_cart(y,m,M,L,g,d,u)
    dy = zeros(4)
    Sy = sin(y[3])
    Cy = cos(y[3])
    D = m*L*L*(M+m*(1-Cy^2));

    dy[1] = y[2];
    dy[2] = (1/D)*(-m^2*L^2*g*Cy*Sy + m*L^2*(m*L*y[4]^2*Sy - d*y[2])) + m*L*L*(1/D)*u;
    dy[3] = y[4];
    dy[4] = (1/D)*((m+M)*m*g*L*Sy - m*L*Cy*(m*L*y[4]^2*Sy - d*y[2])) - m*L*Cy*(1/D)*u +.01*randn();
    return dy
end

m =1
M = 5
L = 2
g = -10
d = 1
s = 1
A = [0.0000001 1 0 0;
    0 -d/M -m*g/M 0;
    0 0 0 1;
    0 -s*d/(M*L) -s*(m+M)*g/(M*L) 0]


B = [0; 1/M; 0; s*1/(M*L)]
#B = randn(4)

Q = [1 0 0 0;
    0 1 0 0;
    0 0 10 0;
    0 0 0 100]
R = 0.1

dt = 100
x0 = [0; 0; pi+0.1; 0]

function simulate_cartpole(A, B,Q,R,x0,T,dt,m,M,L,g,d)
    println("R: $R")
    S = ControlSystems.dare(A,B,Q,R)
    println("S: $S")
    K = inv(((transpose(B)*S*B)+R)) * (transpose(B)*S*A)
    println("K:")
    println(K)
    x = [x0]
    u = [0.0]
    for t in 1:T
        t = Int(t)
        x_t = x[t]
        u_t = -K * (x_t .- [1;0;pi;0])
        push!(u, u_t)
        for tt in 1:dt
            dx = simulate_cart(x_t,m,M,L,g,d,u_t)
            u_t = 0# allow to evolve smoothly from there?
            println(dx)
            x_t  += (dx/dt)
        end
        push!(x, x_t)
    end
    return x,u
end

simulate_cartpole(A,B,Q,R,x0,T,dt,m,M,L,g,d)
