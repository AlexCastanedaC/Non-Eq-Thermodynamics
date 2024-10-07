#=
Prigogine and Lefever wrote an article in 1967 on
symmetry breaking instabilities in dissipative systems. In this paper,
they developed a simple chemical model that showed how a system can become 
unstable and make a transition to an oscillatory state.

The following piece of code uses the Julia Differential Equations package to
test this model known as the Brusselator.
=#

# We must set our problem up, that is, putting the system of ODEs into the Julia lang
function brusselator!(du, u, p, t)
    k_1, k_2, k_3, k_4, A_c, B_c = p
    X, Y = u
    # dX/dt
    du[1] = k_1 * A_c - k_2 * B_c * X + k_3 * X^2 * Y - k_4 * X
    # dY/dt
    du[2] = k_2 * B_c * X - k_3 * X^2 * Y
end

using DifferentialEquations

println("This code allows you to simulate the Brusselator model")
# Defining the parameters with user input
println("Enter the values of the parameters k_1, k_2, k_3, k_4, A_c, B_c separated by spaces: ")
k_1, k_2, k_3, k_4, A_c, B_c = [parse(Float64, x) for x in split(readline())]
# Defining the initial conditions with user input
println("Enter the initial conditions X_0, Y_0 separated by spaces: ")
x_0, y_0 = [parse(Float64, x) for x in split(readline())]

tspan = (0, 20.0)

p = [k_1, k_2, k_3, k_4, A_c, B_c]
u0 = [x_0, y_0]

prob = ODEProblem(brusselator!, u0, tspan, p)
sol = solve(prob)

# Stationary solutions
X_s = (k_1 / k_4) * A_c
Y_s = (k_4 * k_2 * B_c) / (k_3 * k_1 * A_c)

# Saving the plot of the solution
using Plots

plot(sol; idxs = 1, lc=:red, label = "[X](t)", xlabel = "Time", ylabel = "Concentration", title = "Brusselator Model")
plot!(sol; idxs = 2, lc=:blue, label = "[Y](t)")

# Plotting the stationary solution as dashed lines
plot!(t -> X_s; ls=:dash, lc=:red, label = "Stationary [X]" )
plot!(t -> Y_s; ls=:dash, lc=:blue, label = "Stationary [Y]" )

# Saving the plot
savefig("brusselator.png")




