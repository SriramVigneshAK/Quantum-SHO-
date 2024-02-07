#This is a code to calculate the eigenstates of a Simple Harmonic Oscillator
using QuantumOptics
using Plots

#Mass of particle and length of box
m = 1 #mass
k = 1 #Force Constant
n = 1 #Quantum number

#Defining a Basis to represent everything
xmin = -10
xmax = 10
Npoints = 100
b_position = PositionBasis(xmin, xmax, Npoints) #Position Basis
b_momentum = MomentumBasis(b_position) #Momentum Basis

#Defining Operators 
x = position(b_position) #Position Operator in Position Basis
p = momentum(b_momentum) #Momentum Operator in Momentum Basis
px = momentum(b_position) #Momentum Operator in Position Basis

#Defining the Potential 
function V(x)
    v = 0.5*k*x^2
    return v
end

#Plotting the Potential
PE = potentialoperator(b_position, V)
PE = dense(PE)
ptsx = samplepoints(b_position)
#plot(ptsx , V)


#Defining the Hamiltonian
H = LazySum(px^2/2m , PE)
# A more feasible way is to use Fourier transforms(FT), so defining FT operators
Txp = transform(b_position, b_momentum)
Tpx = transform(b_momentum, b_position)
Hkin = LazyProduct(Txp, p^2/2m, Tpx) #A faster method to do Hkin * psi
H = dense(LazySum(Hkin, PE))

#finding eigenstates
E,ψ = eigenstates(H) #Computes all  eigenstates
En,ψn = eigenstates(H,n) #Computes only the required eigenstate

Ep = [ψ[n].data[i] + real(E[n]) for i in 1:length(ψn[1])]#Adding Energies of each level just for representation.
#Plotting the States
#plot!(ptsx, 500*real(ψ[n].data)) #Can be plotted from all eigenstates 
plot(ptsx, [5*real(Ep) , V], title="SHO Box n=$n eigenstate", label=["ψ(x)" "V(x)"])
xlabel!("x")
ylabel!("ψ(x)")
print("The Energy of the corresponding state is  : $(real(E[n])) ")
