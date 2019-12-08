
# Lorenz system
# x_p=σ(y-x)
# y_p=x(ρ-z)-y
# z_p=xy-βz
# Parameters
# σ=5
# β=8/3
# ρ=40
# Initial condition
# x₀=0.1
# y₀=0.1
# z₀=0.1

# Parameter setup
#σ=5
#β=8/3
#ρ=40
a=0.2
b=0.2
c=5.7
# Time setup
LONG=Int(1e5) #step number
δt=1e-3 # step length in seconds
l=LONG-1

# Arrays
x1=zeros(LONG)
y1=zeros(LONG)
z1=zeros(LONG)

x2=zeros(LONG)
y2=zeros(LONG)
z2=zeros(LONG)

t=zeros(LONG)

# Initial conditions
x1[1]=1
y1[1]=0
z1[1]=1

x2[1]=1.4
y2[1]=0.4
z2[1]=1.4

# Simulation
for i=1:l
	# time
	t[i]=i*δt
	# equations of motion
	#x1_p=σ*(y1[i]-x1[i])
	#y1_p=x1[i]*(ρ-z1[i])-y1[i]
	#z1_p=x1[i]*y1[i]-β*z1[i]
	x1_p=-y1[i]-z1[i]
	y1_p=x1[i]+(a*y1[i])
	z1_p=b+(z1[i]*(x1[i]-c))
	#x2_p=σ*(y2[i]-x2[i])
	#y2_p=x2[i]*(ρ-z2[i])-y2[i]
	#z2_p=x2[i]*y2[i]-β*z2[i]
    x2_p=-y2[i]-x2[i]
	y2_p=x2[i]+(a*y2[i])
	z2_p=b+(z2[i]*(x2[i]-c))

	# Euler Integral
	x1[i+1]=x1[i]+δt*x1_p
	y1[i+1]=y1[i]+δt*y1_p
	z1[i+1]=z1[i]+δt*z1_p
	x2[i+1]=x2[i]+δt*x2_p
	y2[i+1]=y2[i]+δt*y2_p
	z2[i+1]=z2[i]+δt*z2_p
end #for

# Plotting
using PyPlot
title("Roessler System")
plot3D(x1[1:l],y1[1:l],z1[1:l])
plot3D(x2[1:l],y2[1:l],z2[1:l])




