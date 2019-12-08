
###################
# The Rössler System      #
# Controller type: RFPT  #
# MIMO System              #
##################

####################
# The equation of motion #
####################


using LinearAlgebra
using PyPlot



#####################
# Functions for FPT #
#####################

function KB(e_int, e, qN_p)
  q_pDes = Λ^2 * e_int + 2 * Λ * e + qN_p
  return q_pDes
end

# The definition of the function of the adaptive deformation
function G_MIMO(r_prev, f_prev, des_now, err_limit)
	Amatr_h = (f_prev - des_now)
	error_norm = norm(Amatr_h) # Frobenius norm
	if error_norm > err_limit
		e_direction = Amatr_h / error_norm
		B_factor= B * tanh(A * error_norm)
		G = (1 + B_factor) * r_prev + B_factor * K * e_direction
  else
    G = r_prev
  end # if
  return G
end

#################
# Control Parameters #
#################
K=1e5
B=-1
A=1.97e-5
Λ=3

#############################################
# Parameters to design a Nominal Trajectory #
#############################################
# Vector of parameters for x, y, z axes respectively
Ampl = [2, 3, 1]
ω = [0.5, 0.75, 1]

#################
# Time variable #
#################
# cycle time in seconds
δt=1e-3
# The lenght of the simulation in seconds
# short run
LONG=Int(5e4)
# long run
# LONG=Int(2.5e4)

##########################
# Exact model parameters #
##########################
#βₑ=8/3
#σₑ=5
#ρₑ=40
a1=0.2
b1=0.2
c1=4
############################
# Approx. model parameters #
############################
#βₐ=7/3
#σₐ=4
#ρₐ=36
a2=0.1
b2=0.1
c2=6.0
# Timestamps
t=zeros(LONG)
# Nominal trajectories
qN = zeros(Float64, LONG, 3)
# 1st derivatives of nominal trajectories
qN_p = zeros(Float64, LONG, 3)
# Adaptive control variables
# Desired responses
qDes_p = zeros(Float64, LONG, 3)
# Actual trajectories
qA = zeros(Float64, LONG, 3)
# Control forces
u = zeros(Float64, LONG, 3)
# Errors: integrated, actual, derivated
errors_int = zeros(Float64, LONG, 3)
errors = zeros(Float64, LONG, 3)
errors_p = zeros(Float64, LONG, 3)
# Responses from the previous step - 1st derivatives of actual trajectories
past_responses = zeros(Float64, LONG, 3)
# Inputs from the previous step
past_inputs = zeros(Float64, LONG, 3)
# Non-adaptive (PID) control variables
# Desired responses
qDes_p_pid = zeros(Float64, LONG, 3)
# Actual trajectories
qA_pid = zeros(Float64, LONG, 3)
# Control forces
u_pid = zeros(Float64, LONG, 3)
# Errors: integrated, actual, derivated
errors_int_pid = zeros(Float64, LONG, 3)
errors_pid = zeros(Float64, LONG, 3)
errors_p_pid = zeros(Float64, LONG, 3)
# Responses from the previous step - 1st derivatives of actual trajectories
past_responses_pid = zeros(Float64, LONG, 3)
# Inputs from the previous step
past_inputs_pid = zeros(Float64, LONG, 3)

#initial conditions
error_limit=1e-20
past_input = zeros(3)
past_response = zeros(3)
error_int = zeros(3)
past_input_pid = zeros(3)
past_response_pid = zeros(3)
error_int_pid = zeros(3)

#Simulation
for i=1:LONG-1
  t[i]=δt*i
  #Nominal Trajectories
  qN[i, :] = [Ampl[j] * sin(ω[j] * t[i]) for j = 1 : 3]
  # 1st derivatives of nominal trajectories
  qN_p[i, :] = [Ampl[j] * ω[j] * cos(ω[j] * t[i]) for j = 1 : 3]

  # Errors
  # Adaptive
  errors_p[i, :] = qN_p[i, :] - past_responses[i, :]
  errors[i, :] = qN[i, :] - qA[i, :]
  global error_int = error_int + δt * errors[i, :]
  errors_int[i, :] = error_int
  #Non-adaptive
  errors_p_pid[i, :] = qN_p[i, :] - past_responses_pid[i, :]
  errors_pid[i, :] = qN[i, :] - qA_pid[i, :]
  global error_int_pid = error_int_pid + δt * errors_pid[i, :]
  errors_int_pid[i, :] = error_int_pid

  qDes_p[i, :] = [KB(errors_int[i, j], errors[i, j], errors_p[i, j]) for j = 1 : 3]
  qDes_p_pid[i, :] = [KB(errors_int_pid[i, j], errors_pid[i, j], errors_p_pid[i, j]) for j = 1 : 3]

  # Deformation
  if i <= 10
    global past_input = qDes_p[i, :]
  else
    global past_input = G_MIMO(past_input, past_response, qDes_p[i, :], error_limit)
  end
  past_inputs[i, :] = past_input

  global past_input_pid = qDes_p_pid[i, :]
  past_inputs_pid[i, :] = past_input_pid

  x_G=past_input[1]
  y_G=past_input[2]
  z_G=past_input[3]

  x = qA[i, 1]
  y = qA[i, 2]
  z = qA[i, 3]

  u_x = x_G - (-y-x)
  u_y = y_G - (x+a1*y)
  u_z = z_G - (z*x+(x-c1)*z)
  u[i, :] = [u_x u_y u_z]

  x_p = (-y-x) + u_x
  y_p = (x+a2*y) + u_y
  z_p = (z*x+(x-c2)*z) + u_z
  global past_response = [x_p, y_p, z_p]
  past_responses[i, :] = past_response
  qA[i + 1, :] = qA[i, :] + δt * past_response

  x_G_pid=past_input_pid[1]
  y_G_pid=past_input_pid[2]
  z_G_pid=past_input_pid[3]

  x_pid = qA_pid[i, 1]
  y_pid = qA_pid[i, 2]
  z_pid = qA_pid[i, 3]

  u_x_pid = x_G_pid - (-y-x)
  u_y_pid = y_G_pid - x_pid * (b1 - z_pid) + y_pid
  u_z_pid = z_G_pid - x_pid * y_pid + b1 * z_pid
  u_pid[i, :] = [u_x_pid u_y_pid u_z_pid]

  x_p_pid = a2 * (y_pid - x_pid) + u_x_pid
  y_p_pid = x_pid * (b2 - z_pid) - y_pid + u_y_pid
  z_p_pid = x_pid * y_pid - b2 * z_pid + u_z_pid
  global past_response_pid = [x_p_pid, y_p_pid, z_p_pid]
  past_responses_pid[i, :] = past_response_pid
  qA_pid[i + 1, :] = qA_pid[i, :] + δt * past_response_pid
end #For

fig_caption = "nominal_realized_trajectories"
fig = figure(fig_caption)
grid("True")
title("Nominal and Realized Trajectories")
# x
# Create the 1st axis of a 3x1 array of axes
subplot(311)
ax1 = gca()
grid1 = grid("True")
ylabel(L"x")
plot(t[3 : LONG - 1], qN[3 : LONG - 1, 1], color="red", linewidth = 2, label="Nominal")
plot(t[3 : LONG - 1], qA[3 : LONG - 1, 1], color="green", linewidth = 2.5, label="Realized Adaptive", linestyle="--")
plot(t[3 : LONG - 1], qA_pid[3 : LONG - 1, 1], color="blue", label="Realized Non-adaptive", linestyle="-.")
# legend(loc="lower right",fancybox="True")
# legend()
# y
subplot(312,sharex=ax1)
ax2 = gca()
grid("True")
ylabel(L"y")
plot(t[3 : LONG - 1], qN[3 : LONG - 1, 2], color="red", linewidth = 2, label="Nominal")
plot(t[3 : LONG - 1], qA[3 : LONG - 1, 2], color="green", linewidth = 2.5, label="Realized Adaptive", linestyle="--")
plot(t[3 : LONG - 1], qA_pid[3 : LONG - 1, 2], color="blue", label="Realized Non-adaptive", linestyle="-.")
# legend(loc="upper right",fancybox="True")
# legend()
# z
subplot(313,sharex=ax2)
ax3 = gca()
grid("True")
xlabel(L"t, [$s$]")
ylabel(L"z")
plot(t[3 : LONG - 1], qN[3 : LONG - 1, 3], color="red", linewidth = 2, label="Nominal")
plot(t[3 : LONG - 1], qA[3 : LONG - 1, 3], color="green", linewidth = 2.5, label="Realized Adaptive", linestyle="--")
plot(t[3 : LONG - 1], qA_pid[3 : LONG - 1, 3], color="blue", label="Realized Non-adaptive", linestyle="-.")
# legend(loc="upper right",fancybox="True")
legend()
tight_layout()
subplots_adjust(hspace=0.0) # Set the vertical spacing between axes
fig[:canvas][:draw]()


fig_caption = "tracking_error"
fig = figure(fig_caption)
grid("True")
title("Tracking Errors vs Time")
# x
# Create the 1st axis of a 3x1 array of axes
subplot(311)
ax1 = gca()
grid1 = grid("True")
ylabel(L"$x^N-x$")
plot(t[3 : LONG - 1], errors[3 : LONG - 1, 1], color="red", linewidth=2, label="Adaptive")
plot(t[3 : LONG - 1], errors_pid[3 : LONG - 1, 1], color="green", linewidth=2.5, linestyle="--", label="Non-adaptive")
# legend(loc="upper left",fancybox="True")
# y
subplot(312,sharex=ax1)
ax2 = gca()
grid("True")
ylabel(L"$y^N-y$")
plot(t[3 : LONG - 1], errors[3 : LONG - 1, 2], color="red", linewidth=2, label="Adaptive")
plot(t[3 : LONG - 1], errors_pid[3 : LONG - 1, 2], color="green", linewidth=2.5, linestyle="--", label="Non-adaptive")
# legend(loc="lower left",fancybox="True")
# z
subplot(313,sharex=ax2)
ax3 = gca()
grid("True")
xlabel(L"t, [$s$]")
ylabel(L"$z^N-z$")
plot(t[3 : LONG - 1], errors[3 : LONG - 1, 3], color="red", linewidth=2, label="Adaptive")
plot(t[3 : LONG - 1], errors_pid[3 : LONG - 1, 3], color="green", linewidth=2.5, linestyle="--", label="Non-adaptive")
# legend(loc="lower left",fancybox="True")
legend()
tight_layout()
subplots_adjust(hspace=0.190) # Set the vertical spacing between axes
fig[:canvas][:draw]()


fig_caption = "phase_trajectories_x"
figure(fig_caption)
grid("True")
title("Phase Trajectories")
xlabel(L"x")
ylabel(L"$\dot{x}$")
plot(qN[3 : LONG - 1, 1], qN_p[3 : LONG - 1, 1], color="red", linewidth = 2, label="Nominal")
plot(qA[3 : LONG - 1, 1], past_responses[3 : LONG - 1, 1], color="green", linewidth = 2.5, linestyle="--", label="Realized Adaptive")
legend(loc="lower left",fancybox="True")
tight_layout()


fig_caption = "phase_trajectories_y"
figure(fig_caption)
grid("True")
title("Phase Trajectories")
xlabel(L"y")
ylabel(L"$\dot{y}$")
plot(qN[3 : LONG - 1, 2], qN_p[3 : LONG - 1, 2], color="red", linewidth = 2, label="Nominal")
plot(qA[3 : LONG - 1, 2], past_responses[3 : LONG - 1, 2], color="green", linewidth = 2.5, linestyle="--", label="Realized Adaptive")
legend(loc="upper left",fancybox="True")
tight_layout()


fig_caption = "phase_trajectories_z"
figure(fig_caption)
grid("True")
title("Phase Trajectories")
xlabel(L"z")
ylabel(L"$\dot{z}$")
plot(qN[3 : LONG - 1, 3], qN_p[3 : LONG - 1, 3], color="red", linewidth = 2, label="Nominal")
plot(qA[3 : LONG - 1, 3], past_responses[3 : LONG - 1, 3], color="green", linewidth = 2.5, linestyle="--", label="Realized Adaptive")
legend(loc="upper left",fancybox="True")
tight_layout()


fig_caption = "control_signal"
figure(fig_caption)
grid("True")
title("Control Signals vs Time")
xlabel(L"t, $[s]$")
ylabel(L"u")
plot(t[3 : LONG - 1], u[3 : LONG - 1, 1], color="red", label=L"$u_1$")
plot(t[3 : LONG - 1], u[3 : LONG - 1, 2], color="green", label=L"$u_1$")
plot(t[3 : LONG - 1], u[3 : LONG - 1, 3], color="blue", label=L"$u_3$")
legend(loc="lower left",fancybox="True")
tight_layout()


fig_caption = "velocities"
fig = figure(fig_caption)
grid("True")
title(L"$1^{st} Time Derivatives$")
# x
# Create the 1st axis of a 3x1 array of axes
subplot(311)
ax1 = gca()
grid1 = grid("True")
ylabel(L"$\dot{x}$")
plot(t[3 : LONG - 1], qN_p[3 : LONG - 1, 1], color="red", linewidth = 2, label=L"\dot{x}^{N}")
plot(t[3 : LONG - 1], past_responses[3 : LONG - 1, 1], color="green", linewidth = 3, label=L"\dot{x}", linestyle="--")
plot(t[3 : LONG - 1], past_inputs[3 : LONG - 1, 1], color="blue", label=L"\dot{x}^{Des}", linestyle="-.")
# legend(loc="lower left",fancybox="True")
# legend()
# y
subplot(312,sharex=ax1)
ax2 = gca()
grid("True")
ylabel(L"$\dot{y}$")
plot(t[3 : LONG - 1], qN_p[3 : LONG - 1, 2], color="red", linewidth = 2, label=L"\dot{y}^{N}")
plot(t[3 : LONG - 1], past_responses[3 : LONG - 1, 2], color="green", linewidth = 3, label=L"\dot{y}", linestyle="--")
plot(t[3 : LONG - 1], past_inputs[3 : LONG - 1, 2], color="blue", label=L"\dot{y}^{Des}", linestyle="-.")
# legend(loc="lower left",fancybox="True")
# legend()
# z
subplot(313,sharex=ax2)
ax3 = gca()
grid("True")
xlabel("t, [s]")
ylabel(L"$\dot{z}$")
plot(t[3 : LONG - 1], qN_p[3 : LONG - 1, 3], color="red", linewidth = 2, label="Nominal")
plot(t[3 : LONG - 1], past_responses[3 : LONG - 1, 3], color="green", linewidth = 3, label="Realized", linestyle="--")
plot(t[3 : LONG - 1], past_inputs[3 : LONG - 1, 3], color="blue", label="Desired", linestyle="-.")
# legend(loc="lower left",fancybox="True")
legend()
tight_layout()
subplots_adjust(hspace=0.0) # Set the vertical spacing between axes
fig[:canvas][:draw]()


# Show plots in noninteractive mode
show()

