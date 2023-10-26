% Model using Neuron class
clear;

% time
tau = 0.2;  tspan = 0:tau:50;
T1=20;

% Izhikevich Neuron params
a=0.03; b=0.25; c=-60;  d=4;

% n - neuron
n = Neuron;
n.id  = 1;
n.voltage = zeros(1, length(tspan));
n.spikes = zeros(1, length(tspan));
n.feedback_var = zeros(1, length(tspan));

n.voltage(1) = -64;
n.feedback_var(1) = b*n.voltage(1);


% input current
input_current = zeros(1, length(tspan));
input_current(T1:T1 + round(20/tau)) = 10;

for ti = 2:length(tspan)
    i = input_current(ti);

    % update voltage
    v_t_minus_1 = n.voltage(ti-1);
    u_t_minus_1 = n.feedback_var(ti-1);

    n = n.update_vars(v_t_minus_1, u_t_minus_1, i, ti);
end


figure
subplot(1,2,1)
plot(n.spikes)

subplot(1,2,2)
plot(input_current)

