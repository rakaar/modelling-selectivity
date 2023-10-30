% Model using Neuron class
clear;

% time
tau = 0.2;  tspan = 0:tau:50;
T1=10;

% Izhikevich Neuron params
a=0.03; b=0.25; c=-60;  d=4;

% n - neuron

neuron_arr = cell(2, 1);
for i = 1:2
    neuron_arr{i,1} = Neuron;
    neuron_arr{i,1}.id  = 1;
    neuron_arr{i,1}.voltage = zeros(1, length(tspan));
    neuron_arr{i,1}.spikes = zeros(1, length(tspan));
    neuron_arr{i,1}.feedback_var = zeros(1, length(tspan));
    
    neuron_arr{i,1}.voltage(1) = -64;
    neuron_arr{i,1}.feedback_var(1) = b*neuron_arr{i,1}.voltage(1);
    
    neuron_arr{i,1}.xr = zeros(1, length(tspan));
    neuron_arr{i,1}.xe = zeros(1, length(tspan));
    neuron_arr{i,1}.xi = zeros(1, length(tspan));
    
    neuron_arr{i,1}.xr(1) = 1;
    
end
% weights initialization
synaptic_weights = zeros(2,2, length(tspan));
% 1 -> 2: pre, post
synaptic_weights(1,2,1) = 1; 
tau_synapse = 50;

% input current
input_current = zeros(1, length(tspan));
input_current(T1:T1 + round(30/tau)) = 10;


% thalamus neurons
thalamus_firing_rate = 5;
thalamus_neurons = cell(2,1);
thalamus_neurons{1,1} = Neuron;
thalamus_neurons{1,1}.id = -1;
thalamus_neurons{1,1}.spikes(51:100) = generate_poisson_spike(thalamus_firing_rate, 50, tau);
[xe,~,~] = generate_xr_xe_xi_from_spike_train(thalamus_neurons{1,1}.spikes(51:100), tau);
thalamus_neurons{1,1}.xe(51:100) = xe;
figure
plot(thalamus_neurons{1,1}.xe, 'LineWidth', 4)
hold on
plot(thalamus_neurons{1,1}.spikes)

thalamus_neurons{2,1} = Neuron;
thalamus_neurons{2,1}.id = -2;
thalamus_neurons{2,1}.spikes(151:200) = generate_poisson_spike(thalamus_firing_rate, 50, tau);
[xe,~,~] = generate_xr_xe_xi_from_spike_train(thalamus_neurons{2,1}.spikes(151:200), tau);
thalamus_neurons{2,1}.xe(151:200) = xe;
figure
plot(thalamus_neurons{2,1}.xe, 'LineWidth', 4)
hold on
plot(thalamus_neurons{2,1}.spikes)

for ti = 2:length(tspan)
    

    for n = 1:2
        if n == 1
            i = input_current(ti);

            v_t_minus_1 = neuron_arr{n,1}.voltage(ti-1);
            u_t_minus_1 = neuron_arr{n,1}.feedback_var(ti-1);
    
            neuron_arr{n,1} = neuron_arr{n,1}.update_vars(v_t_minus_1, u_t_minus_1, i, ti, tau);
        else
            i = synaptic_weights(1,2)*5*neuron_arr{1,1}.xe(ti);
            v_t_minus_1 = neuron_arr{n,1}.voltage(ti-1);
            u_t_minus_1 = neuron_arr{n,1}.feedback_var(ti-1);
    
            neuron_arr{n,1} = neuron_arr{n,1}.update_vars(v_t_minus_1, u_t_minus_1, i, ti, tau);
        end
    end

    % update synaptic weights
    pre_syn_activity = sum(neuron_arr{1,1}.spikes(max(1,ti-50):ti-1));
    post_syn_activity = sum(neuron_arr{2,1}.spikes(max(1,ti-50):ti-1));
    synaptic_weights(1,2,ti) = synaptic_weights(1,2,ti - 1) +   tau*((pre_syn_activity*post_syn_activity)/tau_synapse);
    
end



