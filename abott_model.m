% Model using Neuron class
clear;close all;

% time
tau = 0.2;  tspan = 0:tau:50; tspan_len = length(tspan);
T1=10;

% Izhikevich Neuron params
a=0.03; b=0.25; c=-60;  d=4;

% n - neuron
ai1 = Neuron2(1, tspan_len, b);
ai2 = Neuron2(2, tspan_len, b);
bi = Neuron2(3, tspan_len, b);
ab = Neuron2(4, tspan_len, b);

tau_synapse = 50;

% thalamus neurons
thalamus_firing_rate = 5;
thalamus_neurons = cell(2,1);
thalamus_neurons{1,1} = Neuron;
thalamus_neurons{1,1}.id = -1;
thalamus_neurons{1,1}.spikes(51:100) = generate_poisson_spike(thalamus_firing_rate, 50, tau);
[xe,~,~] = generate_xr_xe_xi_from_spike_train(thalamus_neurons{1,1}.spikes(51:100), tau);
thalamus_neurons{1,1}.xe = zeros(1, length(tspan));
thalamus_neurons{1,1}.xe(51:100) = xe;


thalamus_neurons{2,1} = Neuron;
thalamus_neurons{2,1}.id = -2;
thalamus_neurons{2,1}.spikes(151:200) = generate_poisson_spike(thalamus_firing_rate, 50, tau);
[xe,~,~] = generate_xr_xe_xi_from_spike_train(thalamus_neurons{2,1}.spikes(151:200), tau);
thalamus_neurons{2,1}.xe = zeros(1, length(tspan));
thalamus_neurons{2,1}.xe(151:200) = xe;


tau_re = 0.9; tau_ei = 5.3; tau_ir = 800;
inh_tau_re = 0.9; inh_tau_ei = 5.3; inh_tau_ir = 800;

for ti = 2:length(tspan)
    % inhibitory neuron
    % i = 20*thalamus_neurons{1,1}.xe(ti) + 20*thalamus_neurons{2,1}.xe(ti);
    % v_t_minus_1 = inhibitory_neuron.voltage(ti-1);
    % u_t_minus_1 = inhibitory_neuron.feedback_var(ti-1);
    % inhibitory_neuron = inhibitory_neuron.update_vars(v_t_minus_1, u_t_minus_1, i, ti, tau, inh_tau_re, inh_tau_ei, inh_tau_ir);

    i_thA_to_ai1 = 5*thalamus_neurons{1,1}.xe(ti);
    v_t_minus_1 = ai1.voltage(ti-1);
    u_t_minus_1 = ai1.feedback_var(ti-1);
    ai1 = ai1.update_vars(v_t_minus_1, u_t_minus_1, i_thA_to_ae, ti, tau, tau_re, tau_ei, tau_ir);

    i_ai1_to_ai2 = -5*ae.xe(ti);
    v_t_minus_1 = ai2.voltage(ti-1);
    u_t_minus_1 = ai2.feedback_var(ti-1);
    ai2 = ai2.update_vars(v_t_minus_1, u_t_minus_1, i_ai1_to_ai2, ti, tau, tau_re, tau_ei, tau_ir);

    i_thB_to_bi = 5*thalamus_neurons{2,1}.xe(ti);
    i_ai2_to_bi = -5*ai2.xe(ti);
    i_to_bi = i_thB_to_bi + i_ai2_to_bi;
    v_t_minus_1 = bi.voltage(ti-1);
    u_t_minus_1 = bi.feedback_var(ti-1);
    bi = bi.update_vars(v_t_minus_1, u_t_minus_1, i_to_bi, ti, tau, tau_re, tau_ei, tau_ir);
    
    i_bi_to_ab = -5*bi.xe(ti);
    v_t_minus_1 = ab.voltage(ti-1);
    u_t_minus_1 = ab.feedback_var(ti-1);
    ab = ab.update_vars(v_t_minus_1, u_t_minus_1, i_bi_to_ab, ti, tau, tau_re, tau_ei, tau_ir);


end



