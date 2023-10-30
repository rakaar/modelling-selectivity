function [xr,xe,xi] = generate_xr_xe_xi_from_spike_train(spike_train, tau)
    spike_train_len = length(spike_train);
    xr = zeros(1, spike_train_len);
    xe = zeros(1, spike_train_len);
    xi = zeros(1, spike_train_len);
    

    xr(1) = 1;

    tau_re = 2; tau_ei = 10; tau_ir = 50;


    for t = 2:spike_train_len
       M = spike_train(t);
       xr(t) = xr(t-1) + tau*(-M*xr(t-1)/tau_re + xi(t-1)/tau_ir);
       xe(t) = xe(t-1) + tau*(M*xr(t-1)/tau_re - xe(t-1)/tau_ei);
       xi(t) = xi(t-1) + tau*(xe(t-1)/tau_ei - xi(t-1)/tau_ir); 
    end
end