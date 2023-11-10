% Single Izhikevich neuron class
classdef Neuron2
    
    properties
        id
        voltage
        feedback_var
        spikes 
        pre_synaptic_current
        input_current
        xe
        xr
        xi
        psc
    end

    methods
        function obj = Neuron2(id, tspan_len, b)
            obj.id = id;
            obj.voltage = zeros(1, tspan_len);
            obj.spikes = zeros(1, tspan_len);
            obj.feedback_var = zeros(1, tspan_len);
            
            obj.voltage(1) = -64;
            obj.feedback_var(1) = b*obj.voltage(1);
            
            obj.xr = zeros(1, tspan_len);
            obj.xe = zeros(1, tspan_len);
            obj.xi = zeros(1, tspan_len);

            obj.xr(1) = 1;
        end
        
        function obj = update_vars(obj, v_t_minus_1, u_t_minus_1, i, ti, tau, tau_re, tau_ei, tau_ir)

            obj.input_current(ti) = i;
            
            a=0.03; b=0.25; c=-60;  d=4;
           

            v_t  = v_t_minus_1 + tau*(0.04*v_t_minus_1^2 + 5*v_t_minus_1 + 140 - u_t_minus_1 + i);
            u_t = u_t_minus_1  + tau*a*(b*v_t - u_t_minus_1);
            
            if v_t > 30
                disp('spike')
                obj.voltage(ti) = c;
                obj.feedback_var(ti) = u_t + d;
                obj.spikes(ti) = 1;
            else
                obj.voltage(ti) = v_t;
                obj.feedback_var(ti) = u_t;
            end

            % synapse
            M = obj.spikes(ti);
            obj.xr(ti) = obj.xr(ti-1) + tau*(-M*obj.xr(ti-1)/tau_re + obj.xi(ti-1)/tau_ir);
            obj.xe(ti) = obj.xe(ti-1) + tau*(M*obj.xr(ti-1)/tau_re - obj.xe(ti-1)/tau_ei);
            obj.xi(ti) = obj.xi(ti-1) + tau*(obj.xe(ti-1)/tau_ei - obj.xi(ti-1)/tau_ir);
            

        end % update vars
    end

end