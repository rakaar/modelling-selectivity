% Single Izhikevich neuron class
classdef Neuron
    
    properties
        id
        voltage
        feedback_var
        spikes 
        pre_synaptic_current
        input_current
    end

    methods
        function obj = update_vars(obj, v_t_minus_1, u_t_minus_1, i, ti)
            a=0.03; b=0.25; c=-60;  d=4; tau = 0.2;

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

        end
    end

end