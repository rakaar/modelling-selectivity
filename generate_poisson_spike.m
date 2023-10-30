function poisson_spike = generate_poisson_spike(rate, tspan_len, dt)
        poisson_spike = zeros(1,tspan_len);
        spk_times = [];
        current_time = 0;
        duration = tspan_len/dt;
        while current_time < duration
            isi = -log(1 - rand)/(rate*dt);
            if ceil(isi) == 1
                continue
            end
            current_time = current_time + isi;
            if current_time <= tspan_len
                spk_times = [spk_times current_time];
            end 
        end

        spk_indices = ceil(spk_times);
        poisson_spike(spk_indices) = 1;
end