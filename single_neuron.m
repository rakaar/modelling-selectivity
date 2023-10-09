a=0.03; b=0.25; c=-60;  d=4;
V=-64;  u=b*V;

VV=[];  uu=[];
tau = 0.2;  tspan = 0:tau:200;
T1=20;


spikes = zeros(1, length(tspan));
t_counter = 0;
for t=tspan
    t_counter = t_counter + 1;

    if (t>T1) & (t < T1+20) 
        I=-10;
    else
        I=0;
    end;
    V = V + tau*(0.04*V^2+5*V+140-u+I);
    u = u + tau*a*(b*V-u);
    if V > 30
        disp('spike')
        VV(end+1)=30;
        V = c;
        u = u + d;
        spikes(t_counter) = 1;
    else
        VV(end+1)=V;
    end;
    uu(end+1)=u;
end;
plot(tspan,VV,[0 T1 T1 (T1+5) (T1+5) max(tspan)],-85+[0 0 -5 -5 0 0]);
axis([0 max(tspan) -90 30])
axis off;
title('(M) rebound spike');
