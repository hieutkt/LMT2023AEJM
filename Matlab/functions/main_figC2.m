function main_figC2(xblue,xred)

    alpha = 1/3;
    beta = 0.96;
    phi = 0.5;
    delta = 0.1;
    s_ss = 1;
    s_low = 1.2*s_ss;

    par.alpha = alpha;
    par.beta = beta;
    par.phi = phi;
    par.delta = delta;

    k_vec0= [2;2];

    system = @(k_vec) het_irr_system(k_vec, s_ss, par);

    k_vec = fsolve(system, k_vec0);

    % k1 is irreversible

    T = 10;

    k1_sim = zeros(T+1, 1);
    k2_sim = zeros(T+1, 1);
    k_sim = zeros(T+1, 1);
    y_sim = zeros(T+1, 1);

    k1_sim(1) = k_vec(1);
    k2_sim(1) = k_vec(2);
    k_sim(1) = k1_sim(1)^phi*k2_sim(1)^(1-phi);
    y_sim(1) = s_ss*k_sim(1)^alpha;
    s_sim = ones(T+1,1)*s_ss;
    s_sim(1:end) = s_low;

    mpk1_ss = y_sim(1)/k1_sim(1);
    mpk2_ss = y_sim(1)/k2_sim(1);


    for t = 1:T
        s_prime = s_sim(t+1);


       k_sim(t) = k1_sim(t)^phi*k2_sim(t)^(1-phi);
       y_sim(t) = s_sim(t)*k_sim(t)^alpha;


       system = @(k_vec) het_irr_system(k_vec, s_prime, par);

       k_vec = fsolve(system, k_vec0); 


       if k_vec(1) < (1-delta)*k1_sim(t)

           k1_sim(t+1) = (1-delta)*k1_sim(t);

           k1_prime = k1_sim(t+1);

           system = @(k_vec) het_irr_system_binding(k_vec, s_prime, k1_prime, par);

           k_vec = fsolve(system, k_vec0(2));

           k2_sim(t+1) = k_vec;
       else

           k1_sim(t+1) = k_vec(1);

           k2_sim(t+1) = k_vec(2);



       end


    end


    mpk1_sim = y_sim./k1_sim;
    mpk2_sim = y_sim./k2_sim;

    time =0:1:8;

    figure
    plot(time, ([mpk1_ss; mpk1_sim(1:end-3)])./mpk1_ss, 'Color', xblue,'linewidth',3)
    hold on
    plot(time, [mpk2_ss; mpk2_sim(1:end-3)]./mpk2_ss, '--', 'Color', xred,'linewidth',3)
    legend({'$MPK_1$', '$MPK_2$'},'FontSize',24,'interpreter','latex')
    hold off
    grid on
    xlabel('Time','FontSize',24,'interpreter','latex')
    ylabel('MPK','FontSize',24,'interpreter','latex')
    print('figures/fig_25_C2a','-depsc')
    

    %%
    s_low = 0.8*s_ss;

    par.alpha = alpha;
    par.beta = beta;
    par.phi = phi;
    par.delta = delta;

    k_vec0= [2;2];

    system = @(k_vec) het_irr_system(k_vec, s_ss, par);

    k_vec = fsolve(system, k_vec0);

    % k1 is irreversible

    T = 10;

    k1_sim = zeros(T+1, 1);
    k2_sim = zeros(T+1, 1);
    k_sim = zeros(T+1, 1);
    y_sim = zeros(T+1, 1);

    k1_sim(1) = k_vec(1);
    k2_sim(1) = k_vec(2);
    k_sim(1) = k1_sim(1)^phi*k2_sim(1)^(1-phi);
    y_sim(1) = s_ss*k_sim(1)^alpha;
    s_sim = ones(T+1,1)*s_ss;
    s_sim(1:end) = s_low;

    mpk1_ss = y_sim(1)/k1_sim(1);
    mpk2_ss = y_sim(1)/k2_sim(1);


    for t = 1:T
        s_prime = s_sim(t+1);


       k_sim(t) = k1_sim(t)^phi*k2_sim(t)^(1-phi);
       y_sim(t) = s_sim(t)*k_sim(t)^alpha;


       system = @(k_vec) het_irr_system(k_vec, s_prime, par);

       k_vec = fsolve(system, k_vec0); 


       if k_vec(1) < (1-delta)*k1_sim(t)

           k1_sim(t+1) = (1-delta)*k1_sim(t);

           k1_prime = k1_sim(t+1);

           system = @(k_vec) het_irr_system_binding(k_vec, s_prime, k1_prime, par);

           k_vec = fsolve(system, k_vec0(2));

           k2_sim(t+1) = k_vec;
       else

           k1_sim(t+1) = k_vec(1);

           k2_sim(t+1) = k_vec(2);



       end


    end


    mpk1_sim = y_sim./k1_sim;
    mpk2_sim = y_sim./k2_sim;

    time =0:1:8;

    figure
    plot(time, ([mpk1_ss; mpk1_sim(1:end-3)])./mpk1_ss, 'Color', xblue,'linewidth',3)
    hold on
    plot(time, [mpk2_ss; mpk2_sim(1:end-3)]./mpk2_ss, '--', 'Color', xred,'linewidth',3)
    legend({'$MPK_1$', '$MPK_2$'},'FontSize',24,'interpreter','latex')
    hold off
    grid on
    xlabel('Time','FontSize',24,'interpreter','latex')
    ylabel('MPK','FontSize',24,'interpreter','latex')
    print('figures/fig_25_C2b','-depsc')

end














