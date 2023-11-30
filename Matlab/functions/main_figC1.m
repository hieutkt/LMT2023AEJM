function main_figC1(xblue,xred)


    phi_vec = [0.8 0.2];


    for j = 1:2

        alpha = 1/3;
        beta = 0.96;
        phi = phi_vec(j);
        delta1 = 0.1;
        delta2 = 1;
        s_ss = 1;
        s_low = 1.2*s_ss;

        par.alpha = alpha;
        par.beta = beta;
        par.phi = phi;
        par.delta1 = delta1;
        par.delta2 = delta2;

        k_vec0= [2;2];

        system = @(k_vec) het_depr_system(k_vec, s_ss, par);

        k_vec = fsolve(system, k_vec0);

        % k1 and k2 are irreversible, but k2 fully depreciates so irrelevant

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

        y_ss = y_sim(1);
        k_ss = k_sim(1);
        mpk_ss = y_sim(1)/k_sim(1);
        mpk1_ss = y_sim(1)/k1_sim(1);
        mpk2_ss = y_sim(1)/k2_sim(1);


        for t = 1:T
            s_prime = s_sim(t+1);

           k_sim(t) = k1_sim(t)^phi*k2_sim(t)^(1-phi);
           y_sim(t) = s_sim(t)*k_sim(t)^alpha;

           system = @(k_vec) het_depr_system(k_vec, s_prime, par);

           k_vec = fsolve(system, k_vec0); 


           if k_vec(1) < (1-delta1)*k1_sim(t)

               k1_sim(t+1) = (1-delta1)*k1_sim(t);

               k1_prime = k1_sim(t+1);

               system = @(k_vec) het_depr_system_binding(k_vec, s_prime, k1_prime, par);

               k_vec = fsolve(system, k_vec0(2));

               k2_sim(t+1) = k_vec;
           else

               k1_sim(t+1) = k_vec(1);

               k2_sim(t+1) = k_vec(2);



           end


        end

        mpk_sim = y_sim./k_sim;
        mpk1_sim = y_sim./k1_sim;
        mpk2_sim = y_sim./k2_sim;

        time =0:1:8;



        mpk_pic(:,j) = [mpk_ss ; mpk_sim]./mpk_ss;
    end



    figure
    % subplot(1,2,1)
    plot(time, mpk_pic(1:end-3,1), 'Color', xblue,'linewidth',3)
    hold on
    plot(time, mpk_pic(1:end-3,2), '--', 'Color', xred,'linewidth',3)
    legend({'$\phi = 0.8$', '$\phi = 0.2$'},'FontSize',24,'interpreter','latex')
    grid on
    xlabel('Time','FontSize',24,'interpreter','latex')
    ylabel('MPK','FontSize',24,'interpreter','latex')
    print('figures/fig_24_C1a','-depsc')


    %%

    for j = 1:2
        % 
        % alpha = 1/3;
        % beta = 0.96;
        phi = phi_vec(j);
        % delta1 = 0.1;
        % delta2 = 1;
        % s_ss = 1;
        s_low = 0.8*s_ss;

        par.alpha = alpha;
        par.beta = beta;
        par.phi = phi;
        par.delta1 = delta1;
        par.delta2 = delta2;

        k_vec0= [2;2];

        system = @(k_vec) het_depr_system(k_vec, s_ss, par);

        k_vec = fsolve(system, k_vec0);

        % k1 and k2 are irreversible, but k2 fully depreciates so irrelevant

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

        y_ss = y_sim(1);
        k_ss = k_sim(1);
        mpk_ss = y_sim(1)/k_sim(1);
        mpk1_ss = y_sim(1)/k1_sim(1);
        mpk2_ss = y_sim(1)/k2_sim(1);


        for t = 1:T
            s_prime = s_sim(t+1);

           k_sim(t) = k1_sim(t)^phi*k2_sim(t)^(1-phi);
           y_sim(t) = s_sim(t)*k_sim(t)^alpha;

           system = @(k_vec) het_depr_system(k_vec, s_prime, par);

           k_vec = fsolve(system, k_vec0); 


           if k_vec(1) < (1-delta1)*k1_sim(t)

               k1_sim(t+1) = (1-delta1)*k1_sim(t);

               k1_prime = k1_sim(t+1);

               system = @(k_vec) het_depr_system_binding(k_vec, s_prime, k1_prime, par);

               k_vec = fsolve(system, k_vec0(2));

               k2_sim(t+1) = k_vec;
           else

               k1_sim(t+1) = k_vec(1);

               k2_sim(t+1) = k_vec(2);



           end


        end

        mpk_sim = y_sim./k_sim;
        mpk1_sim = y_sim./k1_sim;
        mpk2_sim = y_sim./k2_sim;

        time =0:1:8;



        mpk_pic(:,j) = [mpk_ss ; mpk_sim]./mpk_ss;
    end


    figure
    % subplot(1,2,2)
    plot(time, mpk_pic(1:end-3,1), 'Color', xblue,'linewidth',3)
    hold on
    plot(time, mpk_pic(1:end-3,2), '--', 'Color', xred,'linewidth',3)
    legend({'$\phi = 0.8$', '$\phi = 0.2$'},'FontSize',24,'interpreter','latex')
    grid on
    xlabel('Time','FontSize',24,'interpreter','latex')
    ylabel('MPK','FontSize',24,'interpreter','latex')
    print('figures/fig_24_C1b','-depsc')


end











