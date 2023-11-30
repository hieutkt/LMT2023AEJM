function [z, k1, k2] = het_irr_system(k_vec, s, par);

k1 = k_vec(1);
k2 = k_vec(2);

k = k1^par.phi*k2^(1-par.phi);

z(1) = 1 - (par.beta*(s*par.alpha*k^(par.alpha -1)*par.phi*k1^(par.phi - 1)*k2^(1-par.phi) + 1-par.delta ));

z(2) = 1 - (par.beta*(s*par.alpha*k^(par.alpha -1)*(1-par.phi)*k1^(par.phi)*k2^(-par.phi) + 1-par.delta ));

end

