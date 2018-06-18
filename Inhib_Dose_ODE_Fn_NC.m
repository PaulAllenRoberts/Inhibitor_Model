function dy = Inhib_Dose_ODE_Fn_NC(~,y,p)

B_Free   = y(1);
B_Bound  = y(2);
A_Free   = y(3);
A_Bound  = y(4);

E        = p.E_init - p.phi_Bac*B_Bound - p.phi_A*A_Bound;

eta = (p.eta_max*E^p.n)/(p.gamma^p.n + E^p.n);

psi_Bac = 0;
psi_A   = 0;

dy = zeros(4,1);

dy(1) = p.r_F*B_Free*(1 - B_Free/p.K_F) + (1 - eta)*(p.r_B/p.h)*B_Bound*(1 - B_Bound/p.K_B)*heaviside(p.K_B - B_Bound)...
    - psi_Bac*B_Free - p.alpha_Bac*p.A_r*B_Free*E + (p.beta_Bac/p.h)*B_Bound;

dy(2) = (1 + (eta - 1)*heaviside(p.K_B - B_Bound))*p.r_B*B_Bound*(1 - B_Bound/p.K_B) + p.alpha_Bac*p.V*B_Free*E...
    - p.beta_Bac*B_Bound - p.delta_B*B_Bound;

dy(3) = -p.alpha_A*p.A_r*A_Free*E + (p.beta_A/p.h)*A_Bound - psi_A*A_Free;

dy(4) = p.alpha_A*p.V*A_Free*E - p.beta_A*A_Bound;