function y_all = NLLF(k1, k2, k3, k4, M, w_h, w_a, zeta_h, zeta_a, A1, A2, A3, A4, H1, H2, H3, H4, psd_fL1, psd_fM1, psd_fL2, psd_fM2, psd_fML, sigma1, sigma2, psd_y_M, f)
y_all1 = 0;
y_all2 = 0;
for w_k = k1:k2
    psd_fL_temp = psd_fL1;
    psd_fM_temp = psd_fM1;
    H_omega_temp = H_omega(w_h, w_a, zeta_h, zeta_a, A1, A2, A3, A4, H1, H2, H3, H4, 2*pi*f(w_k));
    E_psd_y_temp = E_psd_y(H_omega_temp, psd_fL_temp, psd_fM_temp, psd_fML, sigma1, sigma2);
    y0 = real(E_psd_y_temp(1,1)*E_psd_y_temp(2,2))-real(E_psd_y_temp(2,1)*(E_psd_y_temp(1,2)));
    if y0 <= 0
        continue
    end
    y = M * log(y0) + M * real(trace(pinv(E_psd_y_temp) * psd_y_M{w_k}));
    y_all1 = y_all1 + y; 
end

for w_k = k3:k4
    psd_fL_temp = psd_fL2;
    psd_fM_temp = psd_fM2;
    H_omega_temp = H_omega(w_h, w_a, zeta_h, zeta_a, A1, A2, A3, A4, H1, H2, H3, H4, 2*pi*f(w_k));
    E_psd_y_temp = E_psd_y(H_omega_temp, psd_fL_temp, psd_fM_temp, psd_fML, sigma1, sigma2);
    y0 = real(E_psd_y_temp(1,1)*E_psd_y_temp(2,2))-real(E_psd_y_temp(2,1)*(E_psd_y_temp(1,2)));
    if y0 <= 0
        continue
    end
    y = M * log(y0) + M * real(trace(pinv(E_psd_y_temp) * psd_y_M{w_k}));
    y_all2 = y_all2 + y; 
end
y_all = y_all1 + y_all2;
end