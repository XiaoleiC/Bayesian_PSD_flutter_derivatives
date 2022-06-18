function y = E_psd_y(H_omega_temp, psd_fL, psd_fM, psd_fML, sigma1, sigma2)
y = H_omega_temp * [psd_fL, psd_fML; psd_fML, psd_fM] * H_omega_temp' + diag([sigma1, sigma2]);
y(1,2) = 0;
y(2,1) = 0;
y = real(y);
end