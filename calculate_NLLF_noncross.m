function y = calculate_NLLF_noncross(k1, k2, k3, k4, M, w_h, w_a, zeta_h, zeta_a, theta, psd_y_M, f)
A1 = theta(1);
A2 = theta(2);
A3 = theta(3);
A4 = theta(4);
H1 = theta(5);
H2 = theta(6);
H3 = theta(7);
H4 = theta(8);
psd_fL1 = 10^theta(9);
psd_fM1 = 10^theta(10);
psd_fL2 = 10^theta(11);
psd_fM2 = 10^theta(12);
sigma1 = 0;
sigma2 = 0;
y = NLLF(k1, k2, k3, k4, M, w_h, w_a, zeta_h, zeta_a, A1, A2, A3, A4, H1, H2, H3, H4, psd_fL1, psd_fM1, psd_fL2, psd_fM2, 0, sigma1, sigma2, psd_y_M, f);
end