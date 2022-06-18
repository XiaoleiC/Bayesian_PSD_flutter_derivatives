function H_output = H_omega(w_h, w_a, zeta_h, zeta_a, A1, A2, A3, A4, H1, H2, H3, H4, w)
H11R = w_a^2 - A3 - w^2;
H11I = w*(2*zeta_a*w_a - A2);
H12R = H3;
H12I = H2*w;
H21R = A4;
H21I = A1*w;
H22R = w_h^2 - H4 - w^2;
H22I = w*(2*zeta_h*w_h - H1);
%% numerical solutions
Hc = (w_h^2 - H4 - w^2 + (2*zeta_h*w_h - H1)*1i*w)*(w_a^2 - A3 - w^2 + (2*zeta_a*w_a - A2)*1i*w) - (A4 + 1i*A1*w)*(H3 + 1i*H2*w);
H_output= 1/Hc*[H11R + 1i*H11I, H12R + 1i*H12I; H21R + 1i*H21I, H22R + 1i*H22I];
% HcR = real(Hc);
% HcI = imag(Hc);
%% theoretical solutions
% HcR = (w_h^2 - w^2 - H4)*(w_a^2 - w^2 - A3) - (2*zeta_h*w_h*w - w*H1)*(2*zeta_a*w_a*w - w*A2) - A4*H3 + w^2*A1*H2;
% HcI = (w_h^2 - w^2 - H4)*(2*zeta_a*w_a*w - w*A2) + (w_a^2 - w^2 - A3)*(2*zeta_h*w_h*w - w*H1) - (A4*H2 + A1*H3)*w;
% H_output = [((H11R*HcR + H11I*HcI) + (H11I*HcR - H11R*HcI)*1i)/(HcR^2 + HcI^2), ((H12R*HcR + H12I*HcI) + (H12I*HcR - H12R*HcI)*1i)/(HcR^2 + HcI^2);...
%     ((H21R*HcR + H21I*HcI) + (H21I*HcR - H21R*HcI)*1i)/(HcR^2 + HcI^2), ((H22R*HcR + H22I*HcI) + (H22I*HcR - H22R*HcI)*1i)/(HcR^2 + HcI^2)];
end