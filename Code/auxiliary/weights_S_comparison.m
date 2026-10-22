close all
manual10 = 1/((s/M_lat + 10)/(s + A_lat*10));
manual8 = 1/((s/M_lat + wc_lat)/(s + A_lat*wc_lat));

 mkw10 = (makeweight(A_lat, 10, M_lat, 0 , 1));
 mkw8 = (makeweight(A_lat, wc_lat, M_lat, 0 , 1));

 hold on
 grid on


 bode(manual10)
 bode(manual8)

 bode(mkw10)
 bode(mkw8)

 bode(S_2_phi)

 legend('man10', 'man_7', 'mkw10', 'mkw_7', 'S2')