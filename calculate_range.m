% TODO : Find the Bsweep of chirp for 1 m resolution
d_res = 1;
c =  3*10^8;
B = c/(2* d_res);
d_max = 300;

% TODO : Calculate the chirp time based on the Radar's Max Range
%Ts = (2*d_max*B)/(f_s*c);    // another formula;

T_s = (5.5*2*d_max)/c;

IF_array =[0 , 1.1 , 13 , 24 ]*10^6;
calculated_range = c*IF_array*T_s/(2*B)


% TODO : define the frequency shifts 


% Display the calculated range
disp(calculated_range);