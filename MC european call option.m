clear all
close all

S_0=35;
K=35;
sigma=0.2;
r=0.04;
T=0.5;
N=10^5;  
C_sim=zeros(N,1);
D_sim1=zeros(N,1);
D_sim2=zeros(N,1);
ds=0.1;              % Delta S, used to calculate delta
for i=1:N
    W_T=normrnd(0,1);
    S_T=S_0*exp((r-0.5*sigma^2)*T+sqrt(T)*sigma*W_T);
    S_placeholder=exp((r-0.5*sigma^2)*T+sqrt(T)*sigma*W_T);
    C_sim(i,1)=max([S_T-K 0]);
    D_sim1(i,1)=(max([(S_0+ds)*S_placeholder-K 0])-max([S_0*S_placeholder-K 0]))/ds;           %Calculating delta using forward difference approximation
    D_sim2(i,1)=(max([(S_0+ds)*S_placeholder-K 0])-max([(S_0-ds)*S_placeholder-K 0]))/(2*ds);  %Calculating delta using central difference approximation
end
D_sim1=exp(-r*T)*D_sim1;
D_simu1=mean(D_sim1)                               %Monte carlo estimate of the Delta using forward difference approximation
D_sim2=exp(-r*T)*D_sim2;
D_simu2=mean(D_sim2)                               %Monte carlo estimate of the Delta using central difference approximation
C_simu=exp(-r*T)*mean(C_sim);
C_sim=exp(-r*T)*C_sim;                             %Monte carlo estimate of the call option
d_1=(log(S_0/K)+(r+0.5*sigma^2)*T)/(sigma*sqrt(T));
d_2=d_1-sigma*sqrt(T);
normcdf(d_1)
C_exact=S_0*normcdf(d_1)-exp(-r*T)*K*normcdf(d_2); % Theoretical price of call option
std_C=std(C_sim);                                  % Standard deviation of call option
std_D1=std(D_sim1)                                 % Standard deviation of Delta
std_D2=std(D_sim2)                                 %           -||-
e_C=std_C/sqrt(N);                                 % Standard error of call option
e_D1=std_D1/sqrt(N)                                % Standard error of Delta
e_D2=std_D2/sqrt(N)                                %           -||-
abs_error_C=abs(C_simu-C_exact);                   % Absolute value of difference between thoeretical value and simulated
abs_error_D1=abs(D_simu1-normcdf(d_1))             %           -||-
abs_error_D2=abs(D_simu2-normcdf(d_1))             %           -||-

