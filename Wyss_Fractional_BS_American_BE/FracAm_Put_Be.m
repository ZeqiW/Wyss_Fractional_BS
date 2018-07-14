function v=FracAm_Put_Be(alpha, N, M, K, T, r, sigma, S, U)
%U is the initial value when t=0
k=T/M;
%create a matrix to store value for each time step
val_matrix=zeros(N+1,M+1);
%the first column of the matrix is the option value at maturity time
val_matrix(:,1)=U;
%use corr to represent rou(alpha,k)
corr=1/(gamma(2-alpha)*(k^alpha));
w1=frac_w(alpha, 1);
%we calculate the parameter of Vj
%p0;
%we calculate the parameter of Vj prime
%p1;
%we calculate the parameter of Vj double prime
%p2;
for j=1:M
    p0=corr*w1+r;
    p1=-r;
    p2=-(sigma^2)/2;
    common_mult=corr;
    fir_mult=common_mult*(frac_w(alpha, 1));
    sec_mult=0;
    for n=2:j
        sec_mult =sec_mult + frac_w(alpha, n)*(val_matrix(:,j-n+2)-val_matrix(:,j-n+1));
    end
    sec_mult = -common_mult * sec_mult;
    loop_mult=fir_mult * val_matrix(:,j)+sec_mult;
    val_matrix(:,j+1)=fd2_bvpx_am(N,S(1),S(N+1),p0,p1,p2,1,0,K,1,0,0,loop_mult,val_matrix(:,1));
end
v=val_matrix(:,M+1);
%fid = 1;
%pname = 'newpic2.eps';
%figure (fid)
%hold on
%plot(S,val_matrix(:,1));
%plot(S,val_matrix(:,M+1));
%legend('Maturity', 'Present');
%xlabel('Stock Price S');
%ylabel('Put Option Price');
%print('-f1', '-depsc2', pname);


