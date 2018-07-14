function v=FracAm_Put_OS(alpha, N, M, K, T, r, sigma, S, g)
k=T/M;
%create a matrix to store value for each time step
val_matrix=zeros(N+1,M+1);
%the first column of the matrix is the option value at maturity time
val_matrix(:,1)=g;
%use corr to represent rou(alpha,k)
corr=1/(gamma(2-alpha)*(k^alpha));
w1=frac_w(alpha, 1);
psi=zeros(N+1,1);
%psi_0 is a zero vector with length N+1
%we calculate the parameter of Vj
%p0;
%we calculate the parameter of Vj prime
%p1;
%we calculate the parameter of Vj double prime
%p2;
for j=1:M
    %The first step of operator spliting
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
    loop_mult=fir_mult * val_matrix(:,j)+sec_mult+psi;
    Vj_tao=zeros(N+1,1);
    Vj_tao=fd2_bvpx(N,S(1),S(N+1),p0,p1,p2,1,0,K,1,0,0,loop_mult);
    %The second step of operator spliting
    for i=1:N+1
        if(psi(i)<=corr*w1*(Vj_tao(i)-g(i)))
            val_matrix(i,j+1)=Vj_tao(i)-psi(i)/(corr*w1);
            psi(i)=0;
        else
            val_matrix(i,j+1)=g(i);
            psi(i)=(g(i)-Vj_tao(i))*corr*w1+psi(i);
        end
    end
end
v=val_matrix(:,M+1);

%hold on
%plot(S,val_matrix(:,1));
%plot(S,val_matrix(:,M+1));
