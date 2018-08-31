clear all
close all
clc

errvec = zeros(991,1);
hvec = zeros(991,1);
h2vec = zeros(991,1);

for n = 10:1000
    %n = 50
    h = 2/n;
    tend = 2;
    a = 1;
    x = -1+h*(1:n)';
    u0 = exp(-36*x.^2);
    u = u0;
    up = zeros(size(u));
    k = 0.5*h;
    %k = 0.5*h^2;
    nt = ceil(tend/k);
    k = tend/nt; 
    lam = k/h;
    for it = 1:nt
        for i = 1:n
            ip = i+1;
            im = i-1;
            if (ip == n+1) 
                ip = 1;
            end
            if (im == 0) 
                im = n; 
            end                    
            %p = randi(10);
            %if (p<=5)
                up(i) = u(i) - a*lam*(u(i)-u(im)); %upwind
            %end
            %if (p>5)
                %up(i) = u(i) - 0.5*a*lam*(u(ip)-u(im)); %central
            %end    
        end
        u = up;
        %figure(1)
        %plot(x,u,'ro',x,u0,'k')
        %axis([-1 1 -0.1 1.1])
        %drawnow
    end
    errvec(n-9)= max(abs(u-u0));
    hvec(n-9) = h;
    h2vec(n-9) = h^2;
    n
end

figure(1)
loglog(hvec,errvec,'r',hvec,hvec,'b',hvec,h2vec,'k','linewidth',2)
title('Error vs h')
ylabel('ERROR')
xlabel('h')
legend ('Error','h','h^2')
