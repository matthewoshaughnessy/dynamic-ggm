rng(0);
n = 100;
p = 15;
nedge = 15;
E = zeros(p);
iinit = randperm(p,4);
E(1,2) = 1;
E(2,3) = 1;
E(3,4) = 1;
E(4,1) = 1;
for i = 4:nedge
  d = sum(udiag(E(:,1:i-1))~=0)+1;
  pcs = cumsum(d / sum(d));
  r = rand;
  i2 = find(r < pcs, 1);
  E(i,i2) = 1;
end
G = E + E';
plot(graph(G),'layout','force'); grid on; drawnow;
%loglog(sort(sum(E),'descend'),'o');

D = diag(sum(G));
eta = 1.1;
L = eta*D - G;

Lambda = diag(diag(inv(L)));
Theta = sqrt(Lambda)*L*sqrt(Lambda);
Sigma = inv(Theta);
X = mvnrnd(zeros(n,p),Sigma).';
S = 1/n*X*X';

%%
Theta_n = Theta + 0.01*randn(size(Theta));
alpha = 0.1;  eis = 0.1*rand(1,p);
thetani1 = @(i) sum(abs(Theta(i,setdiff(1:p,i))));
thetanis = sum(abs(Theta - diag(diag(Theta))));
thetanis_n = sum(abs(Theta_n - diag(diag(Theta_n))));

l1 = alpha*sum(log(thetanis+eis)) - alpha*sum(log(thetanis_n+eis));
l2 = alpha*sum((thetanis+eis)./(thetanis_n+eis)-1);
lambda = zeros(p,p);
for i = 1:p
  for j = 1:p
    lambda(i,j) = alpha*(1/(thetanis_n(i)+eis(i)) + 1/(thetanis_n(j)+eis(j)));
  end
end
l3 = sum(sum( lambda.*abs(Theta) - diag(diag(lambda.*abs(Theta))) ));

%%
%% derivation of conditional independence in precision matrices
rng(0);
n = 16;
x = randn(n,1);
A = rand(n);
Si = 1/2*(A+A');
ij = randperm(n,2);
i = ij(1);
j = ij(2);
% --- compute line 11 ---
line11 = -1/2*x'*Si*x;

% --- compute line 14 ---
line14 = 0;
for g = 1:n
  for h = g+1:n
    line14 = line14 - Si(g,h)*x(g)*x(h);
  end
end
%for g = 1:n
%  line14 = line14 - 1/2*Si(g,g)*x(g)^2;
%end
% --- compute line 15 ---
line15 = 0;
for g = setdiff(1:n,[i j])
  for h = setdiff(g+1:n,[i j])
    line15 = line15 - Si(g,h)*x(g)*x(h);
  end
end
line15 = line15 - Si(i,j)*x(i)*x(j);
%for g = setdiff(1:n,[i j])
%  line15 = line15 - Si(g,g)*x(g)^2;
%end
%line15 = line15 - 1/2*Si(i,i)*x(i)^2;

%% derivation of (5) in Liu & Ihler 2011
clear;
p = 4;
alpha = 0.01;
e = 0.5*rand(p,1);
Theta = 2*rand(p) - 1;
Thetan = Theta + 0.05*rand(p);
% line 1
line1 = 0;
theta_ni_l1 = @(i) sum(abs(Theta(i,setdiff(1:p,i))));
thetan_ni_l1 = @(i) sum(abs(Thetan(i,setdiff(1:p,i))));
for i = 1:p
  line1 = line1 + (theta_ni_l1(i)+e(i))/(thetan_ni_l1(i)+e(i)) - 1;
end
line1 = alpha*line1;
% line 2
line2 = 0;
for i = 1:p
  line2 = line2 + (theta_ni_l1(i)-thetan_ni_l1(i))/(thetan_ni_l1(i)+e(i));
end
line2 = alpha*line2;
% line 3
line3 = 0;
for i = 1:p
  for j = setdiff(1:p,i)
    line3 = line3 + (abs(Theta(i,j))-abs(Thetan(i,j))) / (thetan_ni_l1(i)+e(i));
  end
end
line3 = alpha*line3;
% line 4
line4 = 0;
for i = 1:p
  for j = setdiff(1:p,i)
    line4 = line4 + abs(Theta(i,j))/(thetan_ni_l1(i)+e(i)) ...
      - abs(Thetan(i,j))/(thetan_ni_l1(i)+e(i));
  end
end
line4 = alpha*line4;
% line 5
line5 = 0;
for i = 1:p
  for j = setdiff(1:p,i)
    line5 = line5 + (1/(thetan_ni_l1(i)+e(i))+1/(thetan_ni_l1(j)+e(j)))*abs(Theta(i,j)) ...
      - 1/(thetan_ni_l1(j)+e(j))*abs(Theta(i,j)) ...
      - 1/(thetan_ni_l1(i)+e(i))*abs(Thetan(i,j));
  end
end
line5 = alpha*line5;


