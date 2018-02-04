%This code runs a simple regression model
%The model is used for expoloratory purposes
%The aim is to see if there is a correaltion between States and number
%of arms recovered
%The data used is from
%https://www.atf.gov/docs/thefthitsontracesbyrecstn-2014163882xlsx/download

clc
clear
format bank
format compact

load firearmdata data
[n k]=size(data);
t=0;
sigma=1;
scale=1;
gibbs=1000;
x=normrnd(0,scale,n,k);
u=normrnd(0,sigma,n,1);
beta=normrnd(0,scale,k,1);
betaset=100;

y=x*beta+u;

tic
for g=1:gibbs;
    t=t+1;
    sigma=sqrt(((y-x*beta)'*(y-x*beta))/(sum(normrnd(0,1,n,1).^2)));
    beta=inv(x'*x)*(x'*y)+chol(sigma^2*inv(x'*x))'*normrnd(0,1,k,1);
    sigmas(t,:)=sigma'
    betas(t,:)=beta'
    where=[t];
end
toc

mean=[mean(betas(t,:))' mean(sigmas(t,:))'];
stddevation=[std(betas(t,:))' std(sigmas(t,:))'];

figure(1)
hist(betas, 100)
title('histogram of beta parameter')

figure(2)
hist(sigmas, 100)
title('histogram of sigma parameter')

