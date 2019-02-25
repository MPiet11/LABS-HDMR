function [mask_pos mask_neg]=mask_gen_2(x0,x1,c,T,p0,p1)

p0=log(p0);
p1=log(p1);

c=1;
Tprob=0.01;

TF=size(x0,2);
n0=size(x0,1);
n1=size(x1,1);


lp0m0=repmat(p0(1,:),TF,1);
lp0m1=repmat(p0(2,:),TF,1);
lp1m0=repmat(p1(1,:),TF,1);
lp1m1=repmat(p1(2,:),TF,1);


lpj0=cat(3,lp0m0,lp0m0,lp0m1,lp0m1);
lpj1=cat(3,lp1m0,lp1m0,lp1m1,lp1m1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%2 SNPs
%c=10;

x0c=1-x0;
x1c=1-x1;

%%%%class 0 a's
aj0i0s0=c+(x0c'*x0c)+c*eye(TF);
aj1i0s0=c+(x0c'*x0)+c*eye(TF);

aj0i1s0=c+(x0'*x0c)+c*eye(TF);
aj1i1s0=c+(x0'*x0)+c*eye(TF);

%%%%class 1 a's
aj0i0s1=c+(x1c'*x1c)+c*eye(TF);
aj1i0s1=c+(x1c'*x1)+c*eye(TF);

aj0i1s1=c+(x1'*x1c)+c*eye(TF);
aj1i1s1=c+(x1'*x1)+c*eye(TF);

%%%%% class 0 ep's
%epj0i0s0=aj0i0s0./(aj0i0s0+aj1i0s0);
%epj1i0s0=aj1i0s0./(aj0i0s0+aj1i0s0);


%epj0i1s0=aj0i1s0./(aj0i1s0+aj1i1s0);
%epj1i1s0=aj1i1s0./(aj0i1s0+aj1i1s0);

epj0i0s0=aj0i0s0./(n0+4*c);
epj1i0s0=aj1i0s0./(n0+4*c);


epj0i1s0=aj0i1s0./(n0+4*c);
epj1i1s0=aj1i1s0./(n0+4*c);

epj0i0s0=epj0i0s0.*(epj0i0s0>Tprob);
epj1i0s0=epj1i0s0.*(epj1i0s0>Tprob);

epj0i1s0=epj0i1s0.*(epj0i1s0>Tprob);
epj1i1s0=epj1i1s0.*(epj1i1s0>Tprob);


%%%%% class 1 ep's

%epj0i0s1=aj0i0s1./(aj0i0s1+aj1i0s1);
%epj1i0s1=aj1i0s1./(aj0i0s1+aj1i0s1);


%epj0i1s1=aj0i1s1./(aj0i1s1+aj1i1s1);
%epj1i1s1=aj1i1s1./(aj0i1s1+aj1i1s1);

epj0i0s1=aj0i0s1./(n1+4*c);
epj1i0s1=aj1i0s1./(n1+4*c);


epj0i1s1=aj0i1s1./(n1+4*c);
epj1i1s1=aj1i1s1./(n1+4*c);

epj0i0s1=epj0i0s1.*(epj0i0s1>Tprob);
epj1i0s1=epj1i0s1.*(epj1i0s1>Tprob);

epj0i1s1=epj0i1s1.*(epj0i1s1>Tprob);
epj1i1s1=epj1i1s1.*(epj1i1s1>Tprob);



eprj0i0=epj0i0s1./epj0i0s0;
eprj1i0=epj1i0s1./epj1i0s0;

eprj0i1=epj0i1s1./epj0i1s0;
eprj1i1=epj1i1s1./epj1i1s0;

for k=1:TF
    
    eprj0i0(k,k)=1;
    eprj1i0(k,k)=1;
    eprj0i1(k,k)=1;
    eprj1i1(k,k)=1;
    
end

elp0=log(cat(3,epj0i0s0,epj0i1s0,epj1i0s0,epj1i1s0));
elp1=log(cat(3,epj0i0s1,epj0i1s1,epj1i0s1,epj1i1s1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%A pos risk approach
mask_pos_p=elp1-elp0-(lpj1-lpj0)-permute(lpj1-lpj0,[2 1 3])>T;


%%%A neg risk approach 
mask_neg_p=elp0-elp1-(lpj0-lpj1)-permute(lpj0-lpj1,[2 1 3])>T;


for i=1:TF
    
    mask_pos_p(i,i,:)=0;
    mask_neg_p(i,i,:)=0;
    
end

%mask_pos_p(:,:,1)=0;
%mask_neg_p(:,:,1)=0;

mask_pos=mask_pos_p;
mask_neg=mask_neg_p;

