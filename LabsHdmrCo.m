function [score0 score1 fa_vec, miss_vec]=LabsHdmrCo(xtr0,xtr1,xts0,xts1,T1,T2,T3,tvec)

TF=size(xtr1,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%binarize SNPs
%%train data
x0b1=(xtr0>0);
x0b2=(xtr0>1);

x1b1=(xtr0>0);
x1b2=(xtr1>1);


x0train=[x0b1 x0b2];
x1train=[x1b1 x1b2];


%%test data
x0b1=(xts0>0);
x0b2=(xts0>1);

x1b1=(xts0>0);
x1b2=(xts1>1);


x0test=[x0b1 x0b2];
x1test=[x1b1 x1b2];

clear x0b1 x0b2 x1b1 x1b2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%some initial values

fav=zeros(length(GFvec),length(tvec));
missv=fav;


p0nb=zeros(3,TF);
p1nb=p0nb;

for i=0:2
    
    if i==0
        
        nc0=sum(x0train==0);
        nc1=sum(x1train==0);
        
    elseif  i==1
        
        nc0=sum(x0train==1);
        nc1=sum(x1train==1);
        
    else
        
        nc0=sum(x0train==1);
        nc1=sum(x1train==1);
        
    end   
    
    p0nb(i+1,:)=nc0;
    p1nb(i+1,:)=nc1;   
    
    
end

%%%%%%%%%%%%%%%
%%now that you hvae counted use it for the binary SNP feature selection too!

a0s=zeros(2,2*TF);
a1s=a0s;

a0s(1,1:TF)=p0nb(1,:);
a0s(2,1:TF)=p0nb(2,:);
a0s(1,1+TF:end)=n0train-p0nb(3,:);
a0s(2,1+TF:end)=p0nb(3,:);

a1s(1,1:TF)=p1nb(1,:);
a1s(2,1:TF)=p1nb(2,:);
a1s(1,1+TF:end)=n1train-p1nb(3,:);
a1s(2,1+TF:end)=p1nb(3,:);

ats=a0s+a1s;

a0s=a0s+c;
a1s=a1s+c;
ats=ats+c;


lpg2=sum(gammaln(a0s))-gammaln(sum(a0s))+sum(gammaln(a1s))-gammaln(sum(a1s))-  (  sum(gammaln(ats))-gammaln(sum(ats))  );


ipv0=a0s./repmat(sum(a0s),2,1);
ipv1=a1s./repmat(sum(a1s),2,1);



[spos sid]=sort(lpg2,'descend');



c=3;  %%%%prior for OBF
T=T3; %%%%%T3
TNB=T1;  %%%T1
TNB2=T2;  %%%%T2


cnt=0;

for GF=GFvec
    
    cnt=cnt+1;
        
    isi=sid(1:GF);
    isiNB=sid(1:GF);
    
    x0d=x0train(:,isi);
    x1d=x1train(:,isi);
    
    x0td=x0test(:,isi);
    x1td=x1test(:,isi);
    
    x0dNB=x0train(:,isiNB);
    x1dNB=x1train(:,isiNB);
    
    x0tdNB=x0test(:,isiNB);
    x1tdNB=x1test(:,isiNB);
    %%%%%%%%%%%%%%%%%%
    
    [mask_pos mask_neg]=mask_gen_2(double(x0d),double(x1d),c,T,ipv0(:,isi),ipv1(:,isi));
    
    
    [x0dd_pos x0dd_neg x1dd_pos x1dd_neg]=dob_data_gen_1(double(x0d),double(x1d),mask_pos,mask_neg);
    
    [x0tdd_pos x0tdd_neg x1tdd_pos x1tdd_neg]=dob_data_gen_1(double(x0td),double(x1td),mask_pos,mask_neg);
    
    x0d2=[x0dNB x0dd_pos x0dd_neg];
    x1d2=[x1dNB x1dd_pos x1dd_neg];
    
    
    x0td2=[x0tdNB x0tdd_pos x0tdd_neg];
    x1td2=[x1tdNB x1tdd_pos x1tdd_neg];
    
    
    %%%%%%%%%%%%%%%%%%%%
    %%train the classifier
    diff_mean=log(mean(x1d2)+10^-10)-log(mean(x0d2)+10^-10);
    
    w=diff_mean.*([abs(diff_mean(1:GF))>TNB  abs(diff_mean(1+GF:end))>TNB2]);
    
    w=w./sqrt(sum(w.^2));
    
    s0=sum(repmat(w,size(x0td2,1),1).*x0td2,2);
    
    s1=sum(repmat(w,size(x1td2,1),1).*x1td2,2);
    
    
    %%%%%%compute false alarm and miss probs
    for i=1:length(tv)
        
        fav(cnt,i)=sum(s0>tv(i))/length(s0);
        missv(cnt,i)=sum(s1<=tv(i))/length(s1);
        
    end
    
    
end

%%%%%%%%%%%%%%%
%%output results

score0=s0;
score1=s1;
fa_vec=fav;
miss_vec=missv;


