function [x0dd_pos x0dd_neg x1dd_pos x1dd_neg]=dob_data_gen_1(x0,x1,mask_pos,mask_neg)

TF=size(x0,2);

x0c=1-x0;
x1c=1-x1;

x03=zeros(TF,TF,4);
%x13=zeros(TF,TF,4);

x0d_pos=zeros(size(x0,1),TF*4);
x0d_neg=zeros(size(x0,1),TF*4);


x1d_pos=zeros(size(x1,1),TF*4);
x1d_neg=zeros(size(x1,1),TF*4);


for i=1:size(x0,1)
    
    cx0=x0(i,:);
    cx0c=x0c(i,:);
    
    x03(:,:,1)=cx0c'*cx0c;
    x03(:,:,2)=cx0'*cx0c;
    x03(:,:,3)=cx0c'*cx0;
    x03(:,:,4)=cx0'*cx0;
    
    x03_pos=squeeze(sum(mask_pos.*x03,1))./squeeze(sum(mask_pos,1)+10^-10);
    x03_neg=squeeze(sum(mask_neg.*x03,1))./squeeze(sum(mask_neg,1)+10^-10);
    
    
    x0d_pos(i,:)=x03_pos(:)';
    x0d_neg(i,:)=x03_neg(:)';
    
    
    
end
   

%%%%%%%%%%%%%%%%%%%%now class 1 but name x0 stands :(
for i=1:size(x1,1)
    
    cx0=x1(i,:);
    cx0c=x1c(i,:);
    
    x03(:,:,1)=cx0c'*cx0c;
    x03(:,:,2)=cx0'*cx0c;
    x03(:,:,3)=cx0c'*cx0;
    x03(:,:,4)=cx0'*cx0;
    
    x03_pos=squeeze(sum(mask_pos.*x03,1))./squeeze(sum(mask_pos,1)+10^-10);
    x03_neg=squeeze(sum(mask_neg.*x03,1))./squeeze(sum(mask_neg,1)+10^-10);
    
    
    x1d_pos(i,:)=x03_pos(:)';
    x1d_neg(i,:)=x03_neg(:)';
    
    
    
    
    
end


   
x0dd_pos=x0d_pos;
x0dd_neg=x0d_neg;

x1dd_pos=x1d_pos;
x1dd_neg=x1d_neg;





    