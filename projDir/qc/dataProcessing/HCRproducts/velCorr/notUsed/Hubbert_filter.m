function Field_F=Hubbert_filter(Field)

n=length(Field);
%[m n]=size(Field);
% Apply Hubbert FIR filter to Field
% Coefficients
b=[1.625807356e-2 2.230852545e-2 2.896372364e-2 3.595993808e-2...
    4.298744446e-2 4.971005447e-2 5.578764970e-2 6.089991897e-2...
    6.476934523e-2 6.718151185e-2 6.800100000e-2 6.718151185e-2...
    6.476934523e-2 6.089991897e-2 5.578764970e-2 4.971005447e-2...
    4.298744446e-2 3.595993808e-2 2.896372364e-2 2.230852545e-2...
    1.625807356e-2];

m=floor(length(b)/2);
for i=m+1:n-m
    Field_F(i)=sum(Field(i-m:i+m).*b);
    % Account for filter power loss (1-sum(b)=0.042)
    Field_F(i)=Field_F(i)*1.042;
    
    %Field_F=filter2(b,Field);
end

for i=1:m
    Field_F(i)=mean(Field(1:i+m));
%     Field_F(i)=sum(Field(1:i+m).*b(length(b)-length(1:i+m)+1:length(b))');
%     Field_F(i)=Field_F(i).*(1+1-sum(b(length(b)-length(1:i+m)+1:length(b))));
end

for i=n-m+1:n
    Field_F(i)=mean(Field(i-m:n));
end




