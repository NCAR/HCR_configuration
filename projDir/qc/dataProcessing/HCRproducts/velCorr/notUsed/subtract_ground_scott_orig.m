
[VG_mean, VG_med]=running_median_lag(VEL_ground,20);

VG_fir(1,1:length(VEL_ground))=Hubbert_filter(VEL_ground);
num_iterations1=10;
for i=2:num_iterations1
    VG_fir(i,1:length(VEL_ground))=Hubbert_filter(VG_fir(i-1,1:length(VEL_ground)));
end

std_VG_mean=std(VG_mean-VEL_ground);
i_outlier=find(abs(VG_mean-VEL_ground)>1.5*std_VG_mean);
VEL_ground_no_outliers=VEL_ground;
VEL_ground_no_outliers(i_outlier)=nan;

[VG_mean_no_outliers,VG_med_no_outliers]=running_median_lag(VEL_ground_no_outliers,50);

VEL_corr_ground_med2=VEL_corr;
VEL_corr_ground_med=VEL_corr;
VEL_corr_ground_mean2=VEL_corr;
VEL_corr_ground_mean=VEL_corr;
VEL_corr_ground_fir=VEL_corr;
% VEL_corr_ground_fir_1=VEL_corr;

VG_fir_tofilter(1:m)=nan;
for i=1:m
    VEL_corr_ground_med2(i,:)=VEL_corr(i,:)-VG_med_no_outliers(i);
    VEL_corr_ground_mean2(i,:)=VEL_corr(i,:)-VG_mean_no_outliers(i);
    
    VEL_corr_ground_med(i,:)=VEL_corr(i,:)-VG_med(i);
    VEL_corr_ground_mean(i,:)=VEL_corr(i,:)-VG_mean(i);
    
    if abs(VEL_ground(i)-VG_fir(1,i))<0.2
        VG_fir_tofilter(i)=VEL_ground(i);
    else
        VG_fir_tofilter(i)=VG_fir(1,i);
    end
end

VG_fir_final(1,1:length(VEL_ground))=Hubbert_filter(VG_fir_tofilter);
num_iterations2=3;
for i=2:num_iterations2
    VG_fir_final(i,1:length(VEL_ground))=Hubbert_filter(VG_fir_final(i-1,1:length(VEL_ground)));
end

VEL_ground_fir(1:m)=nan;
for i=1:m
    VEL_corr_ground_fir(i,:)=VEL_corr(i,:)-VG_fir_final(num_iterations2,i);
    VEL_ground_fir(i)=VEL_corr_ground_fir(i,ground_index(i));
    %VEL_corr_ground_fir_1(i,:)=VEL_corr(i,:)-VG_fir_final(1,i);
end

