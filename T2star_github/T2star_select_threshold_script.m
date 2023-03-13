Echo_1=squeeze(Echoes(:,:,:,:,1));
% Echo_2=squeeze(Echoes(:,:,:,:,2));
% Echo_3=squeeze(Echoes(:,:,:,:,3));
Echo_1=Echo_1(:);
Echo_1=Echo_1(Echo_1>0);

nr_bins=ceil(max(Echo_1))+1;
Echo_1_histogram=hist(Echo_1,nr_bins);
clear Echo_1
Echo_1_cdf=zeros(1,nr_bins);
for bin=2:nr_bins%delete black from the histogram
    Echo_1_cdf(bin)=Echo_1_cdf(bin-1)+Echo_1_histogram(bin);
end
Echo_1_cdf=Echo_1_cdf/Echo_1_cdf(end);
Echo_1_cdf_fig=1;

threshold_selection_fig=4;
figure_size=(screen_height-taskbar_height);
figure(threshold_selection_fig);set(threshold_selection_fig,'Name','Threshold selection',...
    'Position',[screen_left screen_bottom+taskbar_height figure_size figure_size]);

dark_values=find((Echo_1_cdf>0)&(Echo_1_cdf<0.7));
nr_dark_values=numel(dark_values);

subplot(2,2,1)
plot(Echo_1_histogram(2:nr_bins))
% plot(Echo_1_histogram(2:nr_dark_values))
title('Echo\_1\_histogram')

subplot(2,2,2)
plot(Echo_1_cdf)
title('Echo\_1\_cdf')

sliding_offset=ceil(nr_dark_values/10);
sliding_secants=Echo_1_cdf((sliding_offset+dark_values(1)):(sliding_offset+dark_values(end)))-Echo_1_cdf(dark_values(1):dark_values(end));
subplot(2,2,3)
plot(sliding_secants)
title('sliding\_secants')


auto_threshold=find(sliding_secants==min(sliding_secants(:)))+sliding_offset/2;
threshold=auto_threshold(1);
plot_mask=true;
while(plot_mask)
    plot_mask=false;
    %for the background mask 2D image take the selected slice from first echo of the
    %first volume
    vol=1;
    echo=1;
    background_mask=squeeze(Echoes(:,:,slc,vol,echo)>threshold);
    background_mask_fig=3;
    background_mask_name=['background mask(slc,1,1) for threshold=' num2str(threshold)];
    subplot(2,2,4)
    imshow(background_mask,[])
    background_mask_name=['background mask(slc,1,1) for threshold=' num2str(threshold)];
    title(background_mask_name)
    threshold_string=['CBIA>>press RETURN to accept the threshold, or enter a different one: [' num2str(threshold) ']:'];
    figure(threshold_selection_fig);
    manual_threshold=str2num(input(threshold_string,'s'));
    if ~isempty(manual_threshold)
        threshold=manual_threshold;
        plot_mask=true;
    end
end
