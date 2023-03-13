Echoes_T2stars_fig=400;
fh=figure(Echoes_T2stars_fig);set(Echoes_T2stars_fig,'Name',['Echoes_T2stars at ' num2str(row) ',' num2str(col)],...
    'Outerposition',[screen_left screen_bottom+taskbar_height screen_width screen_height-taskbar_height]);
subplot(2,3,1);
echo_face = squeeze(Echoes(:,:,slc,1,1));
imshow(echo_face,[])
title('filtered echo')
hold on
plot(col,row,line_descriptor,...
    'LineWidth',2,...
    'MarkerFaceColor',MarkerFaceColor,...
    'MarkerEdgeColor',MarkerEdgeColor,...
    'MarkerSize',MarkerSize)
hold off

subplot(2,3,2);
% imshow(echo_face,[],'InitialMagnification','fit','Border','tight')
imshow(Noise_mask_slc(:,:,1),[])
title('Noise mask')
hold on
plot(col,row,line_descriptor,...
    'LineWidth',2,...
    'MarkerFaceColor',MarkerFaceColor,...
    'MarkerEdgeColor',MarkerEdgeColor,...
    'MarkerSize',MarkerSize)
hold off

subplot(2,3,3);
imshow(T2star_GS_slc_95prc(:,:,1),[])
title('T2star estimate')
hold on
plot(col,row,line_descriptor,...
    'LineWidth',2,...
    'MarkerFaceColor',MarkerFaceColor,...
    'MarkerEdgeColor',MarkerEdgeColor,...
    'MarkerSize',MarkerSize)
hold off

ref_echo_face = squeeze(Ref_echoes(:,:,slc,1,1));
subplot(2,3,4);imshow(ref_echo_face,[]);
title('reference echo')
% imshow(echo_face,[],'InitialMagnification','fit','Border','tight')

hold on
plot(col,row,line_descriptor,...
    'LineWidth',2,...
    'MarkerFaceColor',MarkerFaceColor,...
    'MarkerEdgeColor',MarkerEdgeColor,...
    'MarkerSize',MarkerSize)
hold off

subplot(2,3,5);
original_4D_pencil=squeeze(Ref_echoes(row,col,slc,:,:));
processed_4D_pencil=squeeze(Echoes(row,col,slc,:,:));
original_4D_pencil=original_4D_pencil';%transpose 2D matrix to permute dimensions
processed_4D_pencil=processed_4D_pencil';%transpose 2D matrix to permute dimensions

y_max=1.1*max(max(original_4D_pencil(:)),max(processed_4D_pencil(:)));
axis([0 max(t_echo) 0 y_max])% axis length must be >0
box on
hold on
plot(t_echo,original_4D_pencil(:), '.r')%
plot(t_echo,processed_4D_pencil(:), '.k')%
hold off
title('original and filtered three echo pencils')

subplot(2,3,6);
T2star_pencil=squeeze(T2star_GS_slc(row,col,:));
% T2star_pencil=squeeze(T2star_GS_slc_95prc(row,col,:));
axis([0 max(t_rf) 0 1.1*max(T2star_pencil)])% axis length must be >0
box on
hold on
plot(t_rf, T2star_pencil, '.r')
title('estimated T2star pencil')
hold off
saveas(Echoes_T2stars_fig,['Echoes_T2stars at ' num2str(row) ',' num2str(col) '.png']);




