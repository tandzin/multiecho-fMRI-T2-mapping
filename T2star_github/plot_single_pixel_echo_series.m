function plot_single_pixel_time_series(original_5D_image, processed_5D_image,slc)

[R C S V E]=size(original_5D_image);
fullscreen = get(0,'ScreenSize');
screen_left=fullscreen(1);
screen_bottom=fullscreen(2);
screen_width=fullscreen(3);
screen_height=fullscreen(4);

% show face image for pixel selection
echo_face = squeeze(original_5D_image(:,:,slc,1));
echo_title = 'echo_face(:,:,slc,1)';
echo_face_fig=100;
taskbar_height=90;
figure(echo_face_fig);
set(echo_face_fig,'Name',echo_title,...
    'Outerposition',[screen_left screen_bottom+taskbar_height screen_width/2 screen_height-taskbar_height]);

imshow(echo_face,[],'InitialMagnification','fit','Border','tight')
set(gcf,'color','blue')
MarkerFaceColor='w';
MarkerEdgeColor = 'g';
marker_shape = 'o';
line_descriptor=['-' MarkerEdgeColor marker_shape];
MarkerSize=10;
pick_another_pixel=true;
while pick_another_pixel
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % pick a pixel for the echoes plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%0
    figure(echo_face_fig)
    uiwait(msgbox('Click a point'));
    [col,row] = ginputc(1,'Color','g');% select column and row of the image in focus
    hold on
    plot(col,row,line_descriptor,...
        'LineWidth',2,...
        'MarkerFaceColor',MarkerFaceColor,...
        'MarkerEdgeColor',MarkerEdgeColor,...
        'MarkerSize',MarkerSize)
    hold off
    row=min(int16(row),int16(R));
    col=min(int16(col),int16(C));
    echo_pencils_fig=200;
    figure(echo_pencils_fig);
    set(echo_pencils_fig,'Name',['echo_pencils at ' num2str(row) ',' num2str(col)],...
        'Outerposition',[screen_left+screen_width/2 screen_bottom+taskbar_height screen_width/2 screen_height-taskbar_height]);
    original_4D_pencil=squeeze(original_5D_image(row,col,slc,:,:));
    processed_4D_pencil=squeeze(processed_5D_image(row,col,slc,:,:));
    original_4D_pencil=original_4D_pencil';%transpose 2D matrix to permute dimensions
    processed_4D_pencil=processed_4D_pencil';%transpose 2D matrix to permute dimensions

    y_max=1.1*max(max(original_4D_pencil(:)),max(processed_4D_pencil(:)));
    y_max=max(y_max,.1);% axis length must be >0
    axis([0 numel(original_4D_pencil) 0 y_max])% axis length must be >0
    box on
    hold on
    plot(original_4D_pencil(:), '.r')
    plot(processed_4D_pencil(:), '.k')
    saveas(echo_pencils_fig,['echo_pencils at ' num2str(row) ',' num2str(col) '.png']);

    pause('on')
    wait_time = 5;
    pause(wait_time)
    commandwindow
    pick_another_pixel=str2num(input('pick_another_pixel:[return/0] ','s'));
    if isempty(pick_another_pixel)
        pick_another_pixel=1;
        if ishandle(echo_pencils_fig)
            close(echo_pencils_fig);
        end
    end
end
saveas(echo_face_fig,['echo_face ' '.png']);

end