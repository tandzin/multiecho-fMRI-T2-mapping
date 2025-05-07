wd=cd;
addpath(wd)
clear variables;close all
% warning off
warning off images:initSize:adjustingMag
warning off MATLAB:MKDIR:DirectoryExists
fullscreen = get(0,'ScreenSize');
screen_left=fullscreen(1);
screen_bottom=fullscreen(2);
screen_width=fullscreen(3);
screen_height=fullscreen(4);
taskbar_height=90;

format shortG
commandwindow
running_operation_string=['CBIA>>...please make sure spm12 toolbox is on your MATLAB path'];
disp(running_operation_string)
% spm_rmpath
% addpath I:\Code\MATLAB_ext_toolboxes\spm12
% SelectT2starDataset
warning off MATLAB:MKDIR:DirectoryExists
running_operation_string=['CBIA>>please enter full directory path containing TV-filtered echoes' newline...
    '      TV_Echo_1.nii, TV_Echo_2.nii and TV_Echo_3.nii' newline... 
    '      such as I:\\Data\\T2star\\4055B\\FCNI1_mag\\cwu\\TV_filtered_echoes\\-16   -6   10   10:'];
data_dir=input(running_operation_string,'s');
data_dir_default='I:\Data\T2star\4055B\FCNI1_mag\cwu\TV_filtered_echoes\-16   -6   10   10';

if isempty(data_dir)
    data_dir=data_dir_default;
end
running_operation_string=['CBIA>>...TV-filtered echoes will be read from' newline...
   '         ' data_dir];
disp(running_operation_string)
running_operation_string=['CBIA>>...reading 4-D TV-filtered echoes TV_Echo_1, TV_Echo_2, TV_Echo_3 '];
disp(running_operation_string)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T2star_get_3_series_of_echoes_script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
running_operation_string=['CBIA>>...concatenating 4-D TV-filtered echoes TV_Echo_1, TV_Echo_2, TV_Echo_3 ' newline...
    '         to 5-D datacube Echoes with echo number as the 5-th dimension'];
disp(running_operation_string)
Echoes=zeros([dims_echocube 3],'single');
Echoes(:,:,:,:,1)=Echo_1;
Echoes(:,:,:,:,2)=Echo_2;
Echoes(:,:,:,:,3)=Echo_3;
clear Echo_1 Echo_2 Echo_3

wd=cd;
cd ..% raw data is in a directory 2 levels higher
cd ..
data_dir=cd;
running_operation_string=['CBIA>>...unfiltered echoes will be read from' newline...
   '         ' data_dir];
disp(running_operation_string)
running_operation_string=['CBIA>>...reading 4-D unfiltered echoes Echo_1, Echo_2, Echo_3 ' newline...
    '         and Echo_1_dicom_header, Echo_2_dicom_header, Echo_3_dicom_header '];
disp(running_operation_string)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T2star_get_3_series_of_echoes_script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(wd);
running_operation_string=['CBIA>>...concatenating 4-D unfiltered echoes Echo_1, Echo_2, Echo_3 ' newline...
    '         to 5-D datacube Ref_echoes with echo number as the 5-th dimension'];
disp(running_operation_string)
Ref_echoes=zeros([dims_echocube 3],'single');
Ref_echoes(:,:,:,:,1)=Echo_1;
Ref_echoes(:,:,:,:,2)=Echo_2;
Ref_echoes(:,:,:,:,3)=Echo_3;
clear Echo_1 Echo_2 Echo_3
running_operation_string=['CBIA>>...reading TE(1),TE(2),TE(3) and TR from dicom headers '];
disp(running_operation_string)
TE(1)  = Echo_1_dicom_header.EchoTime*10^-3;
TE(2)  = Echo_2_dicom_header.EchoTime*10^-3;
TE(3)  = Echo_3_dicom_header.EchoTime*10^-3;
TR = Echo_1_dicom_header.RepetitionTime*10^-3;

subdir=mfilename;
mkdir(subdir);cd(subdir);
[H W D nr_volumes nr_echoes]=size(Echoes);
slc_max = D;
S=size(Echoes,3);
slc_default=int16(S/2);
slc=[];

disp(['CBIA>>maximum slice:' num2str(slc_max)])
commandwindow
slc_string=['CBIA>>pick a z-slice for 2D visualizations, or press RETURN for default: [' num2str(slc_default) ']:'];
slc=str2num(input(slc_string,'s'));
if isempty(slc)
    slc=slc_default;
end
if slc<1 || slc>S
    disp('Error:slice is out of range')
    return
end
slice_name = ['slc' num2str(slc)];
disp(['slc=' num2str(slc)])

%for statistics take all slices from all echoes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T2star_select_threshold_script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% saveastiff options
clear options;
% disp('Lossless LZW compression');
options.compress = 'lzw';
% disp('Overwrite to an existing file');
options.overwrite = true;
% disp('Disable message printing.');
options.message = false;
% disp('16 bit, grayscale image');

running_operation_string=['CBIA>>...interleaving 4-D TV-filtered echoes TV_Echo_1, TV_Echo_2, TV_Echo_3 ' newline...
    '         according to their absolute acquisition time ' newline...
    '         to 4-D datacube interleaved_echoes'];
disp(running_operation_string)

% interleave echoes 1,2,3 of all volumes
interleaved_echoes = zeros(H,W,D,nr_echoes*nr_volumes,'single');
% create time base for volumes
t_rf = zeros(nr_volumes,1);
% create time base for echoes
t_echo = zeros(nr_echoes*nr_volumes,1);
for vol_idx=1:nr_volumes
    t_rf(vol_idx) = TR*(vol_idx-1);
    for echo=1:3
        echo_idx=3*(vol_idx-1)+echo;
        t_echo(echo_idx)=t_rf(vol_idx)+TE(echo);
        interleaved_echoes(:,:,:,echo_idx)=Echoes(:,:,:,vol_idx,echo);
    end
end

nii_structure=make_nii(interleaved_echoes);
filename=['interleaved_echoes.nii'];
save_nii(nii_structure,filename);

% do_limitation=str2num(input('do_limitation:[press RETURN/0] ','s'));
do_limitation=[];
do_limitation=0;
if isempty(do_limitation)
    do_limitation=1;
end
min_T2=TE(1);% first echo time
max_T2=TR;% repetition time


first_echo_time=TE(1);
last_echo_time=TE(3);
display_range=[0 last_echo_time];

running_operation_string=['CBIA>>...calculating T2star values for all volumes of input data '];
disp(running_operation_string)
t=TE';
T_array=repmat(t,[1 H W]);
T_array=permute(T_array,[2 3 1]);
T2star_GS_slc=zeros(H,W,nr_volumes,'single');
T2star_GS=zeros(H,W,D,nr_volumes,'single');
Noise_mask_slc=zeros(H,W,nr_volumes,'single');% Noise_mask_slc == 1 means acknowledged pixels
Noise_mask=zeros(H,W,D,nr_volumes,'single');% Noise_mask == 1 means acknowledged pixels
Echo_array=zeros(H,W,nr_echoes);

tic
for slc_idx = 1:D
    interleaved_echoes_slc=squeeze(interleaved_echoes(:,:,slc_idx,:));
    for vol_idx = 1:nr_volumes%
        layer=1;
        for echo_time=t'
            Echo_array(:,:,layer)  = interleaved_echoes_slc(:,:,layer+3*(vol_idx-1));
            layer=layer+1;
        end
        Noise_mask_slc(:,:,vol_idx)=double...
            ((Echo_array(:,:,1)>Echo_array(:,:,2))&(Echo_array(:,:,2)>Echo_array(:,:,3))&(Echo_array(:,:,3)>0));
        
        Nonzero_mask_array=(Echo_array>0);
        nr_of_iterations=4;
        
        A_i=zeros(H,W,nr_of_iterations);
        B_i=zeros(H,W,nr_of_iterations);
        C_i=zeros(H,W,nr_of_iterations);
        G_i=zeros(H,W,nr_of_iterations);
        
        % golden section initializer
        golden_section=(sqrt(5)-1)/2;
        C_min=zeros(H,W);
        Min_Echo=min(Echo_array,[],3);
        Max_Echo=max(Echo_array,[],3);
        
        layer=3;
        C_max=Echo_array(:,:,layer);
        C_max=C_min;
        C=C_min;
        for iteration=1:nr_of_iterations
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % update b in E=a*exp(-bt)+c
            % weighted linear regression for B and B0 in exponential fit
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % C=zeros(H,W);%!!!!
            C_array=repmat(C,[1 1 size(t)]);
            D_array=Echo_array-C_array;% offset-corrected echo
            F_array=log(double(max(D_array,realmin('double'))));% prevent log from being -Inf, since 0*Inf=NaN
            
            D_array=max(D_array,0);% for weighting, negative offset-corrected echos are replaced with zeros
            Nonzero_mask_array=(D_array>0);
            N=sum(Nonzero_mask_array,3);
            TT_D_array=T_array.*D_array;
            FT_D_array=F_array.*D_array;
            
            D_square=sum(D_array.*D_array,3);
            TT_DD_1=sum(TT_D_array.*D_array,3);
            FT_DD_1=sum(FT_D_array.*D_array,3);
            TT_DD_F=sum(TT_D_array.*FT_D_array,3);
            TT_DD_T=sum(TT_D_array.*TT_D_array,3);
            NUM=TT_DD_1.*FT_DD_1-D_square.*TT_DD_F;
            DEN=D_square.*TT_DD_T-TT_DD_1.*TT_DD_1;
            % DEN_plus=DEN+double(DEN<=0).*realmin('single');%replace denominator with -> 0 wherever denominator==0;
            DEN_plus=DEN.*double(DEN>0)+double(DEN<=0).*realmin('single');%replace denominator with -> 0 wherever denominator==0;
            B=NUM./DEN_plus;%prevents B from being NaN
            B=B.*(DEN~=0)+realmax('single').*(DEN==0);% make T2_GS=0 whenever B infinite
            B=max(B,0);%only non-negative decay rates admissible by the physics
            % most prostate pixels within the body have value <320, so this limitation does not make sense for medical data
            % it only makes sense for artificial randomly generated exponentials
            if do_limitation
                B=min(B,1/min_T2);%limit T2_GS maximum ground truth
                B=max(B,1/max_T2);%limit T2_GS maximum ground truth
            end
            DEN_B0=D_square;
            DEN_B0_plus=max(DEN_B0,realmin('single'));% realmin will prevent division by zero
            % NUM_B0=(FT_DD_1+B.*(DEN>0).*TT_DD_1)+double(DEN<=0);
            NUM_B0=(FT_DD_1+B.*TT_DD_1);
            B0=NUM_B0./DEN_B0_plus;
            B0=min(B0,realmax('single'));%replace Inf by realmax('single') to prevent 0*Inf=NaN
            B0=max(B0,-realmax('single'));
            B0=B0.*(DEN~=0)-realmax('single').*(DEN==0);% DEN==0 pushes A-> zero
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % update a in y=a*exp(-bt)+c
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            A=exp(B0);
            A=min(A,realmax('single'));
            A=A.*(DEN~=0)+0.*(DEN==0);
            
            A_array=repmat(A,[1 1 size(t)]);
            B_array=repmat(B,[1 1 size(t)]);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % update c in y=a*exp(-bt)+c using golden section
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Fit_array=A_array.*exp(-B_array.*T_array)+C_array;
            R_array=Echo_array-Fit_array;% error between the echos and the current fit
            G=sum(R_array.*R_array,3);% cost function values for all pixels
            G=min(G,realmax('single'));
            if(iteration==1)%initialize inner GS point
                G_inner=G;
                x_left=C_min;
                x_right=C_max;
                x_inner=x_right-(x_right-x_left)*golden_section;%will lie to the left of x_probe
                A_i(:,:,iteration)=A;
                B_i(:,:,iteration)=B;
                C_i(:,:,iteration)=C;
                C=x_inner;
            elseif(iteration==2)% initialize probe GS point
                G_inner=G;
                x_probe=x_left+(x_right-x_left)*golden_section;%will lie to the right of x_inner
                A_i(:,:,iteration)=A;
                B_i(:,:,iteration)=B;
                C_i(:,:,iteration)=C;
                C=x_probe;
                G_max=G;
            else % next GS iteration
                G_probe=G;
                G_probe_LT_G_inner=G_probe<G_inner;
                G_probe_GE_G_inner=~G_probe_LT_G_inner;
                A_i(:,:,iteration)=G_probe_LT_G_inner.*A+G_probe_GE_G_inner.* A_i(:,:,iteration-1);
                B_i(:,:,iteration)=G_probe_LT_G_inner.*B+G_probe_GE_G_inner.* B_i(:,:,iteration-1);
                C_i(:,:,iteration)=G_probe_LT_G_inner.*C+G_probe_GE_G_inner.* C_i(:,:,iteration-1);
                
                % inner value for next iteration
                % x_inner is always the minimum
                x_inner_new=G_probe_LT_G_inner.*x_probe+G_probe_GE_G_inner.*x_inner;
                G_inner=G_probe_LT_G_inner.*G_probe+G_probe_GE_G_inner.*G_inner;
                % probe value for next iteration
                % if G_probe<G_inner, extrapolate x_probe beyond x_probe
                % if G_probe>=G_inner, extrapolate x_probe beyond x_inner
                x_probe=G_probe_LT_G_inner.*(x_inner+(x_probe-x_inner)/golden_section)...
                    +G_probe_GE_G_inner.*(x_probe+(x_inner-x_probe)/golden_section);
                % x_probe=min(x_probe,C_max);
                % x_probe=max(x_probe,C_min);
                x_inner=x_inner_new;
                C=x_probe;
            end
            G_i(:,:,iteration)=G_inner;
        end %for iteration=1:nr_of_iterations
        % ' nr_of_iterations='  num2str(nr_of_iterations)]);
        
        % G_GS=G_inner;
        G_last_LT_G_first=G_i(:,:,nr_of_iterations)<G_i(:,:,1);
        G_last_GE_G_first=~G_last_LT_G_first;
        A_GS=G_last_LT_G_first.*A+G_last_GE_G_first.* A_i(:,:,1);
        B_GS=G_last_LT_G_first.*B+G_last_GE_G_first.* B_i(:,:,1);
        C_GS=G_last_LT_G_first.*C+G_last_GE_G_first.* C_i(:,:,1);
        G_GS=G_last_LT_G_first.*G_i(:,:,nr_of_iterations)+G_last_GE_G_first.*G_i(:,:,1);
        G_last_GT_G_first=G_GS>G_i(:,:,1);
        % plot T2_GS for all pixels
        B_GS=max(B_GS,realmin('single'));
        T2_GS=1./B_GS;
        T2_GS=min(T2_GS,realmax('single'));
        % T2_GS=T2_GS.*(NUM~=0)+realmax('single').*(NUM==0);
        % T2_GS=T2_GS.*(DEN~=0)+0.*(DEN==0);
        if do_limitation
            T2_GS=max(T2_GS,min_T2);
            T2_GS=min(T2_GS,max_T2);
        end
        T2star_GS_slc(:,:,vol_idx)= T2_GS;
    end % for vol_idx = 1:nr_volumes%
    T2star_GS(:,:,slc_idx,:)= T2star_GS_slc;
    Noise_mask(:,:,slc_idx,:)= Noise_mask_slc;    
end % for slc_idx = 1:D
ElapsedTime=toc;
disp(['ElapsedTime=' num2str(ElapsedTime)])
Noise_mask(isnan(T2star_GS))=0;%replace NaN with 0
% Noise_mask(isinf(T2star_GS))=0;%replace inf with 0
%replace low magnitude noise with 0

running_operation_string=['CBIA>>...saving noise mask and unmasked/masked T2star values as 4-D .nii files to ' newline...
    '         ' cd];
disp(running_operation_string)

Noise_mask(squeeze(Echoes(:,:,:,:,1))<threshold)=0;%use only echo 1; echoes 2 and 3 are darker
nii_structure=make_nii(Noise_mask);
filename=['Noise_mask.nii'];
save_nii(nii_structure,filename);

nii_structure=make_nii(T2star_GS);
filename=['T2star_GS.nii'];
save_nii(nii_structure,filename);

T2star_GS_masked=T2star_GS.*Noise_mask;
nii_structure=make_nii(T2star_GS_masked);
filename=['T2star_GS_masked.nii'];
save_nii(nii_structure,filename);

running_operation_string=['CBIA>>...extracting 3-D data at your selected slice number '...
    num2str(slc) ' for visualization'];
disp(running_operation_string)

subdir=['slc' num2str(slc)];
[mkdir_status,mkdir_message] = mkdir(subdir);
cd(subdir);
running_operation_string=...
    ['CBIA>>...saving raw/limited/masked 3-D T2stars and noise mask at your selected slice number in'...
    newline cd ];
disp(running_operation_string)

T2star_GS_slc=squeeze(T2star_GS(:,:,slc,:));
T2star_GS_slc_name=[ 'T2star_GS_slc_' num2str(slc) '.tif'];
saveastiff(T2star_GS_slc, T2star_GS_slc_name, options);

Noise_mask_slc=squeeze(Noise_mask(:,:,slc,:));
Noise_mask_slc_name=[ 'Noise_mask_slc_' num2str(slc) '.tif'];
saveastiff(Noise_mask_slc, Noise_mask_slc_name, options);

T2star_GS_slc_masked=T2star_GS_slc.*Noise_mask_slc;
T2star_GS_masked_slc_name=[ 'T2star_GS_masked_slc_' num2str(slc) '.tif'];
saveastiff(T2star_GS_slc_masked, T2star_GS_masked_slc_name, options);
masked_pixel_ratio=1-sum(Noise_mask_slc(:))/numel(Noise_mask_slc(:));
disp(['masked_pixel_ratio = ' num2str(masked_pixel_ratio)])

T2star_GS_slc_95prc=T2star_GS_slc_masked;
[F,X] = ecdf(T2star_GS_slc_95prc(:));
tail=X(F>0.9);
T2star_GS_slc_95prc(T2star_GS_slc_95prc>=tail(1))=tail(1);%replace 95% percentile with 0
T2star_GS_slc_95prc_name=[ 'T2star_GS_limited_slc_' num2str(slc) '.tif'];
saveastiff(T2star_GS_slc_95prc, T2star_GS_slc_95prc_name, options);

running_operation_string=['CBIA>>...To visually inspect echo and T2star time courses at voxels of your choice,' newline...
    '         pick repeatedly those you are interested in.'];
disp(running_operation_string)
pause(5)
% show face image for pixel selection
T2star_title = 'T2star_face(:,:,slc,1)';
T2star_face_fig=100;
figure(T2star_face_fig);set(T2star_face_fig,'Name',T2star_title,...
    'Position',[.5*screen_width+2 105 .5*screen_width-5 screen_height-189])%[left,bottom,width,height];
imshow(T2star_GS_slc_95prc(:,:,1),[],'InitialMagnification','fit','Border','tight')

set(gcf,'color','blue')
%     axis square
MarkerFaceColor='w';
%     MarkerEdgeColor = input('MarkerEdgeColor [r g b c m y k w]:','s');
%     marker_shape = input('marker shape [+ o * . x s d ^ v > < p h ]:','s');
MarkerEdgeColor = 'g';
marker_shape = 'o';
line_descriptor=['-' MarkerEdgeColor marker_shape];
MarkerSize=10;
pick_another_pixel=true;
while pick_another_pixel
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % pick a pixel for the echoes plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure(T2star_face_fig)
    uiwait(msgbox('Click a point'));
    [col,row] = ginputc(1,'Color','g');% select column and row of the image in focus
    %     rng('shuffle');
    %     col=floor(W*rand+1);
    %     row=floor(H*rand+1);
    hold on
    plot(col,row,line_descriptor,...
        'LineWidth',2,...
        'MarkerFaceColor',MarkerFaceColor,...
        'MarkerEdgeColor',MarkerEdgeColor,...
        'MarkerSize',MarkerSize)
    hold off
    row=min(int16(row),int16(W));
    col=min(int16(col),int16(H));
    echo_pencils_fig=200;
    figure(echo_pencils_fig);set(echo_pencils_fig,'Name',['echo_pencils at ' num2str(row) ',' num2str(col)],...
        'Position',[5 105 .5*screen_width-5 screen_height-190])%[left,bottom,width,height];
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
    saveas(echo_pencils_fig,['echo_pencils at ' num2str(row) ',' num2str(col) '.png']);
    
    T2star_pencil_fig=300;
    figure(T2star_pencil_fig);set(T2star_pencil_fig,'Name',['T2star_pencil at ' num2str(row) ',' num2str(col)],...
        'Position',[5 105 .5*screen_width-5 screen_height-190])%[left,bottom,width,height];
    y_pencil=squeeze(T2star_GS_slc_95prc(row,col,:));
    
    y_ax_length=1.1*max(y_pencil);
    if y_ax_length==0
        y_ax_length=TE(3);
    end
    axis([0 max(t_rf) 0 y_ax_length])% axis length must be >0
    box on
    hold on
    plot(t_rf, y_pencil, '.r')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    T2_star_summary_plot_script
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pause('on')
    wait_time = 5;
    pause(wait_time)    
%     pause
    commandwindow
    pick_another_pixel=str2num(input('         pick_another_pixel:[press RETURN for YES/0 for NO] ','s'));
    %     pick_another_pixel=false;
    if isempty(pick_another_pixel)
        pick_another_pixel=1;
        close (T2star_pencil_fig);
        if ishandle(echo_pencils_fig)
            close(echo_pencils_fig);
        end
        if ishandle(Echoes_T2stars_fig)
            close(Echoes_T2stars_fig);
        end
    end
end
commandwindow
