wd=cd;
addpath(wd)
clear variables;close all hidden
warning off
% load('GS_statistics.mat','GS_statistics')
format shortG
commandwindow
running_operation_string=['CBIA>>please make sure spm12 toolbox is on your MATLAB path'];
disp(running_operation_string)
% spm_rmpath
% addpath I:\Code\MATLAB_ext_toolboxes\spm12
fullscreen = get(0,'ScreenSize');
screen_left=fullscreen(1);
screen_bottom=fullscreen(2);
screen_width=fullscreen(3);
screen_height=fullscreen(4);
% SelectT2starDataset
warning off MATLAB:MKDIR:DirectoryExists
data_dir_string=['CBIA>>please enter full directory path containing' newline...
    '      xxx_Echo_N.nii and' newline... 
    '      yyy_Echo_N_dicom_header.mat,' newline...
    '      such as I:\\Data\\T2star\\4055B\\FCNI1_mag\\cwu:'];
data_dir=input(data_dir_string,'s');
data_dir_default='I:\Data\T2star\4055B\FCNI1_mag\cwu';
if isempty(data_dir)
    data_dir=data_dir_default;
end

running_operation_string=['CBIA>>...unfiltered echoes will be read from' newline...
   '         ' data_dir];
disp(running_operation_string)
cd(data_dir)
contents=dir;
running_operation_string=['CBIA>>...reading 4-D echoes Echo_1, Echo_2, Echo_3 ' newline...
    '         and Echo_1_dicom_header, Echo_2_dicom_header, Echo_3_dicom_header '];
disp(running_operation_string)
for fileindex=1:numel(contents)
    filename=contents(fileindex).name;
    namelength=numel(filename);
    k = strfind(filename,'Echo_');
    if ~isempty(k)
        input_variable_name = filename(k:(namelength-4));
        input_variable = genvarname(input_variable_name);
        if (strcmp(filename((namelength-2):namelength),'mat'))
            load(filename);
            input_data = dicom_header;
        elseif (strcmp(filename((namelength-2):namelength),'nii'))
            V=spm_vol(filename);
            input_data = single(spm_read_vols(V));
        else
            continue
        end
        eval([input_variable '= input_data;'])
    end
end
running_operation_string=['CBIA>>...reading TE(1),TE(2),TE(3) and TR from dicom headers '];
disp(running_operation_string)
TE(1)  = Echo_1_dicom_header.EchoTime*10^-3;
TE(2)  = Echo_2_dicom_header.EchoTime*10^-3;
TE(3)  = Echo_3_dicom_header.EchoTime*10^-3;
TR = Echo_1_dicom_header.RepetitionTime*10^-3;
if isequal(size(Echo_1),size(Echo_2)) && isequal(size(Echo_2),size(Echo_3))
    dims_echocube = size(Echo_1);
else
    disp('Error:numbers of echoes do not match')
    return
end

running_operation_string=['CBIA>>...concatenating 4-D echoes Echo_1, Echo_2, Echo_3 ' newline...
    '         to 5-D datacube Echoes with echo number as the 5-th dimension'];
disp(running_operation_string)
Echoes=zeros([dims_echocube 3],'single');
Echoes(:,:,:,:,1)=Echo_1;
Echoes(:,:,:,:,2)=Echo_2;
Echoes(:,:,:,:,3)=Echo_3;
clear Echo_1 Echo_2 Echo_3

subdir='TV_filtered_echoes';
% disp(subdir);
mkdir(subdir);cd(subdir);
close all
% scriptname = mfilename;

S=size(Echoes,3);
slc_default=int16(S/2);
slc=[];
commandwindow
slc_string=['CBIA>>pick a z-slice, or press RETURN for default: [' num2str(slc_default) ']:'];
slc=str2num(input(slc_string,'s'));
if isempty(slc)
    slc=slc_default;
end
if slc<1 || slc>S
    disp('Error:z-slice is out of range')
    return
end

nr_threads=[];
nr_threads=str2num(input('CBIA>>enter number of CPU threads, or press RETURN for 8 threads:','s'));
if isempty(nr_threads)
    nr_threads=8;
end
mu_default=2^-10;
% mu_default=1;
mu_string=['CBIA>>enter log2(mu), or press RETURN for default [' num2str(log2(mu_default)) ']:'];
mu=2^str2num(input(mu_string,'s'));
if isempty(mu)
    mu=mu_default;
end %if isempty(mu)
log2_mu=log2(mu);

beta_default=2^-4;
% beta_default=1;
beta_string=['CBIA>>enter log2(beta), or press RETURN for default [' num2str(log2(beta_default)) ']:'];
beta=2^str2num(input(beta_string,'s'));
if isempty(beta)
    beta=beta_default;
end %if isempty(beta)
log2_beta=log2(beta);

gamma=[];
% gamma=str2num(input('CBIA>>enter gamma, or press RETURN for default:','s'));
if isempty(gamma)
    gamma=1;
end

delta=[];
% delta=str2num(input('CBIA>>enter delta, or press RETURN for default:','s'));
if isempty(delta)
    delta=1;
end

% nr_shrinks_default=100;
% nr_shrinks_default=50;
nr_shrinks_default=10;
nr_shrinks_string=['CBIA>>enter nr_shrinks, or press RETURN for default [' num2str(nr_shrinks_default) ']:'];
nr_shrinks=str2num(input(nr_shrinks_string,'s'));
if isempty(nr_shrinks)
    nr_shrinks=nr_shrinks_default;
end

nr_grads_default=5;
nr_grads_string=['CBIA>>enter nr_grads, or press RETURN for default [' num2str(nr_grads_default) ']:'];
nr_grads=str2num(input(nr_grads_string,'s'));
if isempty(nr_grads)
    nr_grads=nr_grads_default;
end

subdir=[num2str(log2_mu) '   ' num2str(log2_beta) '   ' num2str(nr_shrinks) '   ' num2str(nr_grads)];
[mkdir_status,mkdir_message] = mkdir(subdir);
cd(subdir);
running_operation_string=['CBIA>>...TV filtering results will be stored in' newline...
   '         ' cd];
disp(running_operation_string)

% create Gaussian filter
hsize=2;
sigma=hsize/2;
gauss_filter_2D = fspecial('gaussian', hsize, sigma);% creates a two-dimensional filter h of the specified type.
TV_echoes=zeros(size(Echoes));
% plot_single_pixel_echo_series(Echoes,TV_echoes,slc) %for testing only
set(0, 'DefaultLineLineWidth', 2);
for echo_nr=1:3
    commandwindow
    running_operation_string=['CBIA>>...temporal TV filtering of Echo_' num2str(echo_nr)];
    disp(running_operation_string)
    Img=squeeze(Echoes(:,:,:,:,echo_nr));  
    b=double(Img);
    u_0=double(Img);
    tic
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [u_TV, out]=ROF_ADMM_1D_alg(b,u_0,mu,beta,gamma,delta,nr_shrinks,nr_grads,nr_threads);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ElapsedTime = toc;
    commandwindow
    disp(['ElapsedTime= ' num2str(ElapsedTime)]);
    %save the original 4D array as 4d nifti
    TV_echoes(:,:,:,:,echo_nr)=u_TV;
    nii_structure=make_nii(single(b));
    filename=['Echo' num2str(echo_nr) '.nii'];
    save_nii(nii_structure,filename);
    %save the denoised 4D array as 4d nifti
    nii_structure=make_nii(single(u_TV));
    filename=['TV_Echo_' num2str(echo_nr) '.nii'];
    save_nii(nii_structure,filename);
    
    figure(echo_nr);clf(echo_nr);
    taskbar_height=90;
    set(echo_nr,'Name','comparison',...
    'Outerposition',[screen_left screen_bottom+taskbar_height screen_width screen_height-taskbar_height]);

    subplot(231); imshow(Img(:,:,slc,1),[]);
    title('original image');
    
    subplot(232); imshow(u_TV(:,:,slc,1),[]);
    title('TV denoised image');    
    
    %Gaussian filtering
    u_gaussed = imfilter(double(Img),gauss_filter_2D,'replicate');%boundary padding by border replication
    subplot(233); imshow(u_gaussed(:,:,slc,1),[]);
    title('Gaussian\_filtered\_image');
    
    iter=1:numel(out.w_norm);
    subplot(234); plot(iter,out.w_norm,'g',iter,out.TV_norm,'r')
%     subplot(234); plot(iter,out.w_norm,'g',iter,out.TV_norm,'r','LineWidth',1)
    title('w\_norm (green) TV (red)');
    
    subplot(235); plot(iter,out.data_term,'g')
%     subplot(235); plot(iter,out.data_term,'g','LineWidth',2)
    title('data\_term (green)');
    
    subplot(236); plot(iter,out.P_k,'g',iter,out.Lambda_k,'r')
%     subplot(236); plot(iter,out.P_k,'g',iter,out.Lambda_k,'r','LineWidth',3)
    title('primal (green) dual (red)');
    saveas(gcf,['Echo' num2str(echo_nr) '.png']);
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...% needed for pdfs
    'PaperSize',[screenposition(3:4)]);% needed for pdfs
    saveas(gcf,['Echo' num2str(echo_nr) '.pdf'],'pdf');
    pause(5)
end
commandwindow
running_operation_string=['CBIA>>...saving TV filtering results to' newline...
   '         ' cd];
disp(running_operation_string)

running_operation_string=['CBIA>>...temporal TV filtering of the 5D datacube has been completed.' newline...
    '         To visually compare unfiltered and TV-filtered echoes at voxels of your choice,' newline...
    '         pick repeatedly those you are interested in.'];
disp(running_operation_string)
pause(5)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_single_pixel_echo_series(Echoes,TV_echoes,slc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear Echoes TV_echoes
%cd('I:\Code\matlab\jm\m_files');
% whistle
