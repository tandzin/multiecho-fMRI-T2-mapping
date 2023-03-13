cd(data_dir)
contents=dir;
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
if isequal(size(Echo_1),size(Echo_2)) && isequal(size(Echo_2),size(Echo_3))
    dims_echocube = size(Echo_1);
else
    disp('Error:numbers of echoes do not match')
    return
end
