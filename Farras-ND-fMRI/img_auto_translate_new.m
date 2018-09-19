function img_auto_translate_new(coT1,prefix,RunAtSize,emlopt)

%The array coT1 points to the images that require shifting dcm2nii converted images to match T1
%template

%coT1=Dir2Arr([start_dir filesep 'test' filesep 'ANALYZE_translate_test'],'co*.img');

if nargin<3
    RunAtSize=0;
    emlopt=0;
    if nargin<2
        prefix='backup_';
    end
end


%Set RunAtlasing command

RunAtCmd='RunAtlasing_generic(coT1{t},RunAtSize,1,0)';
%RunAtCmd='RunAtlasing_generic(coT1,RunAtSize,1,0)';

if ischar(coT1), coT1=cellstr(coT1); end
if size(coT1,2)~=1, coT1=coT1'; end

tsp=[-0.5; -0.5872; -0.4066];

try
for t=1:size(coT1,1)
    [d,f,eh]=fileparts(coT1{t});    
    V=spm_vol(coT1{t});
    for i=1:3
        V.mat(i,4)=tsp(i)*V.dim(i)*V.mat(i,i);
    end
    if ~exist([d filesep prefix f '.hdr'],'file')
        copyfile([d filesep f '.hdr'],[d filesep prefix f '.hdr']);
    end
    spm_create_vol(V);
    disp(f);
    if RunAtSize~=0
        eval(RunAtCmd);
    end
end


catch
    if emlopt
        email_status(0,'img_auto_translate',['Job failed at t=' int2str(t)]);
    end
end

if emlopt
    if t==size(coT1,1)
        email_status(1,'img_auto_translate.m',['Job completed successfully.']);
    end
end
end