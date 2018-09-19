function [arr_names,sz]=Dir2Arr(srcpth,extns,omit)
%Dir2Arr sends dir search filenames to matlab array variable (like its pre-
%decessor, DirToArr), but with extended functionality. srcpth, extns and 
%omit can be all be (cell) arrays of text. Returns text array with
%spaced-padding at end



MSL=500;

if nargin > 2
    omit_arr=Dir2Arr(srcpth,omit);
    if isempty(omit_arr)
        omit_list={};
    else
        O=size(omit_arr);
        for o=1:O(1)
            [~,f,e]=fileparts(omit_arr(o,:));
            %[blah,f,e]=fileparts(omit_arr(o,:));
            omit_arr(o,:)=[f e blanks(MSL-length([f e]))];
        end
        omit_list=cellstr(omit_arr);
    end
else
    omit_list={};
end

arr_names=[];

if ~iscellstr(srcpth), srcpth=cellstr(srcpth); end
if ~iscellstr(extns), extns=cellstr(extns); end
    
for p=1:length(srcpth)
    read_dir=srcpth{p};
    for e=1:length(extns)
        G=dir([read_dir filesep extns{e}]);
        for g=1:length(G)
            if ~ismember(G(g).name,omit_list) && ~strcmp(G(g).name(1),'.')
                fname=[read_dir filesep G(g).name];
                arr_names=[arr_names; fname blanks(MSL-length(fname))];
            end
        end    
    end
    %disp(p);
end 


S=size(arr_names); sz=S(1);

disp(['Dir2Arr: ' int2str(sz) ' entries.']);

return;

end

