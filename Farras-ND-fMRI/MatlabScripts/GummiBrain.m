function GummiBrain(ColorIndex,AtlasName,MAP,PlotHemi)
%ColorIndex is the intensity of the color you want to plot for each
%region, you need to provide this. It will have one number per region in
%your atlas.

%AtlasName is the atlased image that you want to plot the brain on 
%(this is the atlased T1 image in the "Atlased" subfolder of the T1scan 
%subdirectory).
at = spm_read_vols(spm_vol(AtlasName));
at(isnan(at)) = 0;
ss = size(at);
roinums = setdiff(unique(at(:)),0);
roinums = roinums(1:90); % farras, 08/06/2012, 21:53

%MAP: the colormap you want to use. just check the help for colormap if you
%want to change this.
if nargin<3
    MAP = colormap(jet(length(roinums))); % colormap(hot(length(roinums)));
end
%If you want to plot the whole brain, use the default setting or input
%PlotHemi as a vector of all ones; If you want to plot just the left or right hemisphere,
%input PlotHemi as a binary vector indicating the regions you want to plot
% with 1's.
if nargin<4
    PlotHemi = ones(length(roinums),1);
end
flag = 0;
if length(unique(PlotHemi)) == 1;
    flag = 1;
end

for i = roinums';
    if flag
        plot3d = at(:)==i;
        plot3d = double(reshape(plot3d,ss));
        myvol_render(plot3d,[1 1 1],MAP,ColorIndex(i))
        hold on;
    elseif PlotHemi(i)
        plot3d = at(:)==i;
        plot3d = double(reshape(plot3d,ss));
        myvol_render(plot3d,[1 1 1],MAP,ColorIndex(i))
        hold on;
    end
end
lightangle(25,30); % was lightangle(45,30)
lighting phong;

function [] = myvol_render(seg,side_len,Map,Color)
% seg is a binary 3D volume to render.  side_len is the aspect ratio of
% each dimension. that we want to render.
[h,w,d] = size(seg);
Ds = smooth3(seg*100);
hiso = patch(isosurface(Ds,5),'FaceColor',Map(Color,:),'EdgeColor','none');
isonormals(Ds,hiso);
view(-60,30); % was view(35,30)
%axis([1 w 1 h 1 d]);
axis off
side_len = side_len([2 1 3]);
daspect(1./side_len);