% function outputFiles = hmri_dartel_import(seg8File, segments, outdir)
function varargout = hmri_dartel_import(seg8File, segments, outvox,outdir)
% Resample external segment images to the SPM DARTEL-import grid.
% seg8file is the output seg8file of the US result of a reference image
% The segment images must have identical voxel geometry to the reference image used to make the seg8File.
% based on spm_preproc_write8.m

%%  read the result from '_seg8.mat'file (= result of spm_preproc8)
 seg8File=load(seg8File);

%% get the segments which need to be imported
for k = 1:numel(segments)
    src{k} = char(segments{k});
    Vsrc{k} = spm_vol(src{k});
    src{k} = src{k}(1:end-2);
%     [unused, src{k}, temp] = fileparts(src{k});
%     src{k}=[src{k} temp]; 
    [unused, src{k}] = fileparts(src{k}); % no need to add .nii to the name
    if numel(Vsrc{k}) ~= 1
        error('Each segment must be a single-volume image: %s', src{k});
    end
end

%% get name of reference
[unused,nam] = fileparts(seg8File.image(1).fname);

% Read essentials from tpm (it will be cleared later)
tpm = seg8File.tpm;
if ~isstruct(tpm) || ~isfield(tpm, 'bg1')
    tpm = spm_load_priors8(tpm);
end
d1      = size(tpm.dat{1});
d1      = d1(1:3);
M1      = tpm.M;

%% get bounding box ==> line 48 - 75 in spm_preproc_write8
% get the bounding box and voxel size of the original segment images and of
% the tpm
%include vx or not?
% vx = NaN;
% [bb,~] = spm_get_bbox(segments{1}, 'fv'); %maybe check if all segments have the same bb?
[bb,vx] = spm_get_bbox(segments{1}, 'fv'); %maybe check if all segments have the same bb?
vx=outvox;

bb=[-78	-112	-70
78	76	85];

vx = mean(vx); % from C:\matlab_work\hMRI-toolbox-master\hmri_run_proc_US.m , line 71
[bb1,vx1] = spm_get_bbox(tpm.V(1), 'old');

bb(~isfinite(bb)) = bb1(~isfinite(bb));
if ~isfinite(vx), vx = abs(prod(vx1))^(1/3); end
bb(1,:) = vx.*round(bb(1,:)/vx); %bounding box approximated in voxel dimensions of the tpm
bb(2,:) = vx.*round(bb(2,:)/vx); %bounding box approximated in voxel dimensions of the tpm

%matrix size of the bounding box using vox
odim    = abs(round((bb(2,1:3)-bb(1,1:3))/vx))+1;

mm  = [[bb(1,1) bb(1,2) bb(1,3)
    bb(2,1) bb(1,2) bb(1,3)
    bb(1,1) bb(2,2) bb(1,3)
    bb(2,1) bb(2,2) bb(1,3)
    bb(1,1) bb(1,2) bb(2,3)
    bb(2,1) bb(1,2) bb(2,3)
    bb(1,1) bb(2,2) bb(2,3)
    bb(2,1) bb(2,2) bb(2,3)]'; ones(1,8)];
vx3 = [[1       1       1
    odim(1) 1       1
    1       odim(2) 1
    odim(1) odim(2) 1
    1       1       odim(3)
    odim(1) 1       odim(3)
    1       odim(2) odim(3)
    odim(1) odim(2) odim(3)]'; ones(1,8)];
mat    = mm/vx3;

%% get d: line 87 in spm_preproc_write8
d = seg8File.image(1).dim(1:3); %res is the reoriented T1w image
[x1,x2] = ndgrid(1:d(1),1:d(2),1);

%% get bsplinc coeffiecients, line 132-136
prm     = [3 3 3 0 0 0];
Coef    = cell(1,3);
Coef{1} = spm_bsplinc(seg8File.Twarp(:,:,:,1),prm);
Coef{2} = spm_bsplinc(seg8File.Twarp(:,:,:,2),prm);
Coef{3} = spm_bsplinc(seg8File.Twarp(:,:,:,3),prm);


%% get y (inverse transform): line 153 - 202, selected lines
y = zeros([seg8File.image(1).dim(1:3),3],'single');
x3 = 1:d(3);
spm_progress_bar('init',length(x3),['Working on ' nam],'Planes completed');

M = M1\seg8File.Affine*seg8File.image(1).mat;
for z=1:length(x3)
    
    [t1,t2,t3] = defs(Coef,z,seg8File.MT,prm,x1,x2,x3,M);
    if exist('y','var')
        % If needed later, save in variable y
        y(:,:,z,1) = t1;
        y(:,:,z,2) = t2;
        y(:,:,z,3) = t3;
    end
    spm_progress_bar('set',z);
    
end
spm_progress_bar('clear');

%% get M0: line 206
M0 = seg8File.image(1).mat;


%% define outpur variables: line 207 & 210
%  rc   = cell(size(tc,1),1);
%  mat0 = M0;


%% get cls: line 211
%   rc   = cell(size(tc,1),1);
for i= 1:numel(Vsrc)
    temp = load_untouch_nii(Vsrc{i}.private.dat.fname);  % the 'tissue classes' which need to be imported, so the files which I want to import?
    %     cls{i}= temp.img*temp.hdr.dime.scl_slope+temp.hdr.dime.scl_inter; %if int can give issues
        cls{i}= single(temp.img)*temp.hdr.dime.scl_slope+temp.hdr.dime.scl_inter;
        maxRange(i) = temp.hdr.dime.glmax*temp.hdr.dime.scl_slope+temp.hdr.dime.scl_inter;
%     cls{i}=uint8(round(255 *temp.img));
end

%  [rc,mat0] = write_imported_tissues(tc(:,[2 6]),d,M0,y,M1,cls,mat,vx,res,pth,nam,odim);
[rc,mat0] = write_imported_tissues(d,M0,y,M1,cls,mat,vx,seg8File,outdir,src,odim,maxRange);

outputFiles=rc;
varargout = {{cls,M0},{rc,M1,mat0}};

%==========================================================================
% function [rc,mat0] = write_imported_tissues(d,M0,y,M1,cls,mat,vx,src,pth,filler,odim)
% (based on spm_preproc_write8.m)
%==========================================================================
function [rc,mat0] = write_imported_tissues(d,M0,y,M1,cls,mat,vx,res,pth,src,odim,maxRange)
% "Imported" tissue class images

% Generate mm coordinates of where deformations map from
x      = affind(rgrid(d),M0);

% Generate mm coordinates of where deformation maps to
y1     = affind(y,M1);

% Procrustes analysis to compute the closest rigid-body
% transformation to the deformation, weighted by the
% interesting tissue classes.
% ind        = find(any(tc,2)); % Saved tissue classes
% [unused,R] = spm_get_closest_affine(x,y1,single(cls{1})/255);
[unused,R] = spm_get_closest_affine(x,y1,single(cls{1})/maxRange(1));
clear x y1

mat0   = R\mat; % Voxel-to-world of original image space
rc     = cell(size(cls));
fwhm   = max(vx./sqrt(sum(res.image(1).mat(1:3,1:3).^2))-1,0.01);
for k1=1:numel(cls)
    nam=src{k1};
    thisMaxRange=maxRange(k1);
    % Low pass filtering to reduce aliasing effects in downsampled images,
    % then reslice and write to disk
    tmp1    = decimate(single(cls{k1}),fwhm);
    tmp2    = zeros(odim,'single');
    for i=1:odim(3)
        tmp2(:,:,i)   = spm_slice_vol(tmp1,M0\mat0*spm_matrix([0 0 i]),odim(1:2),[1,NaN])/thisMaxRange; 
%         tmp2(:,:,i)   = spm_slice_vol(tmp1,M0\mat0*spm_matrix([0 0 i]),odim(1:2),[1,NaN])/255; %second oldest
        %         tmp2(:,:,i)   = spm_slice_vol(tmp1,M0\mat0*spm_matrix([0 0 i]),odim(1:2),[1,NaN]); %oldest
    end
    clear tmp1
    Ni      = nifti;
    Ni.dat  = file_array(fullfile(pth,['rc_', nam, '.nii']),...
        odim,...
        [spm_type('float32') spm_platform('bigend')],...
        0,1,0);
    Ni.mat         = mat;
    Ni.mat_intent  = 'Aligned';
    Ni.mat0        = mat0;
    Ni.mat0_intent = 'Aligned';
    Ni.descrip     = ['Imported Tissue ' num2str(k1)];
    create(Ni);
    Ni.dat(:,:,:) = tmp2;
    
    clear tmp2
end

%==========================================================================
% function [x1,y1,z1] = defs(sol,z,MT,prm,x0,y0,z0,M)
% (from spm_preproc_write8.m)
%==========================================================================
function [x1,y1,z1] = defs(sol,z,MT,prm,x0,y0,z0,M)
iMT = inv(MT);
x1  = x0*iMT(1,1)+iMT(1,4);
y1  = y0*iMT(2,2)+iMT(2,4);
z1  = (z0(z)*iMT(3,3)+iMT(3,4))*ones(size(x1));

% Eliminate NaNs (see email from Pratik on 01/09/23)
x1  = min(max(x1,1),size(sol{1},1));
y1  = min(max(y1,1),size(sol{1},2));
z1  = min(max(z1,1),size(sol{1},3));

x1a = x0    + spm_bsplins(sol{1},x1,y1,z1,prm);
y1a = y0    + spm_bsplins(sol{2},x1,y1,z1,prm);
z1a = z0(z) + spm_bsplins(sol{3},x1,y1,z1,prm);
x1  = M(1,1)*x1a + M(1,2)*y1a + M(1,3)*z1a + M(1,4);
y1  = M(2,1)*x1a + M(2,2)*y1a + M(2,3)*z1a + M(2,4);
z1  = M(3,1)*x1a + M(3,2)*y1a + M(3,3)*z1a + M(3,4);
return;

%==========================================================================
% function x = rgrid(d)
%==========================================================================
function x = rgrid(d)
x = zeros([d(1:3) 3],'single');
[x1,x2] = ndgrid(single(1:d(1)),single(1:d(2)));
for i=1:d(3)
    x(:,:,i,1) = x1;
    x(:,:,i,2) = x2;
    x(:,:,i,3) = single(i);
end
return;

%==========================================================================
% function y1 = affind(y0,M)
%==========================================================================
function y1 = affind(y0,M)
y1 = zeros(size(y0),'single');
for d=1:3
    y1(:,:,:,d) = y0(:,:,:,1)*M(d,1) + y0(:,:,:,2)*M(d,2) + y0(:,:,:,3)*M(d,3) + M(d,4);
end
return;

%==========================================================================
% function dat = decimate(dat,fwhm)
%==========================================================================
function dat = decimate(dat,fwhm)
% Convolve the volume in memory (fwhm in voxels).
lim = ceil(2*fwhm);
x  = -lim(1):lim(1); x = spm_smoothkern(fwhm(1),x); x  = x/sum(x);
y  = -lim(2):lim(2); y = spm_smoothkern(fwhm(2),y); y  = y/sum(y);
z  = -lim(3):lim(3); z = spm_smoothkern(fwhm(3),z); z  = z/sum(z);
i  = (length(x) - 1)/2;
j  = (length(y) - 1)/2;
k  = (length(z) - 1)/2;
spm_conv_vol(dat,dat,x,y,z,-[i j k]);
return;
