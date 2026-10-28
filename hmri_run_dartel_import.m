function out = hmri_run_dartel_import(job)
% Execute the hMRI external-segment import module.

reference = char(job.reference{1});
segments  = job.segments;
outdir    = char(job.outdir{1});

if ~exist(outdir, 'dir')
    mkdir(outdir);
end

spm('Defaults', 'fMRI');
spm_jobman('initcfg');

% use the hMRI specific TPMs.
fn_tpm = hmri_get_defaults('proc.TPM');
% get the number of Gaussians per tissue class
nGauss = hmri_get_defaults('proc.nGauss');

%get bounding box of original data
[bb, vx] = spm_get_bbox(segments{1}, 'fv'); %maybe check if all segments have the same bb?

% Standard SPM unified segmentation.
% the _seg8 file is needed

matlabbatch = {};
matlabbatch{1}.spm.spatial.preproc.channel.vols = {reference};
% matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
% matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
%%hard code the default settings hMRI uses for now, see tbx_scfg_hmri_proc_US.m
%bias cutoff = 'no bias' correction
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = Inf;

for k = 1:6
    matlabbatch{1}.spm.spatial.preproc.tissue(k).tpm = ...
        {[fn_tpm ',' num2str(k)]};
    matlabbatch{1}.spm.spatial.preproc.tissue(k).ngaus = nGauss(k);
    matlabbatch{1}.spm.spatial.preproc.tissue(k).native = [0 0];%[1 1]; I need the result, not the warped tissues
    matlabbatch{1}.spm.spatial.preproc.tissue(k).warped = [0 0];
end
%% just for debugging: make 1 imported segment
% matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 1];
%% end of debugging, MN
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = ...
    [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];
%do not pass along the bounding box, revert to the default bb (as the hMRI
%deafult does)
% matlabbatch{1}.spm.spatial.preproc.warp.bb = bb; % also added by MN
% matlabbatch{1}.spm.spatial.preproc.warp.vox =vx; % also added by MN, but seems to be overwitten during runtime
%use defaut SPM bounding box, voxel size and regularisation (like hMRI does,tbx_scfg_hmri_proc_US.m , bb.def and vox.val)
% for now hardcode voxel size to 1, it seems like I can't pass it on to
% preproc in a clean way?
% !!!! can be made avaialble as user settings later
matlabbatch{1}.spm.spatial.preproc.warp.bb = spm_get_defaults('normalise.write.bb');
outvox = [1 1 1];
matlabbatch{1}.spm.spatial.preproc.warp.vox = mean(outvox);
matlabbatch{1}.spm.spatial.preproc.warp.reg=[0	0	0.100000000000000	0.0100000000000000	0.0400000000000000];
spm_jobman('run', matlabbatch);

% Find the seg8.mat file produced by SPM.
[refFolder, refName] = fileparts(reference);
seg8File = fullfile(refFolder, [refName '_seg8.mat']);

if ~exist(seg8File, 'file')
    error('Expected SPM segmentation info file not found: %s', seg8File);
end

%move the seg8file to the output folder:
movefile(seg8File,  outdir);
seg8File = fullfile(outdir, [refName '_seg8.mat']);

% Use its geometry to resample external segments.
out.outputs = hmri_dartel_import( seg8File, segments, outvox, outdir);

end