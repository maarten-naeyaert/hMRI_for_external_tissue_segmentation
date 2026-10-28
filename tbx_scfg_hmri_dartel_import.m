function dartel_import = tbx_scfg_hmri_dartel_import

reference = cfg_files;
reference.tag = 'reference';
reference.name = 'Reference anatomical image';
reference.filter = 'image';
reference.ufilter = '.*';
reference.num = [1 1];

segments = cfg_files;
segments.tag = 'segments';
segments.name = 'existing segment images';
segments.filter = 'image';
segments.ufilter = '.*';
segments.num = [1 Inf];

outdir = cfg_files;
outdir.tag = 'outdir';
outdir.name = 'Output directory';
outdir.filter = 'dir';
outdir.num = [1 1];

dartel_import = cfg_exbranch;
dartel_import.tag = 'dartel_import';
dartel_import.name = 'Import existing segments for DARTEL';
dartel_import.val = {reference, segments, outdir};
dartel_import.prog = @hmri_run_dartel_import;
dartel_import.help = { ...
    'Resample existing segment images onto the SPM DARTEL-import grid.'};

end