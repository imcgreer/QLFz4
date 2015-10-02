
#Forced Photometry DeepForcedSource.csv column names
#Implied first column is deepForcedSourceId.
FORCEDSOURCE_COLS = [
    'coord_ra', 'coord_decl', 'parent', 'flags_badcentroid', 'centroid_sdss_x', 'centroid_sdss_y',
    'centroid_sdss_xVar', 'centroid_sdss_xyCov', 'centroid_sdss_yVar', 'centroid_sdss_flags', 'centroid_gaussian_x',
    'centroid_gaussian_y', 'centroid_gaussian_xVar', 'centroid_gaussian_xyCov', 'centroid_gaussian_yVar',
    'centroid_gaussian_flags', 'centroid_naive_x', 'centroid_naive_y', 'centroid_naive_xVar', 'centroid_naive_xyCov',
    'centroid_naive_yVar', 'centroid_naive_flags', 'flags_pixel_edge', 'flags_pixel_interpolated_any',
    'flags_pixel_interpolated_center', 'flags_pixel_saturated_any', 'flags_pixel_saturated_center',
    'flags_pixel_cr_any', 'flags_pixel_cr_center', 'flags_pixel_bad', 'shape_sdss_Ixx', 'shape_sdss_Iyy',
    'shape_sdss_Ixy', 'shape_sdss_IxxVar', 'shape_sdss_IxxIyyCov', 'shape_sdss_IxxIxyCov', 'shape_sdss_IyyVar',
    'shape_sdss_IyyIxyCov', 'shape_sdss_IxyVar', 'shape_sdss_flags', 'shape_sdss_centroid_x',
    'shape_sdss_centroid_y', 'shape_sdss_centroid_xVar', 'shape_sdss_centroid_xyCov',
    'shape_sdss_centroid_yVar', 'shape_sdss_centroid_flags', 'shape_sdss_flags_unweightedbad',
    'shape_sdss_flags_unweighted', 'shape_sdss_flags_shift', 'shape_sdss_flags_maxiter',
    'flux_gaussian', 'flux_gaussian_err', 'flux_gaussian_flags', 'flux_gaussian_psffactor',
    'flux_gaussian_flags_psffactor', 'flux_naive', 'flux_naive_err', 'flux_naive_flags',
    'flux_psf', 'flux_psf_err', 'flux_psf_flags', 'flux_psf_psffactor', 'flux_psf_flags_psffactor',
    'flux_sinc', 'flux_sinc_err', 'flux_sinc_flags', 'correctfluxes_apcorr',
    'correctfluxes_apcorr_flags', 'centroid_record_x', 'centroid_record_y',
    'classification_extendedness', 'refFlux', 'refFlux_err', 'objectId', 'coord_raVar',
    'coord_radeclCov', 'coord_declVar', 'exposure_id', 'exposure_filter_id', 'exposure_time',
    'exposure_time_mid', 'cluster_id', 'cluster_coord_ra', 'cluster_coord_decl']

#Coadd Photometry DeepSource.csv column names
#Implied first column is deepSourceId
DEEPSOURCE_COLS = [
    'coord_ra', 'coord_decl', 'parent', 'calib_detected', 'flags_negative', 'deblend_nchild',
    'deblend_deblended_as_psf', 'deblend_psf_center_x', 'deblend_psf_center_y', 'deblend_psf_flux', 'deblend_too_many_peaks',
    'deblend_failed', 'flags_badcentroid', 'centroid_sdss_x', 'centroid_sdss_y', 'centroid_sdss_xVar', 'centroid_sdss_xyCov',
    'centroid_sdss_yVar', 'centroid_sdss_flags', 'centroid_gaussian_x', 'centroid_gaussian_y', 'centroid_gaussian_xVar',
    'centroid_gaussian_xyCov', 'centroid_gaussian_yVar', 'centroid_gaussian_flags', 'centroid_naive_x', 'centroid_naive_y',
    'centroid_naive_xVar', 'centroid_naive_xyCov', 'centroid_naive_yVar', 'centroid_naive_flags', 'flags_pixel_edge',
    'flags_pixel_interpolated_any', 'flags_pixel_interpolated_center', 'flags_pixel_saturated_any',
    'flags_pixel_saturated_center', 'flags_pixel_cr_any', 'flags_pixel_cr_center', 'flags_pixel_bad',
    'shape_sdss_Ixx', 'shape_sdss_Iyy', 'shape_sdss_Ixy', 'shape_sdss_IxxVar', 'shape_sdss_IxxIyyCov',
    'shape_sdss_IxxIxyCov', 'shape_sdss_IyyVar', 'shape_sdss_IyyIxyCov', 'shape_sdss_IxyVar',
    'shape_sdss_flags', 'shape_sdss_centroid_x', 'shape_sdss_centroid_y', 'shape_sdss_centroid_xVar',
    'shape_sdss_centroid_xyCov', 'shape_sdss_centroid_yVar', 'shape_sdss_centroid_flags', 'shape_sdss_flags_unweightedbad',
    'shape_sdss_flags_unweighted', 'shape_sdss_flags_shift', 'shape_sdss_flags_maxiter', 'flux_gaussian', 'flux_gaussian_err',
    'flux_gaussian_flags', 'flux_gaussian_psffactor', 'flux_gaussian_flags_psffactor', 'flux_naive', 'flux_naive_err',
    'flux_naive_flags', 'flux_psf', 'flux_psf_err', 'flux_psf_flags', 'flux_psf_psffactor', 'flux_psf_flags_psffactor',
    'flux_sinc', 'flux_sinc_err', 'flux_sinc_flags', 'multishapelet_psf_inner_1', 'multishapelet_psf_outer_1',
    'multishapelet_psf_ellipse_Ixx', 'multishapelet_psf_ellipse_Iyy', 'multishapelet_psf_ellipse_Ixy',
    'multishapelet_psf_chisq', 'multishapelet_psf_integral', 'multishapelet_psf_flags', 'multishapelet_psf_flags_maxiter',
    'multishapelet_psf_flags_tinystep', 'multishapelet_psf_flags_constraint_r', 'multishapelet_psf_flags_constraint_q',
    'multishapelet_dev_flux', 'multishapelet_dev_flux_err', 'multishapelet_dev_flux_flags', 'multishapelet_dev_psffactor',
    'multishapelet_dev_flags_psffactor', 'multishapelet_dev_ellipse_Ixx', 'multishapelet_dev_ellipse_Iyy',
    'multishapelet_dev_ellipse_Ixy', 'multishapelet_dev_psffactor_ellipse_Ixx', 'multishapelet_dev_psffactor_ellipse_Iyy',
    'multishapelet_dev_psffactor_ellipse_Ixy', 'multishapelet_dev_chisq', 'multishapelet_dev_flags_maxiter',
    'multishapelet_dev_flags_tinystep', 'multishapelet_dev_flags_constraint_r', 'multishapelet_dev_flags_constraint_q',
    'multishapelet_dev_flags_largearea', 'multishapelet_exp_flux', 'multishapelet_exp_flux_err',
    'multishapelet_exp_flux_flags', 'multishapelet_exp_psffactor', 'multishapelet_exp_flags_psffactor',
    'multishapelet_exp_ellipse_Ixx', 'multishapelet_exp_ellipse_Iyy', 'multishapelet_exp_ellipse_Ixy',
    'multishapelet_exp_psffactor_ellipse_Ixx', 'multishapelet_exp_psffactor_ellipse_Iyy',
    'multishapelet_exp_psffactor_ellipse_Ixy', 'multishapelet_exp_chisq', 'multishapelet_exp_flags_maxiter',
    'multishapelet_exp_flags_tinystep', 'multishapelet_exp_flags_constraint_r', 'multishapelet_exp_flags_constraint_q',
    'multishapelet_exp_flags_largearea', 'multishapelet_combo_flux', 'multishapelet_combo_flux_err',
    'multishapelet_combo_flux_flags', 'multishapelet_combo_psffactor', 'multishapelet_combo_flags_psffactor',
    'multishapelet_combo_components_1', 'multishapelet_combo_components_2', 'multishapelet_combo_chisq',
    'correctfluxes_apcorr', 'correctfluxes_apcorr_flags', 'classification_extendedness', 'detect_is_patch_inner',
    'detect_is_tract_inner', 'detect_is_primary', 'coord_raVar', 'coord_radeclCov', 'coord_declVar', 'coadd_id',
    'coadd_filter_id']

#Coadd images DeepCoadd.csv column names
#Implied first column is deepCoaddId
DEEPCOADD_COLS = [
    'tract', 'patch', 'filterId', 'filterName', 'ra', 'decl', 'equinox', 'raDeSys', 'ctype1', 'ctype2',
    'crpix1', 'crpix2', 'crval1', 'crval2', 'cd1_1', 'cd1_2', 'cd2_1', 'cd2_2', 'corner1Ra', 'corner1Decl',
    'corner2Ra', 'corner2Decl', 'corner3Ra', 'corner3Decl', 'corner4Ra', 'corner4Decl', 'fluxMag0',
    'fluxMag0Sigma', 'matchedFwhm', 'measuredFwhm', 'path']

#Single Epoch images Science_Ccd_Exposure.csv column names
#Implied first column is scienceCcdExposureId
SCIENCE_CCD_EXPOSURE_COLS = [
    'run', 'camcol', 'filterId', 'field', 'filterName', 'ra', 'decl',
    'equinox', 'raDeSys', 'ctype1', 'ctype2', 'crpix1', 'crpix2',
    'crval1', 'crval2', 'cd1_1', 'cd1_2', 'cd2_1', 'cd2_2',
    'corner1Ra', 'corner1Decl', 'corner2Ra', 'corner2Decl',
    'corner3Ra', 'corner3Decl', 'corner4Ra', 'corner4Decl',
    'taiMjd', 'obsStart', 'expMidpt', 'expTime',
    'nCombine', 'binX', 'binY', 'fluxMag0', 'fluxMag0Sigma', 'fwhm', 'path']
