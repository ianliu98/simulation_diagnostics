;+
; NAME:
;    c20c_dtos_v2_adjust_sic_pall2007
;
; PURPOSE:
;    This procedure returns the sea ice concentration versus sea surface 
;    temperature function used in Pall (2007) and Pall et alii (2011).
;
; CATEGORY:
;    C20C dtos v2
;
; CALLING SEQUENCE:
;    c20c_dtos_v2_adjust_sic_pall2007, fit_tos_data=fit_tos_data, fit_sic_data=fit_sic_data, fit_lon=fit_lon, fit_lat=fit_lat, fit_time=fit_time, fit_type=fit_type
;
; INPUTS:
;    ---
;
; KEYWORD PARAMETERS:
;    FIT_LAT:  Returns a floating-point vector of length N_FIT_LAT=2 containing 
;        the latitude dimension for FIT_SIC_DATA and FIT_TOS_DATA.  This 
;        contains the mid-points of the two hemispheres, i.e. [-45,45].
;    FIT_LON:  Returns a floating-point vector of length N_FIT_LON=1 containing 
;        the longitude dimension for FIT_SIC_DATA and FIT_TOS_DATA.  This is 
;        mostly a dummy variable but is useful when this procedure is called by 
;        c20c_dtos_v2_adjust_sic.pro.
;    FIT_SIC_DATA:  Returns a 2*N_FIT_LON*N_FIT_LAT*N_FIT_TIME floating point 
;        array defining the sea ice concentration values at the end-points of 
;        the linear fit between full-ice and no-ice cover for both the Northern 
;        and Southern Hemispheres.  In units of %.  The first dimension 
;        corresponds to the 2 points, while N_FIT_LON=N_FIT_TIME=1 and are 
;        included for compatibility with c20c_dtos_v2_adjust_sic.pro.
;    FIT_TIME:  Returns a floating-point vector of length N_FIT_TIME=1 
;        containing the time dimension for FIT_SIC_DATA and FIT_TOS_DATA.  This 
;        os mostly a dummy variable but is useful when this procedure is called 
;        by c20c_dtos_v2_adjust_sic.pro.
;    FIT_TOS_DATA:  Returns a 2*N_FIT_LON*N_FIT_LAT*N_FIT_TIME floating point 
;        array defining the sea surface temperature values at the end-points of 
;        the linear fit between full-ice and no-ice cover for both the Northern 
;        and Southern Hemispheres.  In units of Kelvin.  The first dimension 
;        corresponds to the 2 points, while N_FIT_LON=N_FIT_TIME=1 and are 
;        included for compatibility with c20c_dtos_v2_adjust_sic.pro.
;    FIT_TYPE:  Returns a description for c20c_dtos_v2_adjust_sic.pro on how 
;        to interpret the data returned by this procedure.
;
; OUTPUTS:
;    FIT_LAT, FIT_LON, FIT_SIC_DATA, FIT_TIME, FIT_TOS_DATA, FIT_TYPE
;
; USES:
;    ---
;
; PROCEDURE:
;    This defines the sea ice adjustment function developed in:
;      * Pall, P.  2007.  Constraints on, and attribution of, changes in 
;        extreme precipitation under climate change.  Ph.D. Thesis, St. Cross 
;        College, University of Oxford, 187pp.) and Pall P.
;      * Pall, P., T. Aina, D. A. Stone, P. A. Stott, T. Nozawa, A. G. J. 
;        Hilberts, D. Lohmann, and M. R. Allen.  2011.  Anthropogenic 
;        greenhouse gas contribution to flood risk in England and Wales in 
;        Autumn 2000.  Nature, 470, 382-385.
;    Centre-of-mass values are from Pardeep Pall, personal communication.
;
; EXAMPLE:
;    See c20c_dtos_v2_adjust_sic.pro.
;
; MODIFICATION HISTORY:
;    Written by:  Daithi Stone (dastone@runbox.com), 2017-06-18, as 
;        c20c_adjust_sic_pall.pro.
;    Modified:  DAS, 2017-10-13 (Extracted from c20c_adjust_sic_pall.pro into 
;        c20c_dtos_v2_adjust_sic_pall2007.pro;  added to IDL routine library)
;-

;***********************************************************************

PRO C20C_DTOS_V2_ADJUST_SIC_PALL2007, $
    FIT_TOS_DATA=fit_tos_data, FIT_SIC_DATA=fit_sic_data, $
    FIT_LON=fit_lon, FIT_LAT=FIT_LAT, FIT_TIME=fit_time, $
    FIT_TYPE=fit_type

;***********************************************************************
; Constants

; Define the freezing point
tos_freeze = 271.35
; Define full-ice and no-ice coverage
sic_full = 100.
sic_none = 0.

; Define the centre of mass of the partially ice-covered grid cells
centre_tos_nh = 271.714
centre_sic_nh = 100. * 0.896820
centre_tos_sh = 271.745
centre_sic_sh = 100. * 0.769180

; Define the type of fit
fit_type = 'straight line'

;***********************************************************************
; Calculate the linear functions

; Define the dimensions of the output
; (two points each for the Northern and Southern Hemispheres, no seasons)
n_fit_point = 2
fit_lon = 0
n_fit_lon = n_elements( fit_lon )
fit_lat = [ -45., 45. ]
n_fit_lat = n_elements( fit_lat )
fit_time = ''
n_fit_time = n_elements( fit_time )
; Initialise output arrays
fit_tos_data = !values.f_nan $
    * fltarr( n_fit_point, n_fit_lon, n_fit_lat, n_fit_time )
fit_sic_data = fit_tos_data

; Add the full-ice/freezing point to both hemispheres
fit_sic_data[0,*,*,*] = sic_full
fit_tos_data[0,*,*,*] = tos_freeze
; Calculate the zero-ice point
fit_sic_data[1,*,*,*] = sic_none
fit_tos_data[1,*,0,*] = ( centre_tos_nh - fit_tos_data[0,*,0,*] ) $
    / ( centre_sic_nh - fit_sic_data[0,*,0,*] ) $
    * ( fit_sic_data[1,*,0,*] - fit_sic_data[0,*,0,*] ) + fit_tos_data[0,*,0,*]
fit_tos_data[1,*,1,*] = ( centre_tos_sh - fit_tos_data[0,*,0,*] ) $
    / ( centre_sic_sh - fit_sic_data[0,*,0,*] ) $
    * ( fit_sic_data[1,*,0,*] - fit_sic_data[0,*,0,*] ) + fit_tos_data[0,*,0,*]

;***********************************************************************
; The end

return
END

