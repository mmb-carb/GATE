# Pairing down the Input Files

## METCRO3D

As it happens, we only need the "ZF" paramter from the METCRO3D file, so we can remove all the other parameters to save space. And the small hourly variation of these grid heights does not seem to affect the GATE model, so I average these values over time.  This all speeds up the processing considerably and saves a lot of space. Though, I suppose, it is totally optional.

In case you want to do the same in the future, with your own input files, here are the Linux NCO commands I used to edit this file:

    ncks -v TFLAG,JACOBF,JACOBM,DENSA_J,WHAT_JD,TA,QV,PRES,DENS,ZH,QC,QR,QI,QS,QG -x METCRO3D_2010_12 METCRO3D_2010_12_ZF
    ncra -F -d time,1,45,1 METCRO3D_2010_12_ZF METCRO3D_2010_12_ZF_AVG

