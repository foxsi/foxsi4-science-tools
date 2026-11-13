# What I did

I used Iain Hannah's code to make these files that can be found in the [`fvth_stuff` repo](https://github.com/ianan/fvth_stuff)

In Iain's IDL folder (/Users/kris/Downloads/flares-for-fierce/fvth_stuff-main/idl) I ran:

- `contfile=setup_chianti_cont('CHIANTI',0.1,250,$
    genxfile='chianti_cont_01_250_unity_v101_t41_kris.geny',$
    savfile='chianti_cont_01_250_unity_v101_t41_kris.sav')`
- `linefile=setup_chianti_lines('CHIANTI',0.1,12,ntemp=41,$
    genxfile='chianti_lines_01_12_unity_v101_t41_kris.geny',$
    savfile='chianti_lines_01_12_unity_v101_t41_kris.sav')`

I then just used the `.sav` files.

Other CHIANTI grid files can be found [here](https://soho.nascom.nasa.gov/solarsoft/packages/xray/dbase/chianti/).

- For speed, in the code I think I actually use `chianti_lines_01_12_unity_v901_t41.sav` from the link for the emission line file in `../phot_spec.py` instead of `chianti_lines_01_12_unity_v101_t41_kris.sav` from above.
