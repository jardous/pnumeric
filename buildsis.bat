@REM devices -setdefault @S60_3rd_MR:com.nokia.s60 REM @S60_2nd_FP3:com.nokia.series60
bldmake clean
bldmake bldfiles
abld build gcce urel
makesis -v pnumeric.pkg
signsis pnumeric.sis pnumeric_s60_selfsign.sis c:\keys\selfsign\mycert.cer c:\keys\selfsign\mykey.key tcs60fs
del pnumeric.sis
