f = addfile(file_path, "r")

time_step = 28

figure_name = systemfunc("basename " + file_path + " .nc")

wks = gsn_open_wks("pdf", figure_name)

res = True
res@gsnDraw = False
res@gsnFrame = False
res@cnLinesOn = True
res@cnFillOn = True
res@cnFillPalette = "GMT_panoply"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -3.0e-8
res@cnMaxLevelValF =  3.0e-8
res@cnLevelSpacingF = 2.0e-9
res@mpOutlineOn = False
; res@mpMinLatF = 10.0
; res@mpMaxLatF = 80.0
res@mpCenterLonF = 180.0
res@vpHeightF = 0.5

plot0 = gsn_csm_contour_map(wks, f->pv(time_step,:,:), res)

draw(plot0)
frame(wks)

print("Run the following command to postprocess figure:")
print("pdfcrop " + figure_name + ".pdf && convert -density 300 " + figure_name + "-crop.pdf " + figure_name + ".png")
