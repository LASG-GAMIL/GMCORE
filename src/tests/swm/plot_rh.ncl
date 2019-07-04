f = addfile(file_path, "r")

time_step = 100

figure_name = systemfunc("basename " + file_path + " .nc")

wks = gsn_open_wks("pdf", figure_name)

res = True
res@gsnDraw = False
res@gsnFrame = False
res@cnLinesOn = True
res@cnFillOn = True
res@cnFillPalette = "GMT_panoply"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 8000
res@cnMaxLevelValF = 10500
res@cnLevelSpacingF = 200
res@mpOutlineOn = False

plot0 = gsn_csm_contour_map(wks, f->h(time_step,:,:), res)

res@gsnPolar = "NH"
res@mpMinLatF = 0
res@lbLabelBarOn = False
res@gsnLeftString = ""
res@gsnRightString = ""
res@vpWidthF = 0.45

plot1 = gsn_csm_contour_map(wks, f->h(time_step,:,:), res)

res_panel = True

gsn_panel(wks, (/plot0,plot1/), (/1,2/), res_panel)

print("Run the following command to postprocess figure:")
print("pdfcrop " + figure_name + ".pdf && convert -density 300 " + figure_name + "-crop.pdf " + figure_name + ".png")
