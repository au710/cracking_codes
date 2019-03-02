#### import the simple module from the paraview
from paraview.simple import *
import numpy as np
import os
#from common.constants import Kap_raw,Ga_raw,T_raw,K_raw,G_raw,alp_raw,phi_raw
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()


for fol_name in os.listdir("."):
    if os.path.isdir(fol_name):
        try:
            fil_name=fol_name
            alphapvd = PVDReader(FileName='./'+fol_name+'/d.pvd')
            animationScene1 = GetAnimationScene()

            # update animation scene based on data timesteps
            animationScene1.UpdateAnimationUsingDataTimeSteps()

            # get active view
            renderView1 = GetActiveViewOrCreate('RenderView')
            renderView1.Background = [1,1,1]  #white
            f24LUT = GetColorTransferFunction('f24')
            alphapvdDisplay = Show(alphapvd, renderView1)
            # get opacity transfer function/opacity map for 'f24'
            #f24PWF = GetOpacityTransferFunction('f26')
            # Rescale transfer function
            f24LUT.RescaleTransferFunction(0.0, 1.0)
            # Rescale transfer function
            #f24PWF.RescaleTransferFunction(0.0, 1.0)
            # hide color bar/color legend
            alphapvdDisplay.SetScalarBarVisibility(renderView1, False)
            # reset view to fit data bounds
            renderView1.ResetCamera(-1.0, 1.0, -1.0, 1.0, 0.0, 0.0)
            # current camera placement for renderView1
            renderView1.OrientationAxesVisibility = 0
            # get layout
            viewLayout1 = GetLayout()

            ii=0
            while ii<=8:
                   # save screenshot
                   SaveScreenshot('./'+fil_name+'_%d.png'%(ii), renderView1, ImageResolution=[1170, 927])
                   animationScene1.GoToNext()
                   ii += 1

            animationScene1.GoToLast()
            # destroy alphapvd
            Delete(alphapvd)
            del alphapvd
        except:
            f_err = open('nocrack.txt','a+')
            f_err.write(fol_name+'\n')
