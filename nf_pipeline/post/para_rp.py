#### import the simple module from the paraview
from paraview.simple import *
import numpy as np
from common.constants import Kap_raw,Ga_raw,T_raw,K_raw,G_raw,alp_raw,phi_raw
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

#-----material constants-----
Gc = Ga_raw;
kappa=Kap_raw;

#-------evap parameters--------
Theta=T_raw;
K=K_raw;
G = G_raw;
alp=alp_raw
phi_s0=phi_raw;

fol_name = "SL-theta%.2f-k%.2f-G%.2f-Gc%.2f-kap%.2f-phi%.2f"%(Theta,K,G,Gc,kappa,phi_s0)
fil_name = "fin-theta%.2f-k%.2f-G%.2f-Gc%.2f-kap%.2f-phi%.2f.png"%(Theta,K,G,Gc,kappa,phi_s0)

alphapvd = PVDReader(FileName='./'+fol_name+'/alpha.pvd')
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
renderView1.Background = [1,1,1]  #white
f24LUT = GetColorTransferFunction('f25')
alphapvdDisplay = Show(alphapvd, renderView1)

# get opacity transfer function/opacity map for 'f24'
f24PWF = GetOpacityTransferFunction('f25')

animationScene1.GoToLast()

# Rescale transfer function
f24LUT.RescaleTransferFunction(0.0, 1.0)

# Rescale transfer function
f24PWF.RescaleTransferFunction(0.0, 1.0)

# hide color bar/color legend
alphapvdDisplay.SetScalarBarVisibility(renderView1, False)
# reset view to fit data bounds
renderView1.ResetCamera(-1.0, 1.0, -1.0, 1.0, 0.0, 0.0)

# current camera placement for renderView1
renderView1.OrientationAxesVisibility = 0
# get layout
viewLayout1 = GetLayout()

# save screenshot
SaveScreenshot('./'+fil_name, renderView1, ImageResolution=[1170, 927])

# hide data in view
Hide(alphapvd, renderView1)

# destroy alphapvd
Delete(alphapvd)
del alphapvd

