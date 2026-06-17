import numpy as np
import numpy.linalg as mat
import scipy as sp
import scipy.linalg as smat
# import cvxpy as cp
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

from svo_classes import *

# import sys
# sys.path.append('/Users/dan/Documents/code');
# from drawing import * 
# from mathz import * 
# sys.path.remove('/Users/dan/Documents/code');

# sys.path.append('/Users/dan/Documents/meeko/WORK/SCVX');
# from scvx_examples import *
# sys.path.remove('/Users/dan/Documents/meeko/WORK/SCVX');



#### PENS #### PENS #### PENS #### PENS #### PENS #### PENS #### PENS #### PENS #####
#### PENS #### PENS #### PENS #### PENS #### PENS #### PENS #### PENS #### PENS #####
#### PENS #### PENS #### PENS #### PENS #### PENS #### PENS #### PENS #### PENS #####

pAxes = {'color':[0,0,0], 'alpha':0.2, 'linewidth':1.0, 'linestyle':'-','zorder':1} #,'scale':3}
pNegAaes = {'scale':1.0, 'normalize':True, 'color':[0,0,0], 'alpha':.1, 'linewidth':1.0, 'linestyle':'-'};
pTicks = {'width':0.1,'interval':0.2,'directs':'cyclic'};
pCubeHull = {'facecolor':[0,0,0], 'edgecolor':[0,0,0],'linewidth':1.0, 'edgesolid':False,'alpha':0.02, 'edge_alpha':0.0, 'zorder':10}
pCubeHull2 = {'facecolor':[0,0,0], 'edgecolor':[0,0,0],'linewidth':1.0, 'edgesolid':False, 'edge_alpha':0.0, 'zorder':10}
pCubeEdges = {'alpha':.2, 'color':[0,0,0],'linestyle':'--','linewidth':.5};
pSimplex = {'linewidth':1.,'linestyle':'-','color':[0,0,0],'alpha':0.8,'zorder':10}    
pArrow = {'color':[0,0,0], 'alpha':0.5,'linewidth':1,'linestyle':'-','width':0.05,'head_width':0.15,'overhang':0.05,'shape':'full','zorder':10}
pArrowShadow = {'color':[0,0,0], 'alpha':0.05,'linewidth':1,'linestyle':'-','width':0.05,'head_width':0.15,'overhang':0.05,'shape':'full','zorder':10}
pCircle = {'facecolor':[0,0,1],'edgecolor':[0,0,0],'linewidth':1.0,'alpha':0.2,'zorder':10}
pSphere = {'facecolor':[1,1,1],'edgecolor':[0,0,0],'linestyle':'-','linewidth':1,'alpha':0.0,'edge_alpha':1.0,'fill':True,'edgesolid':True,'zorder':10}
pSphere2 = {'facecolor':[1,1,1],'edgecolor':[0,0,0],'linestyle':'-','linewidth':1,'alpha':0.0,'edge_alpha':0.02,'fill':False,'edgesolid':True,'zorder':10}
pSpherePairs = {'linestyle':'-','linewidth':1,'color':[0,0,0],'alpha':0.5};
pHull = {'facecolor':[0,0,0], 'edgecolor':[0,0,0],'linewidth':1.0, 'edgesolid':True,'alpha':0.02, 'edge_alpha':0.1, 'zorder':10}
pLine = {'alpha':0.1, 'color':[0,0,0],'linestyle':'-','linewidth':2.};

paramsAxes = pAxes.copy()
paramsNegAxes = pNegAxes.copy()
paramsTicks = pTicks.copy()
paramsCubeHull = pCubeHull.copy()
paramsCubeHull2 = pCubeHull2.copy()
paramsCubeEdges = pCubeEdges.copy()
paramsSimplex = pSimplex.copy()
paramsArrow = pArrow.copy()
paramsArrowShadow = pArrowShadow.copy()
paramsCircle = pCircle.copy()
paramsSphere = pSphere.copy()
paramsSphere2 = pSphere2.copy()
paramsSpherePairs = pSpherePairs.copy()
paramsHull = pHull.copy();
paramsLine = pLine.copy();

paramsSphere['facecolor'] = [0,0,1]; paramsSphere['alpha'] = 0.1;
paramsSphere['edgecolor'] = [0,0,1]; paramsSphere['edge_alpha'] = 0.2;
paramsSphere['linestyle'] = '--'; paramsSphere['linewidth'] = 1.;
paramsSpherePairs['color'] = [0,0,1]; paramsSpherePairs['alpha'] = 0.2; 
paramsSpherePairs['linewidth'] = 0.5; paramsSpherePairs['linestyle'] = '--';
paramsTicks['width']=0.05; paramsTicks['interval']=0.2;
paramsTicks['directs'] = 'cyclic'; #paramsTicks['direction'] = [2,0,0]; 
paramsArrow['color']=[0,0,0]; paramsArrow['alpha']=0.4;
# paramsArrow['linewidth']=1.; paramsArrow['linestyle']='-';
paramsArrow['width']=0.02; paramsArrow['head_width']=0.08;


frgbune = [0.5,0.,1.];
frgbuna = [0.7,0.7,0.7];
frgbuso = [1.,0.3,1.];
frgbsamp = [0.,0.,0.];
frgbu1 = [0,0,1];
frgbu2 = [1,0,0];
frgbi = [0,0,0];
frgbdes = [0,0,0];


PENS = {};

PENS['sph1'] = PEN({'frgba':[0,0,1,0.0],'ergba':[0,0,1,0.2],'ls':'-','lw':3})
PENS['sph2'] = PEN({'frgba':[1,0,0,0.0],'ergba':[1,0,0,0.2],'ls':'-','lw':3})

PENS['opt1a'] = PEN({'frgba':[0,0,1,0.1],'ergba':[0,0,1,0.2],'ls':'-','lw':6})
PENS['opt1b'] = PEN({'frgba':[0,0,1,0.1],'ergba':[0,0,1,0.2],'ls':'-','lw':2})
PENS['opt2a'] = PEN({'frgba':[1,0,0,0.1],'ergba':[1,0,0,0.2],'ls':'-','lw':6})
PENS['opt2b'] = PEN({'frgba':[1,0,0,0.1],'ergba':[1,0,0,0.2],'ls':'-','lw':2})



PENS['axes'] = PEN({'frgba':[0,0,0,0],'ergba':[0,0,0,0.2],'ls':'-','lw':1,'zord':1})
PENS['pairs1'] = PEN({'frgba':[0,0,0,0.05],'ergba':[0,0,1,0.2],'ls':'-','lw':1.5})
PENS['pairs2'] = PEN({'frgba':[0,0,0,0.05],'ergba':[1,0,0,0.2],'ls':'-','lw':1.5})
# PENS['tphcurvs'] = PEN({'frgba':[0,0,1,0],'ergba':[0,0,1,0.02],'ls':'-','lw':2}) #### original in papers
PENS['tphcurvs'] = PEN({'frgba':[0,0,1,0],'ergba':[0,0,1,0.02],'ls':'-','lw':4})

PENS['u1'] = PEN({'frgba':frgbu1 + [1],'ergba':[0,0,0,0.2],'ls':'-','lw':1.,'zord':100})
PENS['u2'] = PEN({'frgba':frgbu2 + [1],'ergba':[0,0,0,0.2],'ls':'-','lw':1.,'zord':100})
PENS['une'] = PEN({'frgba':frgbune + [1],'ergba':[0,0,0,0.2],'ls':'-','lw':1.,'zord':100})
PENS['una'] = PEN({'frgba':frgbuna + [1],'ergba':[0,0,0,0.2],'ls':'-','lw':1.,'zord':100})
PENS['uso'] = PEN({'frgba':frgbuso +[1],'ergba':[0,0,0,0.2],'ls':'-','lw':1.,'zord':100})
PENS['usamp'] = PEN({'frgba':frgbsamp +[1],'ergba':[0,0,0,0.2],'ls':'-','lw':1.,'zord':200})
PENS['high1'] = PEN({'frgba':[0.,0.,0.,1],'ergba':[0,0,0,1.],'ls':'-','lw':4.,'zord':10})
PENS['highexp1'] = PEN({'frgba':[0.,0.,0.,1],'ergba':[0,0,0]+[0.3],'ls':'-','lw':2.,'zord':10})
PENS['highexp2'] = PEN({'frgba':[0.,0.,0.,1],'ergba':[0,0,0]+[0.3],'ls':'-','lw':2.,'zord':10})

PENS['Pune'] = PEN({'frgba':frgbune + [0.02],'ergba':frgbune + [0.5],'ls':'-','lw':1.,'zord':10})
PENS['Puna'] = PEN({'frgba':frgbuna + [0.02],'ergba':frgbuna + [0.8],'ls':'-','lw':1.,'zord':10})
PENS['Pu1'] = PEN({'frgba':frgbu1 + [0.02],'ergba':frgbu1 + [0.5],'ls':'-','lw':1.,'zord':10})
PENS['Pu2'] = PEN({'frgba':frgbu2 + [0.02],'ergba':frgbu2 + [0.5],'ls':'-','lw':1.,'zord':10})
PENS['Pinter'] = PEN({'frgba':frgbi + [0.0],'ergba':frgbi + [0.5],'ls':'-','lw':3.,'zord':10})

PENS['fPune'] = PEN({'frgba':frgbune + [0.1],'ergba':frgbune + [0.],'ls':'-','lw':1.,'zord':10})
PENS['fPuna'] = PEN({'frgba':frgbuna + [0.1],'ergba':frgbuna + [0.],'ls':'-','lw':1.,'zord':1})
PENS['fPu1'] = PEN({'frgba':frgbu1 + [0.1],'ergba':frgbu1 + [0.],'ls':'-','lw':1.,'zord':1})
PENS['fPu2'] = PEN({'frgba':frgbu2 + [0.1],'ergba':frgbu2 + [0.],'ls':'-','lw':1.,'zord':1})


# SEE BELOW: PENS['xne'] = PEN({'frgba':frgbune + [0.1],'ergba':frgbune + [1.],'ls':'-','lw':5.,'zord':10})
PENS['xna'] = PEN({'frgba':frgbuna + [0.1],'ergba':frgbuna + [1.],'ls':'-','lw':5.,'zord':10})
PENS['xopt1'] = PEN({'frgba':frgbu1 + [0.1],'ergba':frgbu1 + [1.],'ls':'-','lw':5.,'zord':10})
PENS['xopt2'] = PEN({'frgba':frgbu2 + [0.1],'ergba':frgbu2 + [1.],'ls':'-','lw':5.,'zord':10})
PENS['xdes'] = PEN({'frgba':frgbdes + [0.1],'ergba':frgbdes + [1.],'ls':'--','lw':2.,'zord':10})
PENS['xso'] = PEN({'frgba':frgbuso + [0.1],'ergba':frgbuso + [1.],'ls':'-','lw':2.,'zord':10})


# PENS['usamp'] = PEN({'frgba': [0,0,0,1.],'ergba': [0,0,0,1.],'ls':'-','lw':0.3,'zord':10})
PENS['tsamp'] = PEN({'frgba': [0,0,0,1.],'ergba': [0,0,0,1.],'ls':'-','lw':0.3,'zord':10})
PENS['tblow'] = PEN({'frgba': [0,0,0,1.],'ergba': [0,0,0,1.],'ls':'-','lw':0.3,'zord':10,'msty':'','msz':4})
PENS['tblow1'] = PEN({'frgba': [0,0,0,1.],'ergba': [0,0,0,1.],'ls':'-','lw':0.3,'zord':10,'msty':'x','msz':6})
PENS['tblow2'] = PEN({'frgba': [0,0,0,1.],'ergba': [0,0,0,1],'ls':'-','lw':0.3,'zord':10,'msty':'','msz':4})
PENS['tblow3'] = PEN({'frgba': [0,0,0,1.],'ergba': [0,0,0,1],'ls':'-','lw':0.3,'zord':10,'msty':'','msz':4})
PENS['xblowups'] = PEN({'frgba': [0,0,0,1.],'ergba': [0,0,0,1.],'ls':'-','lw':2,'zord':10})
PENS['xblowups2'] = PEN({'frgba': [0,0,0,1.],'ergba': [0,0,1,1.],'ls':'-','lw':2,'zord':10})
PENS['ablowups'] = PEN({'frgba': [0,0,0,1.],'ergba': [0,0,0,1.],'ls':'-','elw':0.5,'zord':10})
PENS['xbase'] = PEN({'frgba':frgbune + [0.1],'ergba':frgbune + [1.],'ls':'-','lw':5.,'zord':10})
PENS['wblow'] = PEN({'frgba': [0,0,0,1.],'ergba': [0,0,0,1.],'ls':'-','lw':0.3,'zord':10})
PENS['cutoff'] = PEN({'frgba': [0,0,0,1.],'ergba': [0,0,0,1.],'ls':'-','lw':0.3,'zord':10,'msty':'','msz':4})




PENS['xdes'] = PEN({'frgba':[0,0,0,1],'ergba':[0,0,0,0.4],'ls':'--','lw':2.,'zord':10})
# PENS['xne'] = PEN({'frgba':frgbune + [0.1],'ergba':frgbune + [0.3],'ls':'--','lw':4.,'zord':10})

# PENS['basis'] = PEN({'frgba':[0,0,0,1],'ergba':[0,0,0,0.5],'ls':'-','lw':4.,'zord':11})
PENS['basis_scaled'] = PEN({'frgba':[0,0,0,1],'ergba':[0,0,1,0.6],'ls':'-','lw':3.,'zord':10})
PENS['approx'] = PEN({'frgba':[0,0,1,0.5],'ergba':[0,0,1,0.5],'ls':'-','lw':3.,'zord':10})
PENS['exact'] = PEN({'frgba':[0,0,1,1],'ergba':frgbune + [0.7],'ls':'-','lw':4.,'zord':10})

#blig
p1alpha = 1.0;
p2alpha = 1.0;

PENS['xdes1'] = PEN({'frgba':[0,0,0,1],'ergba':[0,0,0,0.4*p1alpha],'ls':'--','lw':2.,'zord':10})
PENS['xdes2'] = PEN({'frgba':[0,0,0,1],'ergba':[0,0,0,0.4*p2alpha],'ls':'--','lw':2.,'zord':10})
nashcolor = [0.5,0.5,0.5];
exactcolor = frgbune
# nashcolor = frgbune;
PENS['xne'] = PEN({'frgba':frgbune + [0.1],'ergba':nashcolor + [0.3],'ls':'-','lw':4.,'zord':10})
PENS['xne1'] = PEN({'frgba':frgbune + [0.1],'ergba':nashcolor + [p1alpha],'ls':'-','lw':8.,'zord':10})
PENS['xne2'] = PEN({'frgba':frgbune + [0.1],'ergba':nashcolor + [p2alpha],'ls':'-','lw':8.,'zord':10})
PENS['basis1'] = PEN({'frgba':[0,0,0,1],'ergba':[0,0,1,p1alpha],'ls':'-','lw':16.,'zord':11}); #None , has priority over basis 
PENS['basis2'] = PEN({'frgba':[0,0,0,1],'ergba':[1,0,0,p2alpha],'ls':'-','lw':16.,'zord':11})  #None , has priority over basis
PENS['exact1'] = PEN({'frgba':[0,0,0,1],'ergba':exactcolor + [p1alpha],'ls':'-','lw':8.,'zord':11}); #None , has priority over basis 
PENS['exact2'] = PEN({'frgba':[0,0,0,1],'ergba':exactcolor + [p2alpha],'ls':'-','lw':8.,'zord':11})  #None , has priority over basis
PENS['approx1'] = PEN({'frgba':[0,0,0,1],'ergba':[0,0,1,p1alpha],'ls':'-','lw':8.,'zord':11}); #None , has priority over basis 
PENS['approx2'] = PEN({'frgba':[0,0,0,1],'ergba':[1,0,0,p2alpha],'ls':'-','lw':8.,'zord':11})  #None , has priority over basis


PENS['basis'] = PEN({'frgba':[0,0,0,1],'ergba':[0,0,0,1.],'ls':'-','lw':8.,'zord':11})