import numpy as np
import numpy.linalg as mat
import scipy as sp
import scipy.linalg as smat
import copy
import random
import itertools

from matplotlib.colors import LinearSegmentedColormap
# import cvxpy as cp

# import osmnx as ox
# import networkx as nx
# import geopandas as gpd
# import pandas as pd
# import peartree as pt #turns GTFS feed into a graph
# import folium
# import gtfs_functions as gtfs

# import matplotlib.pyplot as plt
# from shapely.geometry import Polygon, Point

# from matplotlib.patches import FancyArrow
# from itertools import product 
# from random import sample
# from shapely.geometry import Polygon, Point

# import sklearn as sk
# from sklearn import cluster as cluster
# from sklearn.cluster import OPTICS, cluster_optics_dbscan
# from sklearn import metrics
# from sklearn.cluster import DBSCAN

# from scipy.spatial import ConvexHull, convex_hull_plot_2d

# import time
import warnings
warnings.filterwarnings('ignore')

# import sys
# sys.path.append('/Users/dan/Documents/code');
# from drawing import * 
# from mathz import * 
# sys.path.remove('/Users/dan/Documents/code');




def curveColor2DGrid1(tag,ins={}):

    colorvals = {(0,1): np.array([0,0,1,1]),(1,1): np.array([0,1,0,1]),
              (0,0): np.array([0,0,0,1]),(1,0): np.array([1,0,0,1]),
              'mid1': np.array([1,0,1,1]),'mid2': np.array([1,0,1,1])}    

    version = 'colorgrid';
    divtype = 'polar'
    direction = 'tan'
    if 'version' in ins: version = ins['version']
    if 'divtype' in ins: divtype = ins['divtype'];
    if 'colors' in ins: colorvals = ins['colors']
    if 'direction' in ins: direction = ins['direction']

    if 'xdivs' in ins: th1divs = ins['xdivs'];
    if 'ydivs' in ins: th2divs = ins['ydivs'];
    if 'tdivs' in ins: tdivs = ins['tdivs'];

    alpha_manual = 1; use_alpha_manual = False;
    if 'alpha_manual' in ins: alpha_manual = ins['alpha_manual']; use_alpha_manual = True;

    ph = tag[0];
    ts = [];
    if 'ts' in ins: ts = ins['ts'];
    if len(ts)==0: ts = [tag[1]];
    if len(ts)==1: divtype = 'single';

    ths = np.array([[np.arctan2(tt*np.cos(ph),1),np.arctan2(tt*np.sin(ph),1)] for tt in ts])
    minth1 = np.min(ths[:,0]); minth2 = np.min(ths[:,1]);
    maxth1 = np.max(ths[:,0]); maxth2 = np.max(ths[:,1]);

    # mask1a = minth1 <= th1divs; mask1b = th1divs <= maxth1;
    # mask2a = minth2 <= th2divs; mask2b = th2divs <= maxth2;
    # th1divs = th1divs[mask1a & mask1b];
    # th2divs = th2divs[mask2a & mask2b];

    if divtype == 'square':
        inds = [];
        colors = [];
        tsegs = [];

        for i,th1 in enumerate(th1divs):
            for j,th2 in enumerate(th2divs):
                if (i<len(th1divs)-1) and (j<len(th2divs)-1):
                    th1int = [th1,th1divs[i+1]];
                    th2int = [th2,th2divs[j+1]];

                    mask1a = th1int[0] <= ths[:,0]; mask1b = ths[:,0] <= th1int[1]
                    mask2a = th2int[0] <= ths[:,1]; mask2b = ths[:,1] <= th2int[1]
                    mask = (mask1a & mask1b) & (mask2a & mask2b);
                    temp = list(np.where(mask)[0])

                    if len(temp)>0:
                        if temp[0]  > 0: temp = [temp[0]-1] + temp;
                        if temp[-1] < len(ts)-1: temp = temp + [temp[-1]+1]; 
                    # thint = ths[mask];                
                        inds.append(temp.copy())
                        tsegs.append(ts[temp])
                        scale = (1./(np.pi/2))

                        midth1 = np.mean(th1int)
                        midth2 = np.mean(th2int)


                        if version == 'polargrid':                    
                            ph_temp = np.arctan2(np.tan(midth2),np.tan(midth1));
                            if ph_temp <= np.pi/4: t_temp = np.tan(midth1)/np.cos(ph_temp);
                            else: t_temp = np.tan(midth2)/np.sin(ph_temp);
                            params = {}
                            params['colors'] = colorvals;
                            params['direction'] = direction
                            color = polarColorGrid(ph_temp,t_temp,params=params)

                            if use_alpha_manual: color[3] = alpha_manual; 

                        # elif version == 'colorgrid': color = colorGrid(midth1*scale,midth2*scale);
                        else:
                            params = {}
                            params['colors'] = colorvals;
                            params['direction'] = direction
                            color = colorGrid(midth1*scale,midth2*scale,params=params);

                            if use_alpha_manual: color[3] = alpha_manual; 

                        colors.append(color);

    if divtype == 'polar':
        inds = [];
        colors = [];
        tsegs = [];

        for i,t in enumerate(tdivs):
            if (i<len(tdivs)-1):
                tint = [t,tdivs[i+1]];
                mask1a = tint[0] <= ts; mask1b = ts <= tint[1]
                mask = mask1a & mask1b
                temp = list(np.where(mask)[0])


                if len(temp)>0:
                    if temp[0]  > 0: temp = [temp[0]-1] + temp;
                    if temp[-1] < len(ts)-1: temp = temp + [temp[-1]+1]; 
                    # thint = ths[mask];                
                    inds.append(temp.copy())
                    tsegs.append(ts[temp])
                    scale = (1./(np.pi/2))
                    midt = np.mean(tint);
                    th = np.array([np.arctan2(midt*np.cos(ph),1),np.arctan2(midt*np.sin(ph),1)])
                    midth1 = th[0]; midth2 = th[1];
                    if version == 'polargrid':              
                        params = {}
                        params['colors'] = colorvals;
                        params['direction'] = direction
                        color = polarColorGrid(ph,t,params=params)
                        if use_alpha_manual: color[3] = alpha_manual; 

                    else:
                        params = {}
                        params['colors'] = colorvals;
                        params['direction'] = direction
                        color = colorGrid(midth1*scale,midth2*scale,params=params);
                        if use_alpha_manual: color[3] = alpha_manual; 
                    colors.append(color);

    if divtype == 'theta':
        inds = [];
        colors = [];
        tsegs = [];

        for i,t in enumerate(tdivs):
            if (i<len(tdivs)-1):
                tint = [t,tdivs[i+1]];
                mask1a = tint[0] <= ts; mask1b = ts <= tint[1]
                mask = mask1a & mask1b
                temp = list(np.where(mask)[0])


                if len(temp)>0:
                    if temp[0]  > 0: temp = [temp[0]-1] + temp;
                    if temp[-1] < len(ts)-1: temp = temp + [temp[-1]+1]; 
                    # thint = ths[mask];                
                    inds.append(temp.copy())
                    tsegs.append(ts[temp])

                    if False: scale = (1./(np.pi/2))
                    else: scale = 1;

                    midt = np.mean(tint);
                    if direction == 'theta1': th = np.array([ph,np.arctan2(midt,1)])
                    if direction == 'theta2': th = np.array([np.arctan2(midt,1),ph])
                    
                    midth1 = th[0]; midth2 = th[1];
                    if False: #version == 'polargrid':              
                        params = {}
                        params['colors'] = colorvals;
                        params['direction'] = direction
                        color = polarColorGrid(ph,t,params=params)
                        if use_alpha_manual: color[3] = alpha_manual; 
                    else:
                        params = {}
                        params['colors'] = colorvals;
                        params['direction'] = direction
                        color = colorGrid(midth1*scale,midth2*scale,params=params);
                        if use_alpha_manual: color[3] = alpha_manual; 
                    colors.append(color);


    if divtype == 'single':
        inds = []; tsegs = []; colors = [];
        t = ts[0];
        if version == 'polargrid':
            params = {}
            params['colors'] = colorvals;
            params['direction'] = direction
            color = polarColorGrid(ph,t,params=params)
            if use_alpha_manual: color[3] = alpha_manual; 
        else:
            params = {}
            params['colors'] = colorvals;
            params['direction'] = direction

            if direction == 'tan': th = np.array([np.arctan2(t*np.cos(ph),1),np.arctan2(t*np.sin(ph),1)])
            if direction == 'cot': th = np.array([np.arctan2(t*np.cos(ph),1),np.arctan2(1,t*np.sin(ph))])
            if direction == 'theta1': th = np.array([ph,np.arctan2(t,1)])
            if direction == 'theta2': th = np.array([np.arctan2(t,1),ph])
            th1 = th[0]; th2 = th[1];

            ph0,t0 = np.ar
            color = colorGrid(th1,th2,params=params);
            if use_alpha_manual: color[3] = alpha_manual; 
        colors = color; 

    return inds,tsegs,colors;



def curveColorByCoeffs(tag,ins={}):

    ph = tag[0]; ts = [];
    if 'ts' in ins: ts = ins['ts'];
    if len(ts)==0: ts = [tag[1]];


    colorind = 0;
    if 'colorind' in ins: colorind = ins['colorind'];

    # coeffgridparams[tag]['data'] = GAME1.EXP[tag];
    #coeffgridparams[tag]['data'] = GAME1.EXP[tag]['coeffs_color_data']
    data = ins['data'][ph];
    une = data['une']
    una = data['una']
    reigs = data['reigs']
    rVV = data['rVV']
    rVVi = data['rVVi']
    ninds = np.where(reigs < 0)[0];
    ninds = ninds[::-1];
    bind = ninds[colorind];
    lamj = reigs[bind];
    Vj = rVV[:,bind];
    Wj = rVVi[bind];

    if True:
        Wjnorm = mat.norm(Wj);
        Wj = (1./Wjnorm)*Wj
        Vj = Wjnorm*Vj
    else: Vj = (1./mat.norm(Vj))*Vj  # subtle


    basecolor = [0,0,1,1]; maxval = 1.;
    if 'basecolor' in ins: basecolor = ins['basecolor'];
    if 'maxval' in ins: maxval = ins['maxval']
    if 'tdivs' in ins: tdivs = ins['tdivs'];
    bcolor = np.array(basecolor)
    bcolor_diff = np.array([1,1,1,1]) - bcolor;

    # return inds,tsegs,colors;

    inds = [];
    colors = [];
    tsegs = [];

    for i,t in enumerate(tdivs):
        if (i<len(tdivs)-1):
            tint = [t,tdivs[i+1]];
            mask1a = tint[0] <= ts; mask1b = ts <= tint[1]
            mask = mask1a & mask1b
            temp = list(np.where(mask)[0])

            if len(temp)>0:
                if temp[0]  > 0: temp = [temp[0]-1] + temp;
                if temp[-1] < len(ts)-1: temp = temp + [temp[-1]+1]; 
                # thint = ths[mask];                
                inds.append(temp.copy())
                tsegs.append(ts[temp])
                scale = (1./(np.pi/2))
                midt = np.mean(tint);

                ### main computation 
                val = (1./(1+(lamj/midt)))*(una - une)@Vj
                ### main color choice.  
                val = np.min([maxval,np.abs(val)])
                # color = basecolor[:3] + [val/maxval];
                color = bcolor + (1.-(val/maxval))*bcolor_diff;
                colors.append(color);

    return inds,tsegs,colors;





def segmentsInSquare(dat,typ='slope'):

    if typ == 'slope':
        L = dat['L'];
        ell = dat['ell'];
        bnds = dat['bnds'];
        x1min = bnds[0]; x1max = bnds[1]; x2min = bnds[2]; x2max = bnds[3];
        ROWS = np.array([[1,0],[1,0],[0,1],[0,1]]);
        bs = np.array([x1min,x1max,x2min,x2max]);
        pts = [];
        inds = [];
        for i,row in enumerate(ROWS):
            pt = mat.inv(np.array([[L,-1],list(row)]))@np.array([-ell,bs[i]]);
            if (x1min <= pt[0] <= x1max) and (x2min <= pt[1] <= x2max):
                pts.append(pt)
                inds.append(i)

    if typ == 'pt_direct':
        start = dat['start']
        direct = dat['direct'];
        bnds = dat['bnds'];

        x1min = bnds[0]; x1max = bnds[1]; x2min = bnds[2]; x2max = bnds[3];
        ROWS = np.array([[1,0],[1,0],[0,1],[0,1]]);
        bs = np.array([x1min,x1max,x2min,x2max]);
        pts = [];
        inds = [];
        for i,row in enumerate(ROWS):
            z = (bs[i]-row@start)/(row@direct)
            pt = direct*z + start
            if (x1min <= pt[0] <= x1max) and (x2min <= pt[1] <= x2max):
                pts.append(pt)
                inds.append(i)    
    return pts,inds


def filter_curves(curves,bnd=10,gap=0.2):
    newcurves = []
    for i,curve in enumerate(curves):
        newcurves.append([])
        newcurve = [];
        for j,pt in enumerate(curve):
            if mat.norm(pt)<bnd:        
                if len(newcurve)==0: newcurve.append(pt);
                else:
                    oldpt = newcurve[-1];
                    if mat.norm(pt-oldpt)<gap: newcurve.append(pt);
                    else: newcurves[i].append(np.array(newcurve)); newcurve = [pt]
            else: pass
        if len(newcurve)>0:
            newcurves[i].append(np.array(newcurve))
    return newcurves
    

########

pAxes = {'color':[0,0,0], 'alpha':0.2, 'linewidth':1.0, 'linestyle':'-','head_width':.05,'zorder':1} #,'scale':3}
pNegAxes = {'scale':1.0, 'normalize':True, 'color':[0,0,0], 'alpha':.1, 'linewidth':1.0, 'linestyle':'-'};
pTicks = {'width':0.1,'interval':0.2,'directs':'cyclic'};
pHull = {'facecolor':[0,0,0], 'edgecolor':[0,0,0],'linewidth':1.0, 'edgesolid':False,'alpha':0.02, 'edge_alpha':0.0, 'zorder':10}
pCubeHull = {'facecolor':[0,0,0], 'edgecolor':[0,0,0],'linewidth':1.0, 'edgesolid':True,'alpha':0.04, 'edge_alpha':0.6, 'zorder':10}
pCubeHull1 = {'facecolor':[0,0,0], 'edgecolor':[0,0,0],'linewidth':1.0, 'edgesolid':True,'alpha':0.02, 'edge_alpha':0.02, 'zorder':10}
pCubeHull2 = {'facecolor':[0,0,0], 'edgecolor':[0,0,0],'linewidth':1.0, 'edgesolid':True,'alpha':0.04, 'edge_alpha':0.6, 'zorder':10}
# pCubeHull2 = {'facecolor':[0,0,0], 'edgecolor':[0,0,0],'linewidth':1.0, 'edgesolid':False, 'edge_alpha':0.0, 'zorder':10}
pCubeEdges = {'alpha':.2, 'color':[0,0,0],'linestyle':'--','linewidth':.5};
pSimplex = {'linewidth':1.,'linestyle':'-','color':[0,0,0],'alpha':0.8,'zorder':10}    
pArrow = {'color':[0,0,0], 'alpha':0.5,'linewidth':1,'linestyle':'-','width':0.05,'head_width':0.15,'overhang':0.05,'shape':'full','zorder':10}
pArrowShadow = {'color':[0,0,0], 'alpha':0.05,'linewidth':1,'linestyle':'-','width':0.05,'head_width':0.15,'overhang':0.05,'shape':'full','zorder':10}
pCircle = {'facecolor':[0,0,1],'edgecolor':[0,0,0],'linewidth':1.0,'alpha':0.2,'zorder':10}
pSphere = {'facecolor':[1,1,1],'edgecolor':[0,0,0],'linestyle':'-','linewidth':1,'alpha':0.0,'edge_alpha':1.0,'fill':True,'edgesolid':True,'zorder':10}
pSpherePairs = {'linestyle':'-','linewidth':1,'color':[0,0,0],'alpha':0.5};
pPatch = {'facecolor':[0,0,0], 'edgecolor':[0,0,0],'linewidth':1.0, 'edgesolid':True,'alpha':0.04, 'edge_alpha':0.6, 'zorder':10}
pLine = {'color':[0,0,0],'linewidth':1.0,'linestyle':'-','alpha':0.04,  'zorder':10}


paramsAxes = pAxes.copy()
paramsNegAxes = pNegAxes.copy()
paramsTicks = pTicks.copy()
paramsHull = pHull.copy();
paramsCubeHull = pCubeHull.copy()
paramsCubeHull1 = pCubeHull1.copy()
paramsCubeHull2 = pCubeHull2.copy()
paramsCubeEdges = pCubeEdges.copy()
paramsSimplex = pSimplex.copy()
paramsArrow = pArrow.copy()
paramsArrowShadow = pArrowShadow.copy()
paramsCircle = pCircle.copy()
paramsSphere = pSphere.copy()
paramsSpherePairs = pSpherePairs.copy()
paramsPatch = pPatch.copy();
paramsLine = pLine.copy();


DEFAULT_AXES = np.array([[1.,0.],
                         [0.,1.],
                         [-0.7,-0.7],
                         [0.7,-0.7],
                         [-0.9,0.3],
                         [0.9,0.3],
                         [-1.,0.],
                         [0.,-1.],
                         [-0.3,0.9],
                         [0.3,0.9]])

class PEN:
    def __init__(self,ins={}):
        ### ROOTS
        self.faceRGBA_val = [0,0,0,0.1]; 
        self.edgeRGBA_val = [0,0,0,0.4];
        self.edgelinewidth_val = 1.
        self.facelinewidth_val = 1.
        self.markerstyle_val = ''
        self.markersize_val = 1.
        self.linestyle_val = '-'

        self.scale_val = 1.;
        self.zorder_val = 1.;
        self.normalize_val = False;
        self.interval_val = 1.;
        self.directions_val = list(np.zeros(16));
        self.directions_val[0] = 1;

        self.width_val = 0.1;
        self.headwidth_val = 0.15;
        self.headlength_val = 0.2;
        self.headversion_val = '2d';
        self.headdirection_val = None;

        self.tickwidth_val = 0.15;
        self.ticktype_val = 'b';
        self.tickdirection_val = None;

        self.overhang_val = 0.1;
        self.shape_val = 'arrow';
        self.edgesolid_val = True;
        self.fill_val = True;

        self.consistify();
        self.update(ins=ins);

    def consistify(self):
        # MAYBE OVERKILL 
        self.RGBA = self.faceRGBA_val;
        self.rgba = self.faceRGBA_val;
        self.faceRGBA = self.faceRGBA_val;
        self.edgeRGBA = self.edgeRGBA_val;
        self.frgba = self.faceRGBA_val;
        self.ergba = self.edgeRGBA_val;

        self.color = self.faceRGBA_val[:3];
        self.facecolor = self.faceRGBA_val[:3];
        self.edgecolor = self.edgeRGBA_val[:3];
        self.fc = self.faceRGBA_val[:3];
        self.ec = self.edgeRGBA_val[:3];

        self.alpha = self.faceRGBA_val[3];
        self.facealpha = self.faceRGBA_val[3];
        self.edgealpha = self.edgeRGBA_val[3];
        self.fa = self.faceRGBA_val[3];
        self.ea = self.edgeRGBA_val[3];

        self.linewidth = self.edgelinewidth_val;
        self.edgewidth = self.edgelinewidth_val;
        self.edgelinewidth = self.edgelinewidth_val;
        self.facelinewidth = self.facelinewidth_val;

        self.lw = self.edgelinewidth_val;
        self.fw = self.facelinewidth_val;
        self.ew = self.edgelinewidth_val;
        self.flw = self.facelinewidth_val;
        self.elw = self.edgelinewidth_val;

        self.markerstyle = self.markerstyle_val;
        self.msty = self.markerstyle_val;

        self.markersize = self.markersize_val;
        self.msz = self.markersize_val

        self.linestyle = self.linestyle_val;
        self.ls = self.linestyle_val;


        ############################################################
        ############################################################
        ############################################################
        self.scale = self.scale_val;
        self.zorder = self.zorder_val;
        self.zord = self.zorder_val;
        self.normalize = self.normalize_val;
        self.interval = self.interval_val;
        self.directs = self.directions_val; 
        self.directions = self.directions_val;

        self.width = self.width_val;
        self.wid = self.width_val;
        self.headwidth = self.headwidth_val;
        self.hwid = self.headwidth_val;
        self.headversion = self.headversion_val;
        self.hver = self.headversion_val;
        self.headlength = self.headlength_val;
        self.hlen = self.headlength_val;
        self.headdirection = self.headdirection_val;
        self.hdir = self.headdirection_val;

        self.ttype = self.ticktype_val;
        self.tdir = self.tickdirection_val;
        self.twid = self.tickwidth_val;

        self.overhand = self.overhang_val;
        self.shape = self.shape_val;
        self.edgesolid = self.edgesolid_val;
        self.fill = self.fill_val;


    def update(self,ins={}):

        ## ORDER INDICATES PRIORITY

        if 'color' in ins:
            temp = ins['color'];
            self.faceRGBA_val = [temp[0],temp[1],temp[2],self.faceRGBA_val[3]]; 
            self.edgeRGBA_val = [temp[0],temp[1],temp[2],self.edgeRGBA_val[3]];
        if 'c' in ins:
            temp = ins['c'];
            self.faceRGBA_val = [temp[0],temp[1],temp[2],self.faceRGBA_val[3]]; 
            self.edgeRGBA_val = [temp[0],temp[1],temp[2],self.edgeRGBA_val[3]];
        if 'facecolor' in ins: temp = ins['facecolor']; self.faceRGBA_val = [temp[0],temp[1],temp[2],self.faceRGBA_val[3]]; 
        if 'edgecolor' in ins: temp = ins['edgecolor']; self.edgeRGBA_val = [temp[0],temp[1],temp[2],self.edgeRGBA_val[3]];
        if 'fc' in ins: temp = ins['fc']; self.faceRGBA_val = [temp[0],temp[1],temp[2],self.faceRGBA_val[3]]; 
        if 'ec' in ins: temp = ins['ec']; self.edgeRGBA_val = [temp[0],temp[1],temp[2],self.edgeRGBA_val[3]];         


        if 'alpha' in ins:
            self.faceRGBA_val[3] = ins['alpha'];
            self.edgeRGBA_val[3] = ins['alpha'];
        if 'facealpha' in ins: self.faceRGBA_val[3] = ins['facealpha'];
        if 'edgealpha' in ins: self.edgeRGBA_val[3] = ins['edgealpha'];
        if 'fa' in ins: self.faceRGBA_val[3] = ins['fa'];
        if 'ea' in ins: self.edgeRGBA_val[3] = ins['ea'];

        if 'RGBA' in ins:
            self.faceRGBA_val = ins['RGBA'];
            self.edgeRGBA_val = ins['RGBA'];
        if 'rgba' in ins:
            self.faceRGBA_val = ins['rgba'];
            self.edgeRGBA_val = ins['rgba'];

        if 'faceRGBA' in ins: self.faceRGBA_val = ins['faceRGBA'];
        if 'edgeRGBA' in ins: self.edgeRGBA_val = ins['faceRGBA'];
        if 'fRGBA' in ins: self.faceRGBA_val = ins['fRGBA'];
        if 'eRGBA' in ins: self.edgeRGBA_val = ins['eRGBA'];
        if 'frgba' in ins: self.faceRGBA_val = ins['frgba'];
        if 'ergba' in ins: self.edgeRGBA_val = ins['ergba'];


        ####################################################################################

        if 'linewidth' in ins:
            self.facelinewidth_val = ins['linewidth'];
            self.edgelinewidth_val = ins['linewidth'];
        if 'lw' in ins:
            self.facelinewidth_val = ins['lw'];
            self.edgelinewidth_val = ins['lw'];
        if 'facelinewidth' in ins: self.facelinewidth_val = ins['facelinewidth'];
        if 'facewidth' in ins: self.facelinewidth_val = ins['facewidth'];
        if 'edgelinewidth' in ins: self.edgelinewidth_val = ins['edgelinewidth'];
        if 'edgewidth' in ins: self.edgelinewidth_val = ins['edgewidth'];

        if 'flw' in ins: self.facelinewidth_val = ins['flw'];
        if 'elw' in ins: self.edgelinewidth_val = ins['elw'];
        if 'fw' in ins: self.facelinewidth_val = ins['fw'];
        if 'ew' in ins: self.edgelinewidth_val = ins['ew'];


        if 'linestyle' in ins: self.linestyle_val = ins['linestyle'];
        if 'ls' in ins: self.linestyle_val = ins['ls'];

        if 'markerstyle' in ins: self.markerstyle_val = ins['markerstyle'];
        if 'msty' in ins: self.markerstyle_val = ins['msty'];

        if 'markersize' in ins: self.markersize_val = ins['markersize']
        if 'msz' in ins: self.markersize_val = ins['msz']


        ####################################################################################
        ####################################################################################
        ####################################################################################
        ####################################################################################
        
        if 'scale' in ins: self.scale_val = ins['scale']
        if 'zorder' in ins: self.zorder_val = ins['zorder'];
        if 'zord' in ins: self.zorder_val = ins['zord'];
        if 'normalize' in ins: self.normalize_val = ins['normalize'];
        if 'interval' in ins: self.interval_val = ins['interval'];
        if 'directs' in ins: self.directions_val = ins['directs']
        if 'directions' in ins: self.directions_val = ins['directions'];

        if 'width' in ins: self.width_val = ins['width'];
        if 'wid' in ins: self.width_val = ins['wid'];
        if 'headwidth' in ins: self.headwidth_val = ins['headwidth'];
        if 'hwid' in ins: self.headwidth_val = ins['hwid'];
        if 'headversion' in ins: self.headversion_val = ins['headversion'];
        if 'hver' in ins: self.headversion_val = ins['hver']
        if 'hlen' in ins: self.headlength_val = ins['hlen']
        if 'headlength' in ins: self.headlength_val = ins['headlength'];  
        if 'headdirection' in ins: self.headdirection_val = ins['headdirection'];
        if 'hdir' in ins: self.headdirection_val = ins['hdir'];

        if 'ttype' in ins: self.ticktype_val = ins['ttype'];
        if 'tdir' in ins: self.tickdirection_val = ins['tdir'];
        if 'twid' in ins: self.tickwidth_val = ins['twid'];

        if 'overhang' in ins: self.overhang_val = ins['overhang'];
        if 'shape' in ins: self.shape = self.shape_val = ins['shape']
        if 'edgesolid' in ins: self.edgesolid_val = ins['edgesolid']
        if 'fill' in ins: self.fill_val = ins['fill'];

        self.consistify();


class LQ2GAME:
    def __init__(self,ins={}):

        self.nt = 100; self.dt = 0.01; self.t0 = 0;
        if 'nt' in ins: self.nt = ins['nt'];
        if 'dt' in ins: self.dt = ins['dt'];
        if 't0' in ins: self.t0 = ins['t0'];
        self.ts = self.dt*np.array(list(range(self.nt)))+self.t0
        if 'paramsIndAgents' in ins: self.paramsIndAgents = ins['paramsIndAgents'];


        self.nx = None; self.nu = None;
        self.nu1 = None; self.nu2 = None;
        self.nx1 = None; self.nx2 = None;
        if 'nx' in ins: self.nx = ins['nx']
        if 'nx1' in ins: self.nx1 = ins['nx1']
        if 'nx2' in ins: self.nx2 = ins['nx2']
        if 'nu' in ins: self.nu = ins['nu']
        if 'nu1' in ins: self.nu1 = ins['nu1']
        if 'nu2' in ins: self.nu2 = ins['nu2']

        if self.nu == None: self.nu = self.nu1 + self.nu2;
        if self.nx == None: self.nx = self.nx1 + self.nx2;
        if self.nx1 == None: self.nx1 = int(self.nx/2)
        if self.nx2 == None: self.nx2 = int(self.nx/2)

        self.As = np.eye(self.nx); 
        self.Bs = np.zeros([self.nx,self.nu]);
        self.Bs[:self.nu1,:self.nu1] = np.eye(self.nu1)
        self.Bs[-self.nu2:,-self.nu2:] = np.eye(self.nu2)
        if 'As' in ins: self.As = ins['As'];
        if 'Bs' in ins: self.Bs = ins['Bs'];

        if isinstance(self.As,list): self.As = np.array(self.As); 
        if isinstance(self.Bs,list): self.Bs = np.array(self.Bs);

        if len(self.As.shape)==2: self.As = np.einsum('i,jk',np.ones(self.nt),self.As)
        if len(self.Bs.shape)==2: self.Bs = np.einsum('i,jk',np.ones(self.nt),self.Bs)

        self.B1s = self.Bs[:,:,:self.nu1];
        self.B2s = self.Bs[:,:,self.nu1:];
        if 'B1s' in ins: self.B1s = ins['B1s'];
        if 'B2s' in ins: self.B2s = ins['B2s'];

        if isinstance(self.B1s,list): self.B1s = np.array(self.B1s); 
        if isinstance(self.B2s,list): self.B2s = np.array(self.B2s);

        if len(self.B1s.shape)==2: self.B1s = np.einsum('i,jk',np.ones(self.nt),self.B1s)
        if len(self.B2s.shape)==2: self.B2s = np.einsum('i,jk',np.ones(self.nt),self.B2s)

        self.Q1s = np.eye(self.nx);
        self.Q2s = np.eye(self.nx);
        self.P1s = np.zeros([self.nx,self.nx]);
        self.P2s = np.zeros([self.nx,self.nx]);

        self.R1s = np.eye(self.nu1);
        self.R2s = np.eye(self.nu2);
        self.R12s = np.zeros([self.nu2,self.nu2]);
        self.R21s = np.zeros([self.nu1,self.nu1]);

        if 'Q1s' in ins: self.Q1s = ins['Q1s']
        if 'Q2s' in ins: self.Q2s = ins['Q2s']
        if 'P1s' in ins: self.P1s = ins['P1s']
        if 'P2s' in ins: self.P2s = ins['P2s']
        if 'R1s' in ins: self.R1s = ins['R1s']
        if 'R2s' in ins: self.R2s = ins['R2s']
        if 'R12s' in ins: self.R12s = ins['R12s']
        if 'R21s' in ins: self.R21s = ins['R21s']

        if isinstance(self.Q1s,list): self.Q1s = np.array(self.Q1s); 
        if isinstance(self.Q2s,list): self.Q2s = np.array(self.Q2s); 
        if isinstance(self.P1s,list): self.P1s = np.array(self.P1s); 
        if isinstance(self.P2s,list): self.P2s = np.array(self.P2s); 
        if isinstance(self.R1s,list): self.R1s = np.array(self.R1s); 
        if isinstance(self.R2s,list): self.R2s = np.array(self.R2s); 
        if isinstance(self.R12s,list): self.R12s = np.array(self.R12s); 
        if isinstance(self.R21s,list): self.R21s = np.array(self.R21s); 

        if len(self.Q1s.shape)==2: self.Q1s = np.einsum('i,jk',np.ones(self.nt),self.Q1s)
        if len(self.Q2s.shape)==2: self.Q2s = np.einsum('i,jk',np.ones(self.nt),self.Q2s)
        if len(self.P1s.shape)==2: self.P1s = np.einsum('i,jk',np.ones(self.nt),self.P1s)
        if len(self.P2s.shape)==2: self.P2s = np.einsum('i,jk',np.ones(self.nt),self.P2s)
        if len(self.R1s.shape)==2: self.R1s = np.einsum('i,jk',np.ones(self.nt),self.R1s)
        if len(self.R2s.shape)==2: self.R2s = np.einsum('i,jk',np.ones(self.nt),self.R2s)
        if len(self.R12s.shape)==2: self.R12s = np.einsum('i,jk',np.ones(self.nt),self.R12s)
        if len(self.R21s.shape)==2: self.R21s = np.einsum('i,jk',np.ones(self.nt),self.R21s)

        self.x0 = np.ones(self.nx); 
        self.xdes = np.zeros([self.nt,self.nx]);
        self.x1des = np.zeros([self.nt,self.nx]);
        self.x2des = np.zeros([self.nt,self.nx]);
        self.udes = np.zeros([self.nt,self.nu]);
        self.xbase = np.zeros([self.nt,self.nx]);
        self.u1base = np.zeros([self.nt,self.nu1]);
        self.u2base = np.zeros([self.nt,self.nu2]);

        if 'x0' in ins: self.x0 = ins['x0'];
        if 'xdes' in ins: self.xdes = ins['xdes'];
        # if 'x1des' in ins: self.x1des = ins['x1des'];
        # if 'x2des' in ins: self.x2des = ins['x2des'];
        if 'udes' in ins: self.udes = ins['udes'];

        # self.x1des = self.xdes[:,:self.nx1];
        # self.x2des = self.xdes[:,self.nx1:];
        self.u1des = self.udes[:,:self.nu1];
        self.u2des = self.udes[:,self.nu1:];

        if 'x1des' in ins: self.x1des = ins['x1des']
        if 'x2des' in ins: self.x2des = ins['x2des']
        if 'u1des' in ins: self.u1des = ins['u1des']
        if 'u2des' in ins: self.u2des = ins['u2des']

        if 'xbase' in ins: self.xbase = ins['xbase'];
        if 'u1base' in ins: self.u1base = ins['u1base'];
        if 'u2base' in ins: self.u2base = ins['u2base'];

        self.xdes_flat = self.xdes.reshape(self.nt*self.nx);
        self.x1des_flat = self.x1des.reshape(self.nt*self.nx);
        self.x2des_flat = self.x2des.reshape(self.nt*self.nx);

        self.xbase_flat = self.xbase.reshape(self.nt*self.nx);
        self.u1des_flat = self.u1des.reshape(self.nt*self.nu1);
        self.u2des_flat = self.u2des.reshape(self.nt*self.nu2);        
        self.u1base_flat = self.u1base.reshape(self.nt*self.nu1);
        self.u2base_flat = self.u2base.reshape(self.nt*self.nu2);

        self.obs = {};
        if 'obs' in ins: self.obs = ins['obs'];

        #### EXTENDED NONLINEAR STUFF
        self.nlcosts = {}; self.nlinits = {};
        if 'nlcosts' in ins: self.nlcosts = ins['nlcosts'];
        if 'nlinits' in ins: self.nlinits = ins['nlinits'];
        self.utrajs = {};
        self.gtrajs = {};
        # self.grad = [];
        # self.hess = np.zeros([self.nt,self.nt]);
        self.d1 = self.nt*self.nu1;
        self.d2 = self.nt*self.nu2;
        self.dd = self.d1 + self.d2;

    def gradNL(self,x,u1,u2,params={},shape='flat'):
        if shape == 'flat':
            xx = x.reshape([self.nt,self.nx])
            uu1 = u1.reshape([self.nt,self.nu1])
            uu2 = u2.reshape([self.nt,self.nu2])
            grad1x = [];
            grad2x = []; 
            grad1u1 = [];
            grad2u2 = [];
            for t in range(self.nt):
                temp = {}
                temp['state'] = xx[t];
                temp['control'] = np.hstack([uu1[t],uu2[t]]);
                temp['time'] = t;
                temp['tstep'] = t;
                temp['params'] = self.nlcosts['params'];
                grad1x.append(self.nlcosts['dl1dx'](temp));
                grad2x.append(self.nlcosts['dl2dx'](temp));
                grad1u1.append(self.nlcosts['dl1du1'](temp));
                grad2u2.append(self.nlcosts['dl2du2'](temp));
            grad1x = np.hstack(grad1x);
            grad2x = np.hstack(grad2x);
            grad1u1 = np.hstack(grad1u1);
            grad2u2 = np.hstack(grad2u2);
            gradu1 = grad1u1 + grad1x@self.AA@self.BB1; # + 
            gradu2 = grad2u2 + grad2x@self.AA@self.BB2; # 
            out = np.hstack([gradu1,gradu2]); #{'u1':gradu1,'u2':uu2}
        return out; #

    # def hessNL(self,x,u,shape='flat'):
    #     if shape == 'flat':
    def basicEquilibriaNL(self,typs=['nash']): #,theta1=None,theta2=None):
        ### requires a
        if len(self.nlcosts)>0:
            if 'nash' in typs:
                maxiters = self.nlcosts['maxiters'];
                stepsize = self.nlcosts['stepsize'];
                if len(self.nlinits): u0s = self.nlinits['nash'];
                else: u0s = np.zeros([1,self.dd]);
                self.utrajs['nash'] = np.zeros([len(u0s),maxiters,self.dd]);
                self.gtrajs['nash'] = np.zeros([len(u0s),maxiters,self.dd]);
                for j,u0 in enumerate(u0s):
                    self.utrajs['nash'][j][0] = u0;
                    for k in range(maxiters-1):
                        u1 = self.utrajs['nash'][j][k][:self.d1]
                        u2 = self.utrajs['nash'][j][k][self.d1:]
                        xx = self.AA@(self.BB1@u1 + self.BB2@u2) + self.HH@self.x0
                        grad = self.gradNL(xx,u1,u2)
                        self.gtrajs['nash'][j][k] = grad;
                        self.utrajs['nash'][j][k+1] = self.utrajs['nash'][j][k] - stepsize*grad;

    def buildRollMats(self):
        self.AAi = smat.block_diag(*[np.eye(self.nx) for _ in range(self.nt)])
        for t in range(self.nt)[1:]:
            self.AAi[self.nx*t:self.nx*(t+1),self.nx*(t-1):self.nx*t] = -self.As[t-1];
        self.AA = mat.inv(self.AAi);

        HH = [np.eye(self.nx)];
        for t in range(self.nt):
            temp = self.As[t]@HH[t];
            HH[t] = temp.copy()
            if t < self.nt-1: HH.append(temp.copy())
        self.HH = np.vstack(HH);

        if len(self.Bs)>0: self.BB = smat.block_diag(*list(self.Bs));
        if len(self.B1s)>0: self.BB1 = smat.block_diag(*list(self.B1s));
        if len(self.B2s)>0: self.BB2 = smat.block_diag(*list(self.B2s));

        # if len(self.Qs)>0: self.QQ = smat.block_diag(*self.Qs)
        if len(self.Q1s)>0: self.QQ1 = smat.block_diag(*self.Q1s)
        if len(self.Q2s)>0: self.QQ2 = smat.block_diag(*self.Q2s)
        if len(self.P1s)>0: self.PP1 = smat.block_diag(*self.P1s)
        if len(self.P2s)>0: self.PP2 = smat.block_diag(*self.P2s)
        # if len(self.Rs)>0: self.RR = smat.block_diag(*self.Rs)
        if len(self.R1s)>0: self.RR1 = smat.block_diag(*self.R1s)
        if len(self.R2s)>0: self.RR2 = smat.block_diag(*self.R2s)
        if len(self.R12s)>0: self.RR12 = smat.block_diag(*self.R12s)
        if len(self.R21s)>0: self.RR21 = smat.block_diag(*self.R21s)


        self.xdiff = self.x0@self.HH.T-self.xdes_flat; #.reshape(self.nt*self.nx)
        self.xdrift = self.x0@self.HH.T

        out = {}
        out['FF'] = self.AA; out['FFi'] = self.AAi;
        out['HH'] = self.HH;
        out['GG1'] = self.BB1; out['GG2'] = self.BB2;
        out['QQ1'] = self.QQ1; out['QQ2'] = self.QQ2;
        out['PP1'] = self.PP1; out['PP2'] = self.PP2;
        out['RR1'] = self.RR1; out['RR2'] = self.RR2;
        out['RR12'] = self.RR12; out['RR21'] = self.RR21;
        out['x0'] = self.x0;
        out['xdiff'] = self.xdiff;
        out['xdes'] = self.xdes;
        out['udes'] = self.udes;
        out['u1des'] = self.u1des;
        out['u2des'] = self.u2des;
        out['nT'] = self.nt;
        out['nx'] = self.nx;
        out['nx1'] = self.nx1;
        out['nx2'] = self.nx2;
        return out

class QUADGAME2:
    def __init__(self,ins={}):

        self.build_version = 'from_data'
        if 'version' in ins: self.build_version = ins['version'];

        self.LINQUAD = None
        if 'LQ' in ins: self.LINQUAD = ins['LQ'];
        if self.build_version == 'from_lq':
            self.fromLinQuad();
        if self.build_version == 'from_data':
            self.fromData(ins);
        self.exptags = ['exp1','exp2','exp3','exp4']; #,'exp3','exp4'];

        #### EXTENDED NONLINEAR STUFF   
        self.GRAD = {}
        if 'GRAD' in ins: self.GRAD = ins['GRAD'];
        self.utrajs = {};


    def fromData(self,ins={}):
                    # initialize dimensions
        selfx = 1; self.d2 = 1;
        if 'd1' in ins: self.d1 = ins['d1'];
        if 'd2' in ins: self.d2 = ins['d2'];


        self.A1 = np.eye(self.d1);
        self.A2 = np.eye(self.d2);
        self.B1 = np.zeros([self.d2,self.d1]);
        self.B2 = np.zeros([self.d1,self.d2]);
        self.D1 = np.eye(self.d2);
        self.D2 = np.eye(self.d1);
        self.a1 = np.zeros(self.d1);
        self.b1 = np.zeros(self.d2);
        self.a2 = np.zeros(self.d1);
        self.b2 = np.zeros(self.d2);

        if 'A1' in ins: self.A1 = ins['A1']
        if 'A2' in ins: self.A2 = ins['A2']
        if 'B1' in ins: self.B1 = ins['B1']
        if 'B2' in ins: self.B2 = ins['B2']
        if 'D1' in ins: self.D1 = ins['D1']
        if 'D2' in ins: self.D2 = ins['D2']
        if 'a1' in ins: self.a1 = ins['a1']
        if 'b1' in ins: self.b1 = ins['b1']
        if 'a2' in ins: self.a2 = ins['a2']
        if 'b2' in ins: self.b2 = ins['b2']



        self.d1 = len(self.A1);
        self.d2 = len(self.A2);
        self.M1 = np.block([[self.A1,self.B1.T],[self.B1,self.D1]])
        self.M2 = np.block([[self.D2,self.B2],[self.B2.T,self.A2]])
        self.c1 = np.hstack([self.a1,self.b1]);
        self.c2 = np.hstack([self.b2,self.a2]);


        if 'M1' in ins: self.M1 = ins['M1'];
        if 'M2' in ins: self.M2 = ins['M2'];        
        if 'c1' in ins: self.c1 = ins['c1'];
        if 'c2' in ins: self.c2 = ins['c2'];


        self.center1 = -mat.inv(self.M1)@self.c1
        self.center2 = -mat.inv(self.M2)@self.c2
        if 'center1' in ins: self.center1 = ins['center1'];
        if 'center2' in ins: self.center2 = ins['center2'];
        self.c1 = -self.M1@self.center1
        self.c2 = -self.M2@self.center2


        self.A1 = self.M1[:self.d1,:self.d1];
        self.A2 = self.M2[self.d1:,self.d1:];
        self.D1 = self.M1[self.d1:,self.d1:];
        self.D2 = self.M2[:self.d1,:self.d1];
        self.B1 = self.M1[self.d1:,:self.d1];
        self.B2 = self.M2[:self.d1,self.d1:];
        self.a1 = self.c1[:self.d1]
        self.b1 = self.c1[self.d1:]
        self.a2 = self.c2[self.d1:]
        self.b2 = self.c2[:self.d1]

        ### BUILD REMAINDER
        self.M1inv = mat.inv(self.M1)
        self.M2inv = mat.inv(self.M2)
        self.dd = self.d1 + self.d2;
        self.aa = np.hstack([self.a1,self.a2]);
        self.bb = np.hstack([self.b2,self.b1]);
        self.M = np.block([[self.A1,self.B2],[self.B1,self.A2]]);
        self.N = np.block([[self.D2,self.B1.T],[self.B2.T,self.D1]]);
        self.Minv = mat.inv(self.M);
        
        try: self.Ninv = mat.inv(self.N);
        except: self.Ninv = smat.pinv(self.N); print('WARNING: N not invertible')


    def fromLinQuad(self,ins={}):

        if not(self.LINQUAD==None):
        #     print('blarg')
            AA  = self.LINQUAD.AA;
            AAi = self.LINQUAD.AAi;
            HH  = self.LINQUAD.HH;
            BB1 = self.LINQUAD.BB1;
            BB2 = self.LINQUAD.BB2;
            QQ1 = self.LINQUAD.QQ1;
            QQ2 = self.LINQUAD.QQ2;            
            #### NEW ######
            PP1 = self.LINQUAD.PP1;
            PP2 = self.LINQUAD.PP2;
            ###############
            RR1 = self.LINQUAD.RR1;
            RR2 = self.LINQUAD.RR2;
            RR12 = self.LINQUAD.RR12;
            RR21 = self.LINQUAD.RR21;

            x0 = self.LINQUAD.x0;
            xdiff = self.LINQUAD.xdiff;
            xdes = self.LINQUAD.xdes_flat
            x1des = self.LINQUAD.x1des_flat;
            x2des = self.LINQUAD.x2des_flat;
            xbase = self.LINQUAD.xbase_flat;
            xdrift = self.LINQUAD.xdrift;

            # xdes1 = self.LINQUAD.xdes1_flat
            # xdes2 = self.LINQUAD.xdes2_flat

            ### From other object 
            # self.xdiff = self.x0@self.HH.T-self.xdes_flat; #.reshape(self.nt*self.nx)
            # self.xdrift = self.x0@self.HH.T
            u1base = self.LINQUAD.u1base_flat;
            u2base = self.LINQUAD.u2base_flat;
            u1des = self.LINQUAD.u1des_flat;
            u2des = self.LINQUAD.u2des_flat;

            # print(x1des)
            # print('')
        
            # asdfasdfa
            ### NEW>>> 
            # xconst = AA@BB1@u1base + AA@BB2@u2base + xdrift + xbase - xdes;



            xconst = xbase + xdrift - xdes;
            x1const = xbase + xdrift - x1des;
            x2const = xbase + xdrift - x2des;

            xauto = xbase + xdrift; 

            self.A1 = RR1 + BB1.T@AA.T@(QQ1+PP1)@AA@BB1;
            self.A2 = RR2 + BB2.T@AA.T@(QQ2+PP2)@AA@BB2;
            self.B1 = BB2.T@AA.T@(QQ1+PP1)@AA@BB1;
            self.B2 = BB1.T@AA.T@(QQ2+PP2)@AA@BB2;
            self.D1 = RR12 + BB2.T@AA.T@(QQ1+PP1)@AA@BB2;
            self.D2 = RR21 + BB1.T@AA.T@(QQ2+PP2)@AA@BB1;
            # self.D1 = BB2.T@AA.T@QQ1@AA@BB2;
            # self.D2 = BB1.T@AA.T@QQ2@AA@BB1;

            ### constants - only without obstacles...
            # self.offset1 = (u1base-u1des)@RR1@(u1base-u1des) + xconst@QQ1@xconst
            # self.offset2 = (u2base-u2des)@RR2@(u2base-u2des) + xconst@QQ2@xconst
            self.a1 = (xauto-x1des)@QQ1@AA@BB1 + xauto@PP1@AA@BB1 #+ (u1base-u1des)@RR1;\
            self.b1 = (xauto-x1des)@QQ1@AA@BB2 + xauto@PP1@AA@BB2 #+ (u2base-u2des)@RR12;
            self.a2 = (xauto-x2des)@QQ2@AA@BB2 + xauto@PP2@AA@BB2 #+ (u2base-u2des)@RR2;
            self.b2 = (xauto-x2des)@QQ2@AA@BB1 + xauto@PP2@AA@BB1 #+ (u1base-u1des)@RR21;


            #### ADDING OBSTACLES #### ADDING OBSTACLES #### ADDING OBSTACLES #### ADDING OBSTACLES 
            #### ADDING OBSTACLES #### ADDING OBSTACLES #### ADDING OBSTACLES #### ADDING OBSTACLES 
            #### ADDING OBSTACLES #### ADDING OBSTACLES #### ADDING OBSTACLES #### ADDING OBSTACLES 

            AO1 = np.zeros(self.A1.shape); BO1 = np.zeros(self.B1.shape); DO1 = np.zeros(self.D1.shape)
            aO1 = np.zeros(self.a1.shape); bO1 = np.zeros(self.b1.shape)
            AO2 = np.zeros(self.A2.shape); BO2 = np.zeros(self.B2.shape); DO2 = np.zeros(self.D2.shape);
            aO2 = np.zeros(self.a2.shape); bO2 = np.zeros(self.b2.shape);
            # if len(self.LINQUAD.obs)>0:

            #     nt = self.LINQUAD.nt;
            #     obs1 = self.LINQUAD.obs['p1'];
            #     obs2 = self.LINQUAD.obs['p2'];
            #     vshape1 = np.outer(np.ones(nt),np.array([1.,0.])).reshape(nt*2);
            #     vshape2 = np.outer(np.ones(nt),np.array([0.,1.])).reshape(nt*2);

            #     obs = obs1;
            #     for _,tag in enumerate(obs):
            #         pt = obs[tag]['pt']
            #         SHP = obs[tag]['shape'];
            #         PP = np.kron(np.diag(vshape1),SHP)
            #         xconst = xbase + xdrift - np.kron(vshape1,pt)
            #         AO1 = AO1 - BB1.T@AA.T@PP@AA@BB1
            #         BO1 = BO1 + BB2.T@AA.T@PP@AA@BB1
            #         DO1 = DO1 - BB2.T@AA.T@PP@AA@BB2
            #         aO1 = aO1 - xconst@PP@AA@BB1
            #         bO1 = bO1 - xconst@PP@AA@BB2

            #     obs = obs2;
            #     for _,tag in enumerate(obs):
            #         pt = obs[tag]['pt']
            #         SHP = obs[tag]['shape'];
            #         PP = np.kron(np.diag(vshape2),SHP)
            #         xconst = xbase + xdrift - np.kron(vshape2,pt)
            #         AO2 = AO2 - BB2.T@AA.T@PP@AA@BB2
            #         BO2 = BO2 + BB1.T@AA.T@PP@AA@BB2
            #         DO2 = DO2 - BB1.T@AA.T@PP@AA@BB1
            #         aO2 = aO2 - xconst@PP@AA@BB2
            #         bO2 = bO2 - xconst@PP@AA@BB1


            #     self.A1 = self.A1 + AO1
            #     self.A2 = self.A2 + AO2
            #     self.B1 = self.B1 + BO1
            #     self.B2 = self.B2 + BO2
            #     self.D1 = self.D1 + DO1
            #     self.D2 = self.D2 + DO2

            #     self.a1 = self.a1 + aO1
            #     self.a2 = self.a2 + aO2
            #     self.b1 = self.b1 + bO1
            #     self.b2 = self.b2 + bO2

            ##### FINISHED OBSTACLES...##### FINISHED OBSTACLES...##### FINISHED OBSTACLES...##### FINISHED OBSTACLE
            ##### FINISHED OBSTACLES...##### FINISHED OBSTACLES...##### FINISHED OBSTACLES...##### FINISHED OBSTACLE
            ##### FINISHED OBSTACLES...##### FINISHED OBSTACLES...##### FINISHED OBSTACLES...##### FINISHED OBSTACLE


            self.M1 = np.block([[self.A1,self.B1.T],[self.B1,self.D1]])
            self.M2 = np.block([[self.D2,self.B2],[self.B2.T,self.A2]])

            self.M1inv = mat.inv(self.M1)
            self.M2inv = mat.inv(self.M2)

            self.c1 = np.hstack([self.a1,self.b1]);
            self.c2 = np.hstack([self.b2,self.a2]);

            self.d1 = self.A1.shape[0];
            self.d2 = self.A2.shape[0]
            self.dd = self.d1 + self.d2;
            self.aa = np.hstack([self.a1,self.a2]);
            self.bb = np.hstack([self.b2,self.b1]);
            self.M = np.block([[self.A1,self.B2],[self.B1,self.A2]]);
            self.N = np.block([[self.D2,self.B1.T],[self.B2.T,self.D1]]);
            self.Minv = mat.inv(self.M);

            try: self.Ninv = mat.inv(self.N);
            except: self.Ninv = smat.pinv(self.N); print('WARNING: N not invertible')

    def svoEquilibria(self,theta1=0.,theta2=0.):


        cth1 = np.cos(theta1); sth1 = np.sin(theta1);
        cth2 = np.cos(theta2); sth2 = np.sin(theta2);
        a1 = self.a1; b1 = self.b1;
        a2 = self.a2; b2 = self.b2;
        A1 = self.A1; B1 = self.B1; D1 = self.D1;
        A2 = self.A2; B2 = self.B2; D2 = self.D2;

        cc = np.hstack([cth1*a1+sth1*b2, cth2*a2+sth2*b1]);

        MM = np.block([[cth1*A1+sth1*D2  ,cth2*B2 + sth2*B1.T],
                       [cth1*B1+sth1*B2.T,cth2*A2 + sth2*D1  ]]);

        uth = -cc@mat.inv(MM)
        return uth

    def basicEquilibria(self):        

        self.u1 = -self.c1@self.M1inv;
        self.u2 = -self.c2@self.M2inv;
        self.une  = -self.aa@self.Minv;
        self.una  = -self.bb@self.Ninv;

        self.uso = -(self.c1+self.c2)@mat.inv(self.M1+self.M2)

        # PHIne = smat.block_diag(phi2*np.eye(d1),psi1*np.eye(d2))
        # PHIso = smat.block_diag((phi2-1.)*np.eye(d1),(psi1-1.)*np.eye(d2))


        # self.uneb = -(self.bb+self.une@self.N)@PHIne@mat.inv(self.M+self.N@PHIne);
        # self.cso = self.c1 + self.c2;
        # self.Mso = self.M1 + self.M2;
        # self.Msoinv = mat.inv(self.Mso);

        ##############################
        # self.uso  = -self.cso@self.Msoinv;
        # self.usob = -(self.bb+self.uso@self.N)@PHIso@mat.inv(self.Mso+self.N@PHIso)

        self.xxi_exp1 = self.una - self.une;
        self.xxi_exp2 = self.u2 - self.u1;

        self.une1 = self.une[:self.d1];
        self.une2 = self.une[self.d1:];
        # self.uneb1 = self.uneb[:nu1*nt];
        # self.uneb2 = self.uneb[nu1*nt:]
        # self.uso1 = self.uso[:nu1*nt];
        # self.uso2 = self.uso[nu1*nt:]
        # self.usob1 = self.usob[:nu1*nt];
        # self.usob2 = self.usob[nu1*nt:]


        # eigs,_ = mat.eig(A2+10*D1);
        # eigs,_ = mat.eig(A1+D2);
        self.eigsM1,_ = mat.eig(self.M1);
        self.eigsM2,_ = mat.eig(self.M2);
        self.eigsA1,_ = mat.eig(self.A1);
        self.eigsA2,_ = mat.eig(self.A2);
        self.eigsD1,_ = mat.eig(self.D1);
        self.eigsD2,_ = mat.eig(self.D2);
        self.eigsM1M2,_ = mat.eig(self.M1@self.M2inv);
        self.eigsM2M1,_ = mat.eig(self.M2@self.M1inv);

        eps = 0.0;
        self.eigsM1M2p = self.eigsM1M2[self.eigsM1M2 > eps] 
        self.eigsM1M2m = self.eigsM1M2[self.eigsM1M2 <= eps];
        self.eigsM2M1p = self.eigsM1M2[self.eigsM2M1 > eps] 
        self.eigsM2M1m = self.eigsM1M2[self.eigsM2M1 <= eps];
        self.eigsM1p = self.eigsM1[self.eigsM1 > eps] 
        self.eigsM1m = self.eigsM1[self.eigsM1 <= eps];
        self.eigsM2p = self.eigsM2[self.eigsM2 > eps] 
        self.eigsM2m = self.eigsM2[self.eigsM2 <= eps];


        if not(self.LINQUAD==None):

            AA = self.LINQUAD.AA;
            BB1 = self.LINQUAD.BB1;
            BB2 = self.LINQUAD.BB2;
            HH = self.LINQUAD.HH;
            x0 = self.LINQUAD.x0
            self.xopt1 = AA@np.hstack([BB1,BB2])@self.u1 + HH@x0;
            self.xopt2 = AA@np.hstack([BB1,BB2])@self.u2 + HH@x0;
            self.xne = AA@np.hstack([BB1,BB2])@self.une + HH@x0;
            self.xna = AA@np.hstack([BB1,BB2])@self.una + HH@x0;
            self.xso = AA@np.hstack([BB1,BB2])@self.uso + HH@x0;

        # eigsA1 = svoGame['eigsA1']; eigsA2 = svoGame['eigsA2'];
        # eigsD1 = svoGame['eigsD1']; eigsD2 = svoGame['eigsD2']


    #### EXTENDED NONLINEAR STUFF 
    # def basicEquilibriaNL(self,theta1=None,theta2=None):
    #     ### requires a 
    #     if len(GRAD)>0:
    #         maxiters = self.GRAD['maxiters'];
    #         stepsize = self.GRAD['stepsize'];
    #         df1dx1 = self.GRAD['df1dx1'];
    #         df2dx2 = self.GRAD['df2dx2'];
    #         u0s = self.GRAD['u0s'];
    #         self.utrajs['nash'] = np.zeros([len(u0s),maxiters,self.dd])
    #         for j,u0 in enumerate(u0s):
    #             self.utrajs['nash'][j][0] = u0;
    #             for k in range(maxiters-1):
    #                 unk = self.utrajs['nash'][k]
    #                 grad = np.hstack([df1dx1(unk),df2dx2(unk)])

    def computeFreqStruct(self,version = 'simple1',ins={}):
        pass
    #                 self.utrajs['nash'][k+1] = self.utrajs['nash'][k] - stepsize*grad;

    def computeFitting(self,version='simple1',ins={}):

        if version == 'simple1':
            
            th1 = np.pi/4; th2 = np.pi/4;

            # th1 = (np.pi/2)*np.random.rand(1)[0];
            # th2 = (np.pi/2)*np.random.rand(1)[0];
            if 'th1' in ins: th1 = ins['th1'];
            if 'th2' in ins: th2 = ins['th2'];
            AA  = self.LINQUAD.AA;
            BB1 = self.LINQUAD.BB1;
            BB2 = self.LINQUAD.BB2;
            params = self.LINQUAD.paramsIndAgents
            nt = self.LINQUAD.nt

            tanth1 = np.tan(th1);
            tanth2 = np.tan(th2);
            cotth1 = (1./tanth1)
            cotth2 = (1./tanth2)

            # params = {'r1':1.,'q1':1.,'p1':1.,'r1b':1.,'q1b':1.,'r2':1.,'q2':1.,'p2':1.,'r2b':1.,'q2b':1.}
            r1 = params['r1']; 
            q1 = params['q1'];
            p1 = params['p1'];
            r2 = params['r2'];
            q2 = params['q2'];
            p2 = params['p2'];
            r1b = params['r1b'];
            q1b = params['q1b'];
            r2b = params['r2b'];
            q2b = params['q2b'];

            # print('r1: ', r1)
            # print('r2: ', r2)
            # print('q1: ', q1)
            # print('q2: ', q2)
            # print('p1: ', p1)
            # print('p2: ', p2)
            # print('')
            # print('r1b: ', r1b)
            # print('r2b: ', r2b)
            # print('q1b: ', q1b)
            # print('q2b: ', q2b)


            Qform = np.eye(2);
            Rform = np.eye(2);

            QQ1a = np.kron(np.eye(nt),smat.block_diag(Qform,np.zeros((2,2))));
            QQ1b = np.kron(np.eye(nt),smat.block_diag(np.zeros((2,2)),Qform));
            QQ2a = np.kron(np.eye(nt),smat.block_diag(np.zeros((2,2)),Qform));
            QQ2b = np.kron(np.eye(nt),smat.block_diag(Qform,np.zeros((2,2))));

            MM = BB1.T@AA.T@QQ1a@AA@BB1;
            _,sigs2,VT = smat.svd(MM); V = VT.T
            #sigs2 = sigs*sigs;
            sigs4 = sigs2*sigs2;
            sigs = np.sqrt(sigs2);
            

            PERM1 = np.kron(np.eye(int(self.dd/2)),np.array([1,0]))
            PERM2 = np.kron(np.eye(int(self.dd/2)),np.array([0,1]))
            PERM = np.vstack([PERM1,PERM2]);
            VV = smat.block_diag(V,V)@PERM;
            uth = self.svoEquilibria(theta1=th1,theta2=th2);
            u2 = self.u2; u1 = self.u1;
            M1 = self.M1; M2 = self.M2;
            M1i = mat.inv(self.M1)

            # print('nt = ',nt)
            # print('half dd = ',int(self.dd/2))
            # print('d1 ',self.d1)
            # print('d2 ',self.d2)
            d1 = self.d1
            if False: ### initial test 
                HH1 = np.block([[cotth1*np.eye(d1),np.zeros([d1,d1])],[np.zeros([d1,d1]),1.0*np.eye(d1)]]);
                HH2 = np.block([[1.0*np.eye(d1),np.zeros([d1,d1])],[np.zeros([d1,d1]),cotth2*np.eye(d1)]]);                
                diff = HH1@M1@(uth-u1) - HH2@M2@(u2-uth)
                print(diff)
            elif True: # linear test 
                HH1 = np.block([[cotth1*np.eye(d1),np.zeros([d1,d1])],[np.zeros([d1,d1]),1.0*np.eye(d1)]]);
                uL = VV.T@HH1@M1@(uth-u1);
                uR = VV.T@(u2-uth)
                UR = uR.reshape([self.d1,2]);
            
            if True: 
                YY1 = VV.T@(uth-u1)
                YY2 = VV.T@(uth-u2)
                YY1 = YY1.reshape([self.d1,2]);
                YY2 = YY2.reshape([self.d1,2]);
            # elif False: 
            #     HH1 = np.kron(np.eye(self.d1),np.array([[tanth1,0.],[0.,1.]]))
            #     uR = HH1@VV.T@M1i@(u2-uth)
            #     uL = VV.T@(uth - u1);
            #     UR = uR.reshape([self.d1,2]);

            Ybar1 = np.zeros([0,6]);
            Ybar2 = np.zeros([0,6]);
            for k in range(self.d1):
                chunk1 = np.array([[YY1[k][0],0.],[0.,YY1[k][1]]]);
                chunk2 = sigs2[k]*np.kron(np.eye(2),YY1[k]); 
                Ybar1 = np.vstack([Ybar1,np.hstack([chunk1,chunk2])])
                chunk1 = np.array([[YY2[k][0],0.],[0.,YY2[k][1]]]);
                chunk2 = sigs2[k]*np.kron(np.eye(2),YY2[k]); 
                Ybar2 = np.vstack([Ybar2,np.hstack([chunk1,chunk2])])

            MAT = np.zeros([0,6])
            for k in range(self.d1):
                chunk1 = np.array([[UR[k][0],0.],[0.,UR[k][1]]]);
                chunk2 = sigs2[k]*np.kron(np.eye(2),UR[k]); #np.block([[UR[k],np.zeros(2)],[np.zeros(2),UR[k]]])
                MAT = np.vstack([MAT,np.hstack([chunk1,chunk2])])

            pars1 = np.array([r1*cotth1,r1b,(q1-p1)*cotth1,cotth1*p1,p1,q1b-p1])
            pars2 = np.array([r2b,r2*cotth2,q2b-p2,p2,p2*cotth2,(q2-p2)*cotth2])
            pars  = np.hstack([pars1,pars2]);


            _,temp,_ = mat.svd(MAT)
            thresh = 0.00001

            Ybar = np.hstack([Ybar1,Ybar2]);
            deps1 = [0,2,3,6,8,9];
            deps2 = [1,4,5,7,10,11];


            blah = [];
            for _ in range(20):
                ind1 = random.sample(deps1,1)[0]    
                ind2 = random.sample(deps2,1)[0]
                delinds = [ind1,ind2]; takeinds = [];
                for k in range(12):
                    if not(k in delinds): takeinds.append(k)
                nYbar = np.delete(Ybar,delinds,1)
                parshat = mat.pinv(nYbar)@(-Ybar[:,delinds]@pars[delinds])
                diff = parshat - pars[takeinds];
                blah.append(mat.norm(diff))
            #print(blah)

            # print(mat.matrix_rank(Ybar[:,[2,3,8,9]]))
            # print(mat.matrix_rank(Ybar[:,[4,5,10,11]]))
            # print('')
            # print(mat.matrix_rank(Ybar[:,[0,2,3,8,9]]))
            # print(mat.matrix_rank(Ybar[:,[4,5,7,10,11]]))
            # print('')
            # print(mat.matrix_rank(Ybar[:,[2,3,6,8,9]]))
            # print(mat.matrix_rank(Ybar[:,[1,4,5,10,11]]))
            print('')
            print(mat.matrix_rank(Ybar[:,[0,2,3,6,8,9]]))
            print(mat.matrix_rank(Ybar[:,[1,4,5,7,10,11]]))
            # asdfasdf

            nYbar = [];
            nYbar.append(Ybar[:,0])
            nYbar.append(Ybar[:,2])
            nYbar.append(Ybar[:,3])
            nYbar.append(Ybar[:,4]-Ybar[:,5]);
            nYbar.append(Ybar[:,7])
            nYbar.append(Ybar[:,8]+Ybar[:,9])
            nYbar.append(Ybar[:,10])
            nYbar.append(Ybar[:,11])
            nYbar = np.vstack(nYbar).T;
            # for k in range(nYbar.shape[1]):
            #     temp = np.delete(nYbar,k,1)
            #     print(mat.matrix_rank(temp))

            npars1 = np.array([r1*cotth1,(q1-p1)*cotth1,cotth1*p1,p1])
            npars2 = np.array([r2*cotth2,p2,p2*cotth2,(q2-p2)*cotth2])
            npars = np.hstack([npars1,npars2]);

            print('reduced shape: ',nYbar.shape)
            print('reduced rank: ',mat.matrix_rank(nYbar))


            BLOB1 = Ybar[:,deps1]
            BLOB2 = Ybar[:,deps2];
            # _,_,temp = mat.svd(BLOB1); print(np.round(temp[-1],4));
            # _,_,temp = mat.svd(BLOB2); print(np.round(temp[-1],4));

            YBR1 = [];
            YBR1.append(Ybar1[:,0])
            YBR1.append(Ybar1[:,2])
            YBR1.append(Ybar1[:,3])
            YBR1.append(Ybar1[:,4] - Ybar1[:,5])
            YBR1 = np.vstack(YBR1).T

            YBR2 = [];
            YBR2.append(Ybar2[:,1])
            YBR2.append(-Ybar2[:,2]+Ybar2[:,3])
            YBR2.append(Ybar2[:,4])
            YBR2.append(Ybar2[:,5])
            YBR2 = np.vstack(YBR2).T
            


            if False: 
                print('')
                print('rank test:')
                print('rank1: ',mat.matrix_rank(BLOB1))
                print('rank2: ',mat.matrix_rank(BLOB2))

                scl1 = 1.
                scl2 = 1.

                print('test1: ', np.sum(np.abs(scl1*BLOB1@pars[deps1])))
                print('test2: ', np.sum(np.abs(scl2*BLOB2@pars[deps2])))




            # _,temp,_ = mat.svd(Ybar1);
            # print(np.sum(temp>thresh))
            # _,temp,_ = mat.svd(Ybar2);
            # print(np.sum(temp>thresh))
            
            blah = []
            val = 6;
            numfound = 0.; 
            if False: 
                for k in range(100):
                    Ybar = np.hstack([Ybar1,Ybar2]);
                    takeinds = random.sample(list(range(12)),val)
                    nYbar = Ybar[:,takeinds]
                    _,temp,_ = mat.svd(nYbar);
                    rnk = np.sum(temp>thresh)
                    if rnk < val: print(np.sort(takeinds)); numfound = numfound + 1;
                    blah.append(np.sum(temp > thresh))
                    #### reverse test
                    delinds = random.sample(list(range(12)),2)
                    nYbar = np.delete(Ybar,delinds,1)
                    # delinds = [int(12*np.random.rand(1)[0]) for _ in range(3)];
                    # for k in delinds: Ybar = np.delete(Ybar,k,1)
                    _,temp,_ = mat.svd(nYbar);
                    #blah.append(np.sum(temp > thresh))            
                print('num found: ',numfound)
                
            ### PRINTING TEST 1
            if False: 
                parshat2 = mat.pinv(MAT)@uL
                diff = parshat2 - pars2
                print(diff)

            if False: 
                diff = Ybar1@pars1 + Ybar2@pars2
                print(diff) # should be 0

            if False: 
                print('testing for removing player 1 variables...')
                for k in range(6):
                    chunk = np.delete(Ybar1,k,1)
                    MAT = np.hstack([chunk,Ybar2]);
                    parshat = mat.pinv(MAT)@(-Ybar1[:,k]*pars1[k])
                    pars = np.hstack([pars1,pars2])
                    pars = np.delete(pars,k,0)
                    diff = pars - parshat
                    _,temp,_ = mat.svd(MAT);
                    print('rank: ',np.sum(temp > thresh))
                    # print(temp)
                    # print(mat.norm(diff));
                print('')
                print('testing for removing player 2 variables...')
                for k in range(6):
                    chunk = np.delete(Ybar2,k,1)
                    MAT = np.hstack([Ybar1,chunk]);
                    parshat = mat.pinv(MAT)@(-Ybar2[:,k]*pars2[k])
                    pars = np.hstack([pars1,pars2])
                    pars = np.delete(pars,k,0)
                    diff = pars - parshat
                    _,temp,_ = mat.svd(MAT);
                    print('rank: ',np.sum(temp > thresh))
                    # print(temp)                    
                    # print(mat.norm(diff));

            print('doing extra random shit...')

            h2 = (u2-u1)@VV; ht = (uth-u1)@VV;
            h2 = h2.reshape([self.d1,2]);
            ht = ht.reshape([self.d1,2]);
            H2 = []; Ht = [];
            for k in range(self.d1):
                H2.append(np.kron(np.eye(2),h2[k]));
                Ht.append(np.kron(np.eye(2),ht[k]));
            H2 = np.vstack(H2)
            Ht = np.vstack(Ht)
            Hbar = np.hstack([-H2,Ht])
            xi = np.array([1,0,0,1]);

            S4 = np.kron(np.diag(sigs4),np.eye(2));
            S2 = np.kron(np.diag(sigs2),np.eye(2));
            HH = np.hstack([S4@Hbar,S2@Hbar,Hbar]);            
            hh = np.vstack([S4@Ht@xi,-S2@Ht@xi]).T@np.array([q1*p1,r1*p1])

            Xi4_2 = np.array([[0.0,   -p1,  p1*np.tan(th1), -np.tan(th1)],
                              [-p1,    p1, -p1*np.tan(th1),          0.0],
                              [0.0, q1-p1,  p1*np.tan(th1),          0.0],
                              [ r1,-q1+p1, -p1*np.tan(th1), -np.tan(th1)]]);

            Xi2_2 = np.array([[  0.0, 0.0,   0.0, -np.tan(th1)],
                              [  0.0, 0.0,   -p1,          0.0],
                              [  0.0,  r1,   0.0,          0.0],
                              [   r1, -r1, q1-p1, -np.tan(th1)]]);

            Xi4_t = np.outer(xi,np.array([q1-p1,-q1,  0.0,-np.tan(th1)]));
            Xi2_t = np.outer(xi,np.array([   r1,-r1,q1-p1,-np.tan(th1)]));
            
            Xi4 = np.vstack([Xi4_2,Xi4_t]);
            Xi2 = np.vstack([Xi2_2,Xi2_t]);
            Xi0 = np.array([np.hstack([np.array([0.,0.,0.,r1]),xi*r1])]).T;

            phi4 = np.array([q2*(1./np.tan(th2)),p2*(1./np.tan(th2)),p2,q2*p2*(1./np.tan(th2))])
            phi2 = np.array([q2*(1./np.tan(th2)),p2*(1./np.tan(th2)),r2*(1./np.tan(th2)),r2*p2*(1./np.tan(th2))])
            phi0 = r2*(1./np.tan(th2));

            phitru = np.hstack([phi4,phi2,phi0])

            QQ = HH@smat.block_diag(Xi4,Xi2,Xi0);
            phihat = mat.pinv(QQ)@hh

            phi4hat = phihat[:4];
            phi2hat = phihat[4:8];
            phi0hat = phihat[8];

            print('difference...')

            print(mat.norm(QQ@phitru - hh))
            print('')
            print('shape of Xi4: ',Xi4.shape,' with rank ',mat.matrix_rank(Xi4))
            print('shape of Xi2: ',Xi2.shape,' with rank ',mat.matrix_rank(Xi2))
            print('shape of HH: ',HH.shape, ' with rank ',mat.matrix_rank(HH))
            print('shape of QQ: ',QQ.shape, ' with rank ',mat.matrix_rank(QQ));

            print('')
            print('COMPUTING PARAMETERS...')
            print('')
            print('actual parameters: ')
            print('phi4: ', np.round(np.real(phi4),3))
            print('phi2: ', np.round(np.real(phi2),3))
            print('phi0: ', np.round(np.real(phi0),3))

            print('')
            print('estimated parameters: ')
            print('phi4: ', np.round(np.real(phi4hat),3))
            print('phi2: ', np.round(np.real(phi2hat),3))
            print('phi0: ', np.round(np.real(phi0hat),3))


            
        # return out


    # NNinv = mat.inv(NN);
    # nt1 = data['nph']; nt2 = data['nt']; puff1 = 0.001; puff2 = 0.001;
    # nt1h = int(nt1/2); nt2h = int(nt2/2)
    # maxt2 = 200;
    # #t1s = np.linspace(-np.pi/2+puff1,np.pi/2-puff1,nt1)
    # theta0 = 0.+puff1;
    # theta1 = np.pi/4;
    # theta2 = np.pi/2-puff1;
    # # theta2 = np.pi/3;
    # # theta3 = np.pi/2-puff1;
    # midt2 = 1.;
    # # t1s = np.linspace(theta0,theta1,nt1)
    # #t1s = np.linspace(theta1,theta2,nt1)
    # t1s_a = np.linspace(theta0,theta1,int(nt1/2))
    # t1s_b = np.linspace(theta1,theta2,int(nt1/2))
    # # # t1s_c = np.linspace(theta2,theta3,int(nt1/3))
    # t1s = np.hstack([t1s_a,t1s_b])
    # t2s_a = np.linspace(midt2,maxt2,int(nt2/2))
    # t2s_b = (1./t2s_a)[::-1]
    # t2s = np.hstack([t2s_b,t2s_a])


    def blowupsSVO(self,ins={}):
        self.nph = 10; self.nt = 20;
        if 'nph' in ins: self.nph = ins['nph'];
        if 'nt' in ins: self.nt = ins['nt'];

        bunchparams =  {'nt':10,'bunch1':True,'bunch2':True,'maxpreth1':2.,'maxpreth2':2.}
        bunchparams0 = {'nt':10,'bunch1':True,'bunch2':True,'maxpreth1':2.,'maxpreth2':2.}
        if 'bunchparams' in ins: bunchparams = ins['bunchparams'];
        if 'bunchparams0' in ins: bunchparams0 = ins['bunchparams0'];

        scale_evecs = 1.; 
        if 'scale_evecs' in ins: scale_evecs = ins['scale_evecs'];


        phbnds = [0.,np.pi/2];
        tbnds = (0,100);
        # tbnds0 = [0,100];

        if 'tbnds' in ins:
            tbnds = ins['tbnds']
        self.phbnds = phbnds;
        self.tbnds = tbnds;



        center1 = np.array([0.2,1.0]); center2 = np.array([1.0,0.2])

        dph = 0.2; 
        twin = [0.2,0.3]; phwin = [np.pi/4-dph,np.pi/4+dph];

        safe_buff = 0.05;
        sing_buff = 0.01;
        show_safe = True;
        show_sing = True;
        intstoshow = [0,1,2,3,4,5,6]; 
        show_samplepts_th1th2 = True;

        phwin = phwin
        twin =  twin; #[0.1,0.3]; #[0.3,1.3];  found: [2.2,2.4]; [2.8,3.]; [3.7,4.1]
        sing_num = 10;
        show_u1u2 = False;
        show_nean = True;
        nph = 10; nt = 100;

        ins2 = {};
        ins2['safe_buff'] = safe_buff
        ins2['sing_buff'] = sing_buff
        ins2['show_safe'] = show_safe
        ins2['show_sing'] = show_sing
        ins2['intstoshow'] = intstoshow
        ins2['show_samplepts_th1th2'] = show_samplepts_th1th2

        ins2['phbnds'] = phbnds
        ins2['tbnds'] = tbnds
        ins2['phwin'] = phwin
        ins2['twin'] =  twin; #[0.1,0.3]; #[0.3,1.3];  found: [2.2,2.4]; [2.8,3.]; [3.7,4.1]
        ins2['sing_num'] = sing_num
        # ins['rollMats'] = rollMats;
        # ins['svoGame'] = svoGame;
        ins2['show_u1u2'] = show_u1u2;
        ins2['show_nean'] = show_nean
        ins2['nph'] = nph;
        ins2['nt'] = nt;
        ins2['lims'] = {
                # 0:[-1,1,-1,1],
                1:[-0.2,1.1*np.pi/2,-0.2,1.1*np.pi/2]}



        colormap = 'viridis'
        # viridis, plasma, inferno, magma, cividis
        # spring, summer, autumn, winter
        # cool, hot, wistia
        # binary
        # hsv, gist_rainbow, rainbow, jet
        # good ones... 'viridis_r','plasma_r','cool_r'

        # ins['x0'] = x02player
        # ins['name'] = figname;

        ins2['aspect_equal'] = True;
        ins2['transparent'] = True;
        ins2['cmap'] = 'viridis'
        ins2['cmap1'] = 'winter'; #'plasma'; #
        ins2['cmap2'] = 'cool_r'; #'magma'; #

        # def blowups_lqgame(nph,nt,data):
        self.svo_safe_buff = ins2['safe_buff'];
        self.svo_sing_buff = ins2['sing_buff'];
        self.svo_singnum = ins2['sing_num'];
        self.svo_show_safe = ins2['show_safe'];
        self.svo_show_sing = ins2['show_sing'];

        self.phbnds = ins2['phbnds'];
        # self.tsbnds = ins['tsbnds'];
        # # mint = 0.;
        # # maxt = 100;
        # midt = 1; 
        # ts1 = np.linspace(midt,maxt,int(nt/2));
        # ts = np.hstack([(1./ts1)[::-1],ts1]);

        self.tbuff = 0.01
        self.mint = ins2['tbnds'][0]
        self.maxt = ins2['tbnds'][1];

        # print(self.mint)
        # print(self.maxt)
        
        self.minph = ins2['phbnds'][0]
        self.maxph = ins2['phbnds'][1];
        self.phbuff = 0.001
        self.ts = atan_linspace(self.mint+self.tbuff,self.maxt-self.tbuff,self.nt); #,invert=True)
        self.phs = np.linspace(self.minph+self.phbuff,self.maxph-self.phbuff,self.nph)

        if 'ts' in ins: self.ts = ins['ts']; self.nt = len(self.ts);
        if 'phs' in ins: self.phs = ins['phs']; self.nph = len(self.phs);
        self.nt4 = 6;        
        if 'nt4' in ins: self.nt4 = ins['nt4'];

        self.probpts = []
        self.fullpts = [];

        self.mineigA1 = np.min(self.eigsA1);
        self.mineigA2 = np.min(self.eigsA2);
        self.mineigD1 = np.min(self.eigsD1);
        self.mineigD2 = np.min(self.eigsD2);

        self.A1sqrt = smat.sqrtm(self.A1)
        self.A2sqrt = smat.sqrtm(self.A2)
        self.A1sqrti = mat.inv(self.A1sqrt);
        self.A2sqrti = mat.inv(self.A2sqrt);
        temp1,_ = mat.eig(self.A1sqrti.T@self.D2@self.A1sqrti)
        temp2,_ = mat.eig(self.A2sqrti.T@self.D1@self.A2sqrti)
        self.lammin1 = np.min(np.real(temp1))
        self.lammin2 = np.min(np.real(temp2))

        cot1 = -self.lammin1;
        cot2 = -self.lammin2;
        self.maxtheta1 = np.arctan2(1.,cot1);
        self.maxtheta2 = np.arctan2(1.,cot2);

        # print('max theta 1: ',self.maxtheta1)
        # print('max theta 2: ',self.maxtheta2)

        # # mineigA1 + t*np.cos(ph)*mineigD2 = 0
        # # mineigA2 + t*np.sin(ph)*mineigD1 = 0;
        # t1cut = np.abs(mineigA1)/(np.cos(ph)*np.abs(mineigD2));
        # t2cut = np.abs(mineigA2)/(np.sin(ph)*np.abs(mineigD1));
        # t = np.min([t1cut,t2cut]);
        # theta1 = np.arctan2(t*np.cos(ph),1)
        # theta2 = np.arctan2(t*np.sin(ph),1)
        # pdcutoffs2.append([theta1,theta2]);

        # t1cut = -1./(lammin1*np.cos(ph))
        # t2cut = -1./(lammin2*np.sin(ph))
        # t = np.abs(np.min([t1cut,t2cut]));
        # theta1 = np.arctan2(t*np.cos(ph),1)
        # theta2 = np.arctan2(t*np.sin(ph),1)
        # pdcutoffs.append([theta1,theta2]);
        # tcutoff = t;


        # print(x0)
        # print(HH)

        out = {};
        out['exp1'] = {}; out['exp2'] = {}; out['exp3'] = {}; out['exp4'] = {}

        ########################################################################################        
        ########################################################################################

        self.probpts = [];
        self.pdcutoffs = [];
        self.pdcutoffs2 = [];

        self.trajs = {}
        self.trajs['exp1'] = []
        self.trajs['exp2'] = []
        self.EXP = {};


        for tag in self.exptags:
            self.EXP[tag] = {};
            EXP = self.EXP[tag];

            EXP['tbuff'] = self.tbuff
            EXP['nt'] = self.nt;
            EXP['mint'] = self.mint; #data['tbnds'][0]
            EXP['maxt'] = self.maxt; #data['tbnds'][1];
            EXP['ts'] = self.ts; #atan_linspace(mint+tbuff,maxt-tbuff,nt); #,invert=True)
            EXP['minph'] = self.minph; #data['phbnds'][0]
            EXP['maxph'] = self.maxph; #data['phbnds'][1];
            EXP['phbuff'] = self.phbuff
            EXP['phs'] = self.phs; #np.linspace(self.minph+self.phbuff,self.maxph-self.phbuff,self.nph)            
            EXP['nph'] = self.nph;

            EXP['trajs'] = {};
            EXP['by_phs'] = {};            
            EXP['tbnds'] = tbnds;

            EXP['Pphs'] = {}
            EXP['Vs'] = {};
            EXP['Vis'] = {};
            EXP['Hph'] = {};
            EXP['Hphinv'] = {};
            EXP['eigs'] = {};
            EXP['neigs'] = {};
            EXP['fullpts'] = {};
            EXP['probts'] = {}
            EXP['probpts'] = {};
            EXP['probVis'] = {};
            EXP['probxVis'] = {};

            EXP['pdcutoffs'] = {};
            EXP['pdcutoffs2'] = {};
            EXP['fullpts'] = np.zeros([EXP['nph'],EXP['nt'],2]);

            for i,ph in enumerate(EXP['phs']):
                for j,t in enumerate(EXP['ts']):
                    theta1 = np.arctan2(t*np.cos(ph),1)
                    theta2 = np.arctan2(t*np.sin(ph),1)
                    EXP['fullpts'][i,j] = np.array([theta1,theta2])


            EXP['Mths'] = {}
            EXP['Nths'] = {}
            EXP['aths'] = {}
            EXP['bths'] = {}
            EXP['Mthis'] = {};
            EXP['uths'] = {};
            EXP['xxiths'] = {};                    


            EXP['pre_sing_trajs'] = {};
            EXP['post_sing_trajs'] = {};

            EXP['trajs'] = {};
            EXP['coeffs'] = {};
            EXP['ctrls'] = {};
            EXP['computed_curves'] = {}

            EXP['trajs_full'] = {}
            EXP['coeffs_full'] = {}
            EXP['ctrls_full'] = {};
            EXP['computed_curves_full'] = {}            

            EXP['coeffs_color_data'] = {};


                # self.trajs[ph][interval]['ts','ctrls']
                # self.coeffs
                # self.ctrls

            # EXP['Vs'] = {ph:np.zeros([self.dd,self.dd],dtype='complex') for ph in EXP['phs']}
            # EXP['Vis'] = {ph:np.zeros([self.dd,self.dd],dtype='complex') for ph in EXP['phs']}
            # EXP['Hph'] = {ph:np.zeros([self.dd]) for ph in EXP['phs']}
            # EXP['Hphinv'] = {ph:np.zeros([self.dd]) for ph in EXP['phs']};
            # EXP['eigs'] = {ph:np.zeros([self.dd],dtype='complex') for ph in EXP['phs']};
            # EXP['neigs'] = {}; #[[] for _ in range(EXP['nph'])];
            # EXP['EIGS'] = np.zeros([EXP['nph'],self.dd],dtype='complex');

            # # trajs = {}; 
            # # pre_sing_trajs = {};
            # post_sing_trajs = {};
            tbuff = 0.1;

            EXP['computed_curves']['ucurves'] = {}
            EXP['computed_curves']['thcurves'] = {}
            EXP['computed_curves']['tcurves'] = {}

            EXP['computed_curves_full']['ucurves'] = {}
            EXP['computed_curves_full']['thcurves'] = {}
            EXP['computed_curves_full']['tcurves'] = {}

            for ph in EXP['phs']:

                EXP['trajs'][ph] = {};
                pdcutoffs = [];
                pdcutoffs2 = [];
                EXP['pre_sing_trajs'][ph] = {};
                EXP['post_sing_trajs'][ph] = {};

                t1cut = -1./(self.lammin1*np.cos(ph))
                t2cut = -1./(self.lammin2*np.sin(ph))
                tcutoff = np.min([t1cut,t2cut]);
                if tcutoff < 0: tcutoff = 1000000.;
                # tcutoff = np.abs();
                theta1 = np.arctan2(tcutoff*np.cos(ph),1)
                theta2 = np.arctan2(tcutoff*np.sin(ph),1)
                pdcutoffs.append([theta1,theta2]);

                EXP['pdcutoffs'][ph] = np.array([theta1,theta2])

                # mineigA1 + t*np.cos(ph)*mineigD2 = 0
                # mineigA2 + t*np.sin(ph)*mineigD1 = 0;
                # t1cut = np.abs(mineigA1)/(np.cos(ph)*np.abs(mineigD2));
                # t2cut = np.abs(mineigA2)/(np.sin(ph)*np.abs(mineigD1));
                # t = np.min([t1cut,t2cut]);
                # theta1 = np.arctan2(t*np.cos(ph),1)
                # theta2 = np.arctan2(t*np.sin(ph),1)
                # pdcutoffs2.append([theta1,theta2]);


                EXP['Hph'][ph] = np.hstack([np.cos(ph)*np.ones(self.d1),np.sin(ph)*np.ones(self.d2)]);
                EXP['Hphinv'][ph] = np.hstack([(1./np.cos(ph))*np.ones(self.d1),(1./np.sin(ph))*np.ones(self.d2)]);
                Hph = np.diag(EXP['Hph'][ph]);
                Hphinv = np.diag(EXP['Hphinv'][ph]);
                # Ha = np.sqrt(2)*np.diag([np.cos(t1),np.sin(t1)]);
                # Hainv = mat.inv(Ha)                    
                # Ha = np.sqrt(2)*np.diag([np.cos(t1),np.sin(t1)]);
                # Hainv = mat.inv(Ha)


                Mth = None; Nth = None;
                ath = None; bth = None;
                Mthi = None; uth = None; xxith = None;

                if tag == 'exp1':
                    eigs,VV = mat.eig(self.M@Hphinv@self.Ninv)
                    VVi = mat.inv(VV)
                    xxiVV = self.xxi_exp1@VV;
                    ubase = self.une;

                elif tag == 'exp2':
                    eigs,VV = mat.eig(self.M1@Hphinv@self.M2inv)
                    VVi = mat.inv(VV)
                    xxiVV = self.xxi_exp2@VV;
                    ubase = self.u1

                elif tag == 'exp3':
                    # Ha = np.diag([np.cos(t1),np.sin(t1)]);    
                    Mth = np.hstack([np.cos(ph)*self.M1[:,:self.d1] + np.sin(ph)*self.M2[:,:self.d1],self.M2[:,self.d1:]]);
                    ath = np.hstack([np.cos(ph)*self.c1[:self.d1] + np.sin(ph)*self.c2[:self.d1],self.c2[self.d1:]]);
                    Nth = np.hstack([np.zeros([self.dd,self.d1]),self.M1[:,self.d1:]]);
                    bth = np.hstack([np.zeros([self.d1]),self.c1[self.d1:]])

                    Mthi = mat.inv(Mth);
                    uth = -ath@Mthi;
                    xxith = -bth@Mthi + ath@Mthi@Nth@Mthi;

                    # Hainv = mat.inv(Ha)
                    eigs,VV = mat.eig(Nth@Mthi)
                    VVi = mat.inv(VV)
                    xxiVV = xxith@VV
                    ubase = uth;

                elif tag == 'exp4':
                    # Ha = np.diag([np.cos(t1),np.sin(t1)]);    
                    Mth = np.hstack([self.M1[:,:self.d1],np.cos(ph)*self.M2[:,self.d1:] + np.sin(ph)*self.M1[:,self.d1:]]);
                    ath = np.hstack([self.c1[:self.d1],np.cos(ph)*self.c2[self.d1:] + np.sin(ph)*self.c1[self.d1:]]);
                    Nth = np.hstack([self.M2[:,:self.d1],np.zeros([self.dd,self.d2])]);
                    bth = np.hstack([self.c2[:self.d1],np.zeros(self.d2)])
                    Mthi = mat.inv(Mth);
                    uth = -ath@Mthi;
                    xxith = -bth@Mthi + ath@Mthi@Nth@Mthi;
                    # Hainv = mat.inv(Ha)
                    eigs,VV = mat.eig(Nth@Mthi)
                    VVi = mat.inv(VV);
                    xxiVV = xxith@VV
                    ubase = uth;

                mask1 = np.real(eigs)<0.0
                mask2 = np.abs(np.imag(eigs))<0.000001
                mask = mask1 & mask2; #np.real(eigs)<0.000 & np.abs(np.imag(eigs))<0.0000001;
                EXP['eigs'][ph] = eigs.copy();
                sortinds = np.argsort(np.real(eigs))[::-1];
                sorteigs = eigs[sortinds];
                sortV = VV[:,sortinds];
                sortVi = VVi[sortinds]
                sort_xxiVV = xxiVV[sortinds];

                EXP['neigs'][ph] = eigs[mask]

                tbnds_full = [0.];
                tbreaks = [0.];
                probts = [];
                probpts = [];
                probVis = {};
                probVs = {};
                realparts = np.array([]);


                # coeffs = np.zeros([len(t4s),self.dd],dtype='complex');
                # if tag == 'exp1' or tag == 'exp2':
                #     for j,lam in enumerate(eigs):
                #         coeffs[:,j] = np.array([xxiVV[j]/((lam/t4)+ 1) for t4 in t4s],dtype='complex');
                # elif tag == 'exp3' or tag == 'exp4':
                #     for j,lam in enumerate(eigs):
                #         coeffs[:,j] = np.array([xxiVV[j]*(1./(1./t4 + lam)) for t4 in t4s],dtype='complex');


                track_subspace = False;

                ### FUCK - WHO WROTE THIS SHIT
                found_inds = {};
                conj_found = {}
                for j,eig in enumerate(sorteigs):
                    if  (np.real(eig) < 0.0) and (np.abs(np.imag(eig)) < 0.000001): #and (magimagpart < 0.00001): #(
                        t = np.abs(np.real(eig));
                        tt = np.round(t,6)
                        already_found = False;
                        if tt in probVis: already_found = True;

                        if not(already_found):
                            #realparts = np.hstack([realparts,np.real(eig)])
                            if t < tcutoff: tbreaks.append(t);
                            tbnds_full.append(t)
                            theta1 = np.arctan2(t*np.cos(ph),1)
                            if tag == 'exp1': theta2 = np.arctan2(t*np.sin(ph),1);
                            if tag == 'exp2': theta2 = np.arctan2(1,t*np.sin(ph));
                            probts.append(tt);
                            probpts.append([theta1,theta2]);
                            if not(track_subspace):
                                inds = np.isclose(np.real(sorteigs),np.real(eig))
                                vv = sort_xxiVV[inds]@sortVi[inds]
                                vv = scale_evecs*(1./mat.norm(vv))*vv;
                                probVis[tt] = [vv];

                        #### NOT USED?? ##### 
                        if track_subspace:
                            # magimagpart = mat.norm(sortVi[j])
                            vv = sortVi[j]; vv = (1./mat.norm(vv))*vv;
                            magimagpart = mat.norm(np.imag(vv))
                            if magimagpart > 0.00001: #00000001:
                                remaineigs = sorteigs[j+1:]
                                inds = np.where(np.isclose(remaineigs,np.conj(eig)))[0]
                                
                                if len(inds) > 0:
                                    ### basically only these...
                                    j2 = j+1+inds[0];
                                    # vv1 = np.imag(sortVi[j]);
                                    # vv2 = np.imag(sortVi[j2]);
                                    # print(mat.norm(vv1+vv2))


                                    if True: #track_subspace:
                                        vv1 = np.real(sortVi[j])
                                        vv2 = np.imag(sortVi[j])
                                        WW = np.vstack([vv1,vv2])
                                        WW = smat.sqrtm(mat.inv(WW@WW.T))@WW
                                        vv1 = WW[0]; vv2 = WW[1];
                                        vv = vv1; #sortVi[j] + sortVi[j2];
                                    else:
                                        # print('glarg')
                                        vv1 = sortVi[j];
                                        vv2 = sortVi[j2];
                                        coeff1 = sort_xxiVV[j]
                                        coeff2 = sort_xxiVV[j2]
                                        vv = np.real(vv1*coeff1 + vv2*coeff2);
                                        vv2 = np.ones(len(vv))
                                        # print(vv2)
                                        # vv2 = vv.copy();
                                        # vv = vv1*coeff1 + vv2*coeff2;

                                    # vv2 = sortVi[j2];
                                    # print(sorteigs[j])
                                    # print(sorteigs[j+1])
                                    # MM  = np.vstack([np.real(vv1),np.imag(vv1),np.real(vv2),np.imag(vv2)])
                                    # _,sigs,_ = smat.svd(MM);
                                    # print(sigs)
                                    # print(mat.norm(np.imag(vv1)+np.imag(vv2)))

                                    vv = (1./mat.norm(vv))*vv
                                    conj_found[np.conj(eig)] = vv2;
                                else: 
                                    if eig in conj_found:
                                        vv = conj_found[eig];
                                    else:
                                        ### NONE OF THESE...
                                        vv = sortVi[j]
                                        vv = remove_phase(vv)
                                        vv = (1./mat.norm(vv))*vv
                            else:
                                vv = (1./mat.norm(vv))*vv

                            # print('blurb: ',mat.norm(np.imag(vv)))

                            # badpart = mat.norm(np.imag(vv))
                            # if badpart > 0.1: print('warning bad part found.')
                            vv = scale_evecs*np.real(vv)

                            if not(tt in probVis): probVis[tt] = [vv]
                            else: probVis[tt].append(vv);

                            
                            # probVis.append(vv) #np.real(sortVi[j]))
                probts = np.array(probts);
                probpts = np.array(probpts);
                # probVis = np.array(probVis);


                # if tcutoff > tbreaks[-1]: tbreaks.append(tcutoff);
                # else:
                tbreaks.append(tbreaks[-1]+1000);

                tbnds_full.append(self.maxt)
                EXP['probts'][ph] = probts.copy()                
                EXP['probpts'][ph] = probpts.copy()
                EXP['probVis'][ph] = probVis.copy()
                EXP['tbnds'] = tbnds;
                EXP['tbnds_full'] = tbnds_full

                EXP['Vs'][ph] = VV;
                EXP['Vis'][ph] = VVi
                EXP['Pphs'][ph] = np.real(VV@VV.conj().T)


                if True: 
                    ##### EXTRA ADDED TO COLOR COEFFICIENTS - ONLY FOR THAT, COULD BE DONE BETTER
                    ##### EXTRA ADDED TO COLOR COEFFICIENTS - ONLY FOR THAT, COULD BE DONE BETTER
                    temp = {}
                    temp['une'] = self.une;
                    temp['una'] = self.una;
                    temp['eigs'] = EXP['eigs'][ph];
                    temp['V'] = VV;
                    temp['Vi'] = VVi;


                    LAMnew,VVnew,other = eigDeComplexify(temp['eigs'],VV,ins={'realeps':0.0001});
                    VVinew = mat.inv(VVnew);
                    mostlyreal = other['mostlyreal']
                    reigs = np.diag(LAMnew);
                    reigs = reigs[mostlyreal]
                    VVnew = VVnew[:,mostlyreal]
                    VVinew = VVinew[mostlyreal]
                    inds = np.argsort(reigs);
                    # inds = inds[::-1];
                    reigs = reigs[inds];
                    rVV = VVnew[:,inds];
                    rVVi = VVinew[inds];

                    temp['reigs'] = reigs
                    temp['rVV'] = rVV
                    temp['rVVi'] = rVVi
                    EXP['coeffs_color_data'][ph] = temp
                    ##################################################################################


                EXP['coeffs'][ph] = {}
                EXP['trajs'][ph] = {}
                EXP['coeffs_full'][ph] = {}
                EXP['trajs_full'][ph] = {}

                # if tag == 'exp3' or tag == 'exp4':
                EXP['Mths'][ph] = Mth
                EXP['Nths'][ph] = Nth
                EXP['aths'][ph] = ath
                EXP['bths'][ph] = bth
                EXP['Mthis'][ph] = Mthi
                EXP['uths'][ph] = uth
                EXP['xxiths'][ph] = xxith


                tbnds = EXP['tbnds'];
                tbuff = 0.01;

                intervals = [tbnds];
                for k,_ in enumerate(tbreaks):
                    if k<len(tbreaks)-1: intervals.append((tbreaks[k],tbreaks[k+1]));
                # intervals.insert(0,tbnds);

                for k,interval in enumerate(intervals):
                    t1 = interval[0]; t2 = interval[1]

                    bunch1 = True; bunch2 = True;
                    if k==0: bunch1 = False;
                    if k==len(tbnds)-2: bunch2 = False;

                    if k==0:
                        bparams = bunchparams0;
                    else:
                        bparams = bunchparams;
                        bparams['bunch1'] = bunch1;
                        bparams['bunch2'] = bunch2;


                    t4s = []
                    if show_safe:
                        mask1 = t1 + safe_buff < EXP['ts']; mask2 = EXP['ts'] < t2-safe_buff;
                        t4s = EXP['ts'][mask1 & mask2];
                    if show_sing:

                        ### NEW CODE 
                        scale4 = np.min([np.cos(ph),np.sin(ph)]);
                        thbnds = ts2ths(interval,scale4,version='tan');
                        th1 = thbnds[0]; th2 = thbnds[1];

                        nt4 = bparams['nt']
                        nth1 = int((nt4/2)*(th2 - th1)); nth2 = nth1;
                        thtemp = edgebunch_linspace(th1,th2,nth1=nth1,nth2=nth2,mid=0.5,params = bparams);
                        t4s = ths2ts(thtemp,scale4,version='tan');

                        #th4s = ts2ths(t4s,scale,version='tan');

                        if tag == 'exp1': th4s = np.array([np.array([np.arctan2(t4*np.cos(ph),1),np.arctan2(t4*np.sin(ph),1)]) for t4 in t4s]);
                        if tag == 'exp2': th4s = np.array([np.array([np.arctan2(t4*np.cos(ph),1),np.arctan2(1,t4*np.sin(ph))]) for t4 in t4s]);
                        if tag == 'exp3': th4s = np.array([np.array([ph,np.arctan2(t4,1)]) for t4 in t4s])
                        if tag == 'exp4': th4s = np.array([np.array([np.arctan2(t4,1),ph]) for t4 in t4s])

                        #### OLD CODE
                        # t4s_pre = atan_linspace(t1+self.svo_sing_buff,t1+self.svo_safe_buff,self.svo_singnum);#,extra='inverse')
                        # t4s_post = atan_linspace(t2-self.svo_safe_buff,t2-self.svo_sing_buff,self.svo_singnum);#,extra='inverse')
                        # t4s = np.hstack([t4s_pre,t4s,t4s_post])

                    coeffs = np.zeros([len(t4s),self.dd],dtype='complex');

                    if tag == 'exp1' or tag == 'exp2':
                        for j,lam in enumerate(eigs):
                            coeffs[:,j] = np.array([xxiVV[j]/((lam/t4)+ 1) for t4 in t4s],dtype='complex');
                    elif tag == 'exp3' or tag == 'exp4':
                        for j,lam in enumerate(eigs):
                            coeffs[:,j] = np.array([xxiVV[j]*(1./(1./t4 + lam)) for t4 in t4s],dtype='complex');

                    if k == 0:
                        EXP['coeffs_full'][ph][interval] = coeffs.copy();
                        EXP['trajs_full'][ph][interval] = {}
                        EXP['trajs_full'][ph][interval]['ctrls'] = np.real(ubase + coeffs@EXP['Vis'][ph]);
                        EXP['trajs_full'][ph][interval]['ths'] = th4s; #np.real(ubase + coeffs@EXP['Vis'][ph]);
                        EXP['trajs_full'][ph][interval]['ts'] = t4s;
                        EXP['computed_curves_full']['ucurves'][(ph,interval)] = EXP['trajs_full'][ph][interval]['ctrls']
                        EXP['computed_curves_full']['thcurves'][(ph,interval)] = EXP['trajs_full'][ph][interval]['ths']
                        EXP['computed_curves_full']['tcurves'][(ph,interval)] = EXP['trajs_full'][ph][interval]['ts']

                    else: 
                        EXP['coeffs'][ph][interval] = coeffs.copy();
                        EXP['trajs'][ph][interval] = {}
                        EXP['trajs'][ph][interval]['ctrls'] = np.real(ubase + coeffs@EXP['Vis'][ph]);
                        EXP['trajs'][ph][interval]['ths'] = th4s; #np.real(ubase + coeffs@EXP['Vis'][ph]);
                        EXP['trajs'][ph][interval]['ts'] = t4s;

                        EXP['computed_curves']['ucurves'][(ph,interval)] = EXP['trajs'][ph][interval]['ctrls']
                        EXP['computed_curves']['thcurves'][(ph,interval)] = EXP['trajs'][ph][interval]['ths']
                        EXP['computed_curves']['tcurves'][(ph,interval)] = EXP['trajs'][ph][interval]['ts']

    def integrateCtrls(self,ins={}):

        if self.build_version == 'from_lq':
            AA  = self.LINQUAD.AA;
            BB1 = self.LINQUAD.BB1;
            BB2 = self.LINQUAD.BB2;
            BB = np.hstack([BB1,BB2]);
            HH = self.LINQUAD.HH;
            for tag in self.exptags:
                EXP = self.EXP[tag]
                if tag == 'exp1': ubase = self.une;
                if tag == 'exp2': ubase = self.u1;

                EXP['computed_curves']['xcurves'] = {}
                for ph in EXP['phs']:
                    for k,interval in enumerate(EXP['coeffs'][ph]):
                        uu = EXP['trajs'][ph][interval]['ctrls']
                        if len(uu) > 0:
                            xx = (AA@BB@uu.T).T + HH@self.LINQUAD.x0;
                            EXP['trajs'][ph][interval]['xs'] = xx.copy();
                            EXP['computed_curves']['xcurves'][(ph,interval)] = xx.copy();

                    uu = EXP['probVis'][ph]; xx = {};
                    if len(uu) > 0:
                        for tt in uu:
                            xx[tt] = (AA@BB@(np.array(uu[tt])+ubase).T).T + HH@self.LINQUAD.x0;
                        EXP['probxVis'][ph] = xx.copy();


def diagtan_linspace(vec1,vec2,nph=100,buff=0.01,version='tan'):
    diff = vec2 - vec1;
    ths = vec1 + np.outer(np.linspace(buff,1-buff,nph),diff);
    if version == 'tan':
        phs = np.array([np.arctan2(np.tan(th[0]),np.tan(th[1])) for th in ths]);
    return phs


def atan_linspace(t1,t2,nt,scale=1,buff=0,extra=None):
    bot = np.arctan2(t1*scale,1);
    top = np.arctan2(t2*scale,1);
    
    if extra==None: ts = np.linspace(bot+buff,top-buff,nt)
    elif extra=='inverse': diff = top - buff - (bot+buff); ts = bot+buff+(diff/np.linspace(1,nt,nt))
    elif extra=='inverse_r': diff = top - buff - (bot+buff); ts = bot+buff+(diff/np.linspace(1,nt,nt))[::-1]
    
    ts = np.array([np.tan(t)/scale for t in ts])
    return ts


def eigDeComplexify(eigs,V,Vi=[],ins={}):
    ### warning: assumes complex eigenvalues come in sequential pairs!
    # if len(Vi)==0: Vi = mat.inv(V);
    realeps = 0.00001
    if 'realeps' in ins: realeps = ins['realeps'];

    U2 = (1./np.sqrt(2))*np.array([[1,1j],[1,-1j]]).T
    U2i = mat.inv(U2);
    eiglist = []; Vlist = [];
    i = 0; max_imag = 0;
    mostlyreal = [];
    while i < len(eigs):
        eig = eigs[i];
        if eig.imag == 0:
            eiglist.append(np.real(eig));
            Vlist.append(V[:,i]);
            mostlyreal.append(True)
            i = i+1;
        else: 
            eig1 = eigs[i]; eig2 = eigs[i+1]
            max_imag = np.max([max_imag,np.abs(eig1.imag)]);
            D2 = np.array([[eig1,0.],[0.,eig2]]);
            eiglist.append(np.real(U2@D2@U2i));
            Vlist.append(np.real((V[:,[i,i+1]]@U2i).T))
            if np.abs(eig1.imag) > realeps: mostlyreal = mostlyreal + [False,False];
            else: mostlyreal = mostlyreal + [True,True];
            i = i+2;

    mostlyreal = np.array(mostlyreal)
    other = {}
    other['eigs'] = eiglist;
    other['evecs'] = Vlist;
    other['maximag'] = max_imag;
    other['mostlyreal'] = mostlyreal;
    VV = np.vstack(V).T;
    LAM = smat.block_diag(*eigs);
    return LAM,VV,other
    # return neweigs,newV,out    

def ts2ths(ts,scale,version='tan'):
    # th = atan(t*scale)
    # th = acot(t*scale)
    ths = ts;
    if version == 'tan': ths = np.array([np.arctan2(t*scale,1.) for t in ts]);
    if version == 'cot': ths = np.array([np.arctan2(1.,t*scale) for t in ts]);
    return ths


def edgebunch_linspace(th1,th2,nth1=10,nth2=10,mid=0.5,params = {}):
    bunch1 = False;
    bunch2 = False;
    maxpreth1 = 2;
    maxpreth2 = 2;
    if 'bunch1' in params: bunch1 = params['bunch1'];
    if 'bunch2' in params: bunch2 = params['bunch2'];
    if 'maxpreth1' in params: maxpreth1 = params['maxpreth1'];
    if 'maxpreth2' in params: maxpreth2 = params['maxpreth2'];
    
    thmid = mid*th1 + (1-mid)*th2;
    dth1 = (1-mid)*np.abs(th1-th2); dth2 = mid*np.abs(th1-th2);

    half1 = np.linspace(0,dth1,nth1);
    half2 = np.linspace(0,dth2,nth2);

    if bunch1: half1 = dth1*np.tanh(np.linspace(0,maxpreth1,nth1));
    if bunch2: half2 = dth2*np.tanh(np.linspace(0,maxpreth2,nth2));

    half1 = half1[1:]
    half1 = -half1[::-1] + thmid;
    half2 = half2 + thmid;
    ths = np.hstack([half1,half2]);
    return ths    


def ths2ts(ths,scale,version='tan'):
    # t*scale = tan(th)
    # t*scale = cot(th)
    ts = ths; 
    if version == 'tan': ts = np.array([np.tan(th)/scale for th in ths])
    if version == 'cot': ts = np.array([(1./np.tan(th))/scale for th in ths])
    return ts




class POLYTOPE:
    def __init__(self,ins={}):
        self.nodes = [];
        if 'nodes' in ins: self.nodes = ins['nodes'];
        if isinstance(self.nodes,list): self.nodes = np.array(self.nodes);
        self.faces = [];
        if 'faces' in ins: self.faces = ins['faces'];
        self.edges = [];
        #### ONLY WORKS FOR 2D
        if len(self.faces)>0:
            for face in self.faces:
                for i,node in enumerate(face):
                    node1 = node; 
                    if i<len(face)-1: node2 = face[i+1];
                    elif i==len(face)-1: node2 = face[0]
                    self.edges.append(np.array([node1,node2]))

        elif len(self.nodes)>0:
            for i,node in enumerate(self.nodes):
                node1 = node; 
                if i<len(self.nodes)-1: node2 = self.nodes[i+1];
                elif i==len(self.nodes)-1: node2 = self.nodes[0];
                self.edges.append(np.array([node1,node2]))
        if 'edges' in ins: self.edges = ins['edges'];
        self.nx = len(self.nodes[0])

    
    def ptsInside(self,pts):
        poly = self.nodes;
        filtered_pts,inds = ptsInNHull(pts,poly.T)
        return filtered_pts,inds

    def intersectSegPoly(self,seg,eps=-0.0001):
        out = []
        if self.nx == 2: out = intersectSegPoly2D(seg,self.nodes,eps=eps)
        if self.nx == 3: out = intersectSegPoly3D(seg,self.faces,eps=eps)
        return out

    def intersectPlanePoly(self,plane,eps=-0.0001):
        out = [];
        if self.nx==3: out = intersectPlanePoly3D(plane,self.edges,sort=True,eps=eps);
        return out


class AXIS:
    def __init__(self,ins={}):
        self.PENS = {};
        self.PENS['main'] = PEN()
        if 'pen_main' in ins: self.PENS['main'] = ins['pen_main'];
        if 'pens' in ins: self.PENS = {**self.PENS,**ins['pens']};
        self.nx = 2;
        if 'nx' in ins: self.nx = ins['nx'];
        self.AXES = DEFAULT_AXES.copy()
        if 'AXES' in ins: self.AXES = ins['AXES'];
        self.AXES = self.AXES[:self.nx];
        self.loc = np.zeros(2);
        self.scale = 1.
        if 'loc' in ins: self.loc = ins['loc'];
        if 'ascale' in ins: self.scale = ins['ascale'];
        self.negscale = 0.2 * self.scale;
        if 'negscale' in ins: self.negscale = ins['negscale'];

        self.headwid = 0.05
        self.negheadwid=0;
        if 'headwid' in ins: self.headwid = ins['headwid'];
        if 'negheadwid' in ins: self.negheadwid = ins['negheadwid']
        self.normalize = False; 
        if 'normalize' in ins: self.normalize = ins['normalize']

    def draw(self,ax):
        draw_axes(ax,self.loc,self.AXES,normalize=self.normalize,head_width=self.headwid,scale=self.scale,params={'pen':self.PENS['main']})
        draw_axes(ax,self.loc,-self.AXES,normalize=self.normalize,head_width=self.negheadwid,scale=self.negscale,params={'pen':self.PENS['main']})





def draw_axes(ax, origin, AXES, normalize=False, color='black', alpha=0.5, linewidth=1.0, linestyle='-',
              scale=1.0, negscale=None,show=True,zorder=10,head_width=0.1,params={}):


    if 'normalize' in  params: normalize = params['normalize'];
    if 'edgecolor' in params: color = params['edgecolor'];
    elif 'color' in params: color = params['color'];
    if 'edgealpha' in params: alpha = params['edgealpha'];
    elif 'alpha' in params: alpha = params['alpha'];

    if 'linewidth' in params: linewidth = params['linewidth'];
    if 'linestyle' in params: linestyle = params['linestyle'];
    if 'scale' in params: scale = params['scale'];
    if 'negscale' in params: scale = params['negscale'];
    if 'show' in params: show = params['show'];
    if 'zorder' in params: zorder = params['zorder'];
    if 'headwidth' in params: head_width = params['headwidth'];


    #### PEN UPDATES...
    if 'pen' in params:
        pen = params['pen'];
        facecolor = pen.frgba[:3];
        # edgecolor = pen.ergba[:3];
        color = pen.ergba[:3];
        alpha = pen.ergba[3];
        edge_alpha = pen.ergba[3];
        linestyle = pen.ls
        linewidth = pen.lw

    if not(isinstance(scale, list)): scale = list(scale*np.ones(len(AXES)));
    else: scale = list(scale);

    if show:
        for k in range(np.size(AXES, 0)):
            axis = AXES[k]; negaxis = -axis;
            if normalize:
                axis = 1.0 / mat.norm(axis) * axis
            #axis = scale * axis
            ax.arrow((origin[0]), (origin[1]), (scale[k]*axis[0]), (scale[k]*axis[1]),
                      color=color, linewidth=linewidth, linestyle=linestyle,alpha=alpha,
                      length_includes_head=True, head_width=head_width,zorder=zorder)   




class PTS:
    def __init__(self,ins={}):

        self.id = None
        if 'id' in ins:  self.id = ins['id']
        self.pts = {};
        if 'pts' in ins: self.pts = ins['pts'];
        self.PENS = {};
        self.PENS['main'] = PEN()
        if 'pen_main' in ins: self.PENS['main'] = ins['pen_main'];
        if 'pens' in ins: self.PENS = {**self.PENS,**ins['pens']};

        self.nx = 2;
        if 'nx' in ins: self.nx = ins['nx'];
        self.AXES = DEFAULT_AXES.copy()
        if 'AXES' in ins: self.AXES = ins['AXES'];
        self.AXES = self.AXES[:self.nx];
        self.loc = np.zeros(2);
        self.scale = 1.
        if 'loc' in ins: self.loc = ins['loc'];
        if 'scale' in ins: self.scale = ins['scale'];


        self.LEGS = {};
        if 'legs' in ins: self.LEGS = ins['legs']

        self.curveColorFunc = None;
        self.curveColorInputs = None;
        self.useColorFunc = False; 
        if 'curveColorFunc' in ins: self.curveColorFunc = ins['curveColorFunc']; self.useColorFunc = True; 
        if 'curveColorInputs' in ins: self.curveColorInputs = ins['curveColorInputs'];        

        self.rad = 0.1;
        if 'rad' in ins: self.rad = ins['rad'];

        self.bndpolytopes = {};
        if 'bnds' in ins: self.addBndPolytopes(ins['bnds']);

        self.PTS = {};
        for i,pt in enumerate(self.pts):
            self.PTS[pt] = {};
            self.PTS[pt]['pt'] = self.pts[pt]
            self.PTS[pt]['filtered'] = {};
            # if 'pens' in ins:
            #     if curve in ins['pens']: self.CURVE[curve]['pen'] = ins['pens'][curve];
            #     else: self.CURVE[curve]['pen'] = self.PENS['main']

        self.drawtype = 'patch'
        if 'drawtype' in ins: self.drawtype = ins['drawtype']
        if 'dtyp' in ins: self.drawtype = ins['dtyp']

    def addBndPolytopes(self,polys):
        self.bndpolytopes = {**self.bndpolytopes,**polys}


    def addPts(self,pts={},overwrite=True,ins={}):

        for i,pt in enumerate(pts):
            if overwrite:
                self.pts[pt] = pts[pt]
                self.PTS[pt] = {};
                self.PTS[pt]['pt'] = self.pts[pt]
                self.PTS[pt]['filtered'] = {};
                if 'pens' in ins:
                    if pt in ins['pens']:
                        self.PENS[pt] = ins['pens'][pt]
                
            else: 
                if not(pt in self.pts):
                    self.pts[pt] = pts[pt]
                    self.PTS[pt] = {};
                    self.PTS[pt]['pt'] = self.pts[pt]
                    self.PTS[pt]['filtered'] = {};
                    if 'pens' in ins:
                        if pt in ins['pens']:
                            self.PENS[pt] = ins['pens'][pt]


    def filterPts(self,polytag,pt_tags=[],params={}):

        filterstyle='remove';
        if 'filterstyle' in params: filterstyle = params['filterstyle'];
        if 'nbndpts' in params: nbndpts = params['nbndpts'];
        if not(polytag==None):
            POLY = self.bndpolytopes[polytag];
            if len(pt_tags)==0: pt_tags = list(self.PTS)
            for tag in pt_tags:
                pt = self.PTS[tag]['pt'];
                if filterstyle == 'remove':
                    filtered_pts,_ = POLY.ptsInside([pt])
                    if len(filtered_pts)>0: self.PTS[tag]['filtered'][polytag] = filtered_pts[0];
                if filterstyle=='boundary':
                    filtered_pts,_ = POLY.ptsInside([pt])
                    if len(filtered_pts)>0: self.PTS[tag]['filtered'][polytag] = filtered_pts[0];
                    else: self.PTS[tag]['filtered'][polytag] = snap2boundary2D(pt,POLY.nodes,npts = nbndpts)

    def drawPts(self,ax,poly=None,pts=[],typ='filtered'):
        if len(pts)==0: pt_tags = list(self.PTS);
        else: pt_tags = pts;
        for i,tag in enumerate(pt_tags):
            self.draw(ax,tag,poly=poly,typ=typ);

    def draw(self,ax,tag,poly=None,typ='filtered',ins={}):
        pt = self.PTS[tag]['pt']
        # if typ == 'pt':
        pen = self.PENS['main'];
        if tag in self.PENS: pen = self.PENS[tag];
        draw_point = True;
        if typ == 'pts':
            pt = self.PTS[tag]['pt'];
        if typ == 'filtered':
            draw_point = False;
            if poly in self.PTS[tag]['filtered']:
                draw_point = True;
                pt = self.PTS[tag]['filtered'][poly];

        scale = self.scale; rad = self.rad; loc = self.loc; 
        if 'scale' in ins: scale = ins['scale'];
        if 'rad' in ins: rad = ins['rad'];
        # if 'loc' in params: loc = params['loc'];
        

        if draw_point: 
            if self.useColorFunc:
                _,_,tempcolors = self.curveColorFunc(tag,ins=self.curveColorInputs);
                if isinstance(tempcolors,list): tempcolor = tempcolors[0];
                else: tempcolor = tempcolors
            else: tempcolor = pen.frgba[:3]
            # pen = PEN(ins={'fc':tempcolor});
            # print(self.useColorFunc)
            if self.drawtype == 'patch':
                params = {};
                params['facecolor'] = tempcolor
                params['edgecolor'] = pen.ergba[:3];
                params['facealpha'] = pen.frgba[3];
                params['edgealpha'] = pen.ergba[3];
                params['linestyle'] = pen.ls;
                params['linewidth'] = pen.lw
                params['edgesolid'] = True;
                params['fill'] = True;
                params['zorder'] = pen.zord

                if 'absloc' in ins: pt = ins['absloc']
                pt2 = self.scale*pt@self.AXES;
                if 'absloc' in ins: pt2 = np.zeros(2); loc = ins['absloc'];
                draw_patch(ax,[pt2],np.eye(2),rad=self.rad,shift=loc,params=params)
            if self.drawtype == 'marker':
                pt2 = self.scale*pt@self.AXES + self.loc
                if 'absloc' in ins: pt2 = ins['absloc'];
                ax.plot(pt2[0],pt2[1],color=pen.ergba[:3],alpha=pen.ergba[3],marker=pen.msty,markersize=pen.msz,zorder=pen.zord)

    def addToLegend(self,ax,ins={}):

        leg = 'basic';
        if 'leg' in ins: leg = ins['leg'];
        tags = [];
        if 'tags' in ins: tags = ins['tags'];
        lloc = np.zeros(2); lsize = 1.;
        if 'loc' in ins: lloc = ins['loc'];
        if 'sz' in ins: lsize = ins['sz'];
        if isinstance(lloc,list): lloc = np.array(lloc)

        info = ins['info'];
        ### defaults
        rloc = np.zeros(2);
        rsz = 1.;
        rann = 'tag';
        rann_loc = np.zeros(2);
        rann_sz = 1.;
        ### dicts
        rlocs = info['locs'];
        tags = list(rlocs);
        if 'tags' in info: tags = info['tags'];
        rsizs = info['sizes'];
        anns = info['anns'];
        alocs = info['alocs'];
        asizs = info['asizes'];
        for tag in tags:
            ####
            rloc = rlocs[tag];
            arloc = alocs[tag];
            if isinstance(rloc,list): rloc = np.array(rloc);
            if isinstance(arloc,list): arloc = np.array(arloc);

            loc = lloc + lsize * rloc;
            sz = lsize * rsizs[tag];
            rad = lsize * rsizs[tag] * self.rad;
            self.draw(ax,tag,typ='pts',ins={'absloc':loc});
            #'scale':sz,'absloc':loc,'rad':rad});
            ####
            aloc = lloc + lsize * arloc;
            asz = lsize * asizs[tag];
            writeAnn(ax,text=anns[tag],loc=aloc,sz=asz)
        #'leg':self.name,'loc':self.loc,'scale':self.scale})                


class CURVES:
    def __init__(self,ins={}):
        self.curves = {};
        self.tcurves = {};
        # self.curve_pen_range_tags = {};
        if 'curves' in ins: self.curves = ins['curves'];
        if 'tcurves' in ins: self.tcurves = ins['tcurves'];

        # if 'curve_pen_range_tags' in ins: self.curve_pen_range_tags = ins['curve_pen_range_tags'];

        self.PENS = {};
        self.PENS['main'] = PEN()
        if 'pen_main' in ins: self.PENS['main'] = ins['pen_main'];

        self.curveColorFunc = None;
        self.curveColorInputs = None;
        self.useColorFunc = False; 
        if 'curveColorFunc' in ins: self.curveColorFunc = ins['curveColorFunc']; self.useColorFunc = True; 
        if 'curveColorInputs' in ins: self.curveColorInputs = ins['curveColorInputs'];

        self.curveColorGrad = None;
        self.useColorGrad = False;
        if 'curveColorGrad' in ins: self.curveColorGrad = ins['curveColorGrad']
        if 'useColorGrad' in ins: self.useColorGrad = ins['useColorGrad'];

        self.LEGS = {};
        if 'legs' in ins: self.LEGS = ins['legs']

        self.CURVES = {};
        for i,curve in enumerate(self.curves):
            self.CURVES[curve] = {};
            self.CURVES[curve]['pts'] = self.curves[curve];
            self.CURVES[curve]['ts'] = [];
            self.CURVES[curve]['filtered'] = {};
            self.CURVES[curve]['filtered_ts'] = {};
            self.CURVES[curve]['segmented'] = {};
            self.CURVES[curve]['segmented_ts'] = {};
            # self.CURVES[curve]['pen_range_tag'] = None
            if len(self.tcurves)>0: self.CURVES[curve]['ts'] = self.tcurves[curve];
            # if curve in self.curve_pen_range_tags: self.CURVES[curve]['pen_range_tag'] = self.curve_pen_range_tags[curve];
            if 'pens' in ins:
                if curve in ins['pens']: self.CURVES[curve]['pen'] = ins['pens'][curve];
                else: self.CURVES[curve]['pen'] = self.PENS['main']


        # self.PEN_RANGES = {};
        # if 'pen_ranges' in ins: self.PEN_RANGES = ins['pen_ranges'];

        self.nx = 2; 
        if 'nx' in ins: self.nx = ins['nx'];
        self.AXES = DEFAULT_AXES.copy()
        if 'AXES' in ins: self.AXES = ins['AXES'];
        self.AXES = self.AXES[:self.nx];
        self.loc = np.zeros(2);
        self.scale = 1.
        if 'loc' in ins: self.loc = ins['loc'];
        if 'scale' in ins: self.scale = ins['scale'];

        self.use_shadow = False;
        if 'use_shadow' in ins: self.use_shadow = ins['use_shadow']

        self.bndpolytopes = {};
        if 'bnds' in ins: self.addBndPolytopes(ins['bnds']);

        self.PENDICT = {};
        if 'pendict' in ins: self.PENDICT = ins['pendict'];

    def addBndPolytopes(self,polys):
        self.bndpolytopes = {**self.bndpolytopes,**polys}

    def addCurves(self,curves={},tcurves={}):
        for i,curve in enumerate(curves):
            self.curves[curve] = curves[curve]
            if curve in tcurves: self.tcurves[curve] = tcurves[curve];
            self.CURVES[curve] = {};
            self.CURVES[curve]['pts'] = curves[curve];
            self.CURVES[curve]['ts'] = [];
            self.CURVES[curve]['filtered'] = {};
            self.CURVES[curve]['filtered_ts'] = {};
            self.CURVES[curve]['segmented'] = {};
            self.CURVES[curve]['segmented_ts'] = {};
            if curve in tcurves: self.CURVES[curve]['ts'] = tcurves[curve];
            if 'pens' in ins:
                if curve in ins['pens']: self.CURVES[curve]['pen'] = ins['pens'][curve];
                else: self.CURVES[curve]['pen'] = self.PENS['main']

    ###########################################################################
    ##### FILTERS - FILTERS - FILTERS - FILTERS - FILTERS - FILTERS - FILTERS - 
    ###########################################################################

    def filterCurve(self,curve,polytag):
        pts = self.CURVES[curve]['pts'];
        ts = self.CURVES[curve]['ts'];
        POLY = self.bndpolytopes[polytag];
        filtered_pts,inds = POLY.ptsInside(pts)
        filtered_ts = ts;
        if len(ts)>0: filtered_ts = ts[inds]
        self.CURVES[curve]['filtered'][polytag] = filtered_pts;
        self.CURVES[curve]['filtered_ts'][polytag] = filtered_ts;
    def segmentCurve(self,curve,polytag,eps=0.01):
        filtered_pts = self.CURVES[curve]['filtered'][polytag]
        filtered_ts = self.CURVES[curve]['filtered_ts'][polytag]
        segmented_pts,segmented_inds = filter_for_gaps(filtered_pts,eps=eps);   
        segmented_ts = [];
        for i,inds in enumerate(segmented_inds):
            if len(filtered_ts)>0: segmented_ts.append(filtered_ts[inds])
            else: segmented_ts.append(filtered_ts);
        self.CURVES[curve]['segmented'][polytag] = segmented_pts;
        self.CURVES[curve]['segmented_ts'][polytag] = segmented_ts;
    def filterAndSegment(self,polytag,curves=[],eps=0.02):
        if len(curves)==0: curves = list(self.CURVES.keys())
        for curve in curves:
            self.filterCurve(curve,polytag);
            self.segmentCurve(curve,polytag,eps=eps)

    ##############################################################################################
    ##### PLOTTING - PLOTTING - PLOTTING - PLOTTING - PLOTTING - PLOTTING - PLOTTING - ###########
    ##############################################################################################


    def drawCurve(self,ax,curve,poly=None,typ='segmented',params={}):
        if typ == 'pts':
            segs = [self.CURVES[curve]['pts']]
            tsegs = [self.CURVES[curve]['ts']]
        if typ == 'filtered':
            segs = [self.CURVES[curve]['filtered'][poly]];
            tsegs = [self.CURVES[curve]['filtered_ts'][poly]];
        if typ == 'segmented':
            segs = self.CURVES[curve]['segmented'][poly];
            tsegs = self.CURVES[curve]['segmented_ts'][poly];

        for i,seg in enumerate(segs):
            ### version 1: defaults
            tseg = [];
            if len(tsegs)>0: tseg = tsegs[i]
            pts = self.scale*seg@self.AXES + self.loc;
            chunks = [pts];

            ### version 2: use color funcs
            if self.useColorFunc:
                inds,_,tempcolors = self.curveColorFunc(curve,ins={**self.curveColorInputs,**{'ts':tseg}}); 
                chunks = [self.scale*seg[ind]@self.AXES + self.loc for ind in inds];

            if self.useColorGrad: ### should be subsumed into useColorFunc
                nchunks = self.curveColorGrad['nchunks']; colorlist = self.curveColorGrad['colorlist'];
                dchunk = int(len(pts)/nchunks);
                chunks = []; tempcolors = [];
                for j in range(nchunks-1): chunks.append(pts[j*dchunk:(j+1)*dchunk+1]);
                chunks.append(pts[(nchunks-1)*dchunk:]);
                for j in range(nchunks): tempcolors.append(colorlist[int((j/nchunks)*len(colorlist))]);

            for j,chunk in enumerate(chunks):
                # pen = self.CURVES[curve]['pen'];
                # color = pen.ergba[:3]; alpha = pen.ergba[3]; lw = pen.lw; ls = pen.ls;
                lw = self.PENS['main'].lw; ls = self.PENS['main'].ls; zord = self.PENS['main'].zord;
                msty = self.PENS['main'].msty; msz = self.PENS['main'].msz; 
                if self.useColorFunc or self.useColorGrad:
                    color = tempcolors[j][:3]; #alpha = self.PENS['main'].ergba[3]
                    alpha = tempcolors[j][3];
                else:
                    pen = self.CURVES[curve]['pen']
                    color = pen.ergba[:3]
                    alpha = pen.ergba[3]
                    lw = pen.lw; ls = pen.ls; zord = pen.zord;
                    msty = pen.msty; msz = pen.msz;
                # else:
                #     # color = [0.,0.,0.]
                #     color = self.PENS['main'].ergba[:3];
                #     alpha = self.PENS['main'].ergba[3];
                    # alpha = 1.                
                if self.use_shadow: ax.plot(chunk[:,0],chunk[:,1],color=[0,0,0],alpha=alpha,lw=1.2*lw,ls=ls,zorder=zord);
                if len(msty) > 0: ax.plot(chunk[:,0],chunk[:,1],color=color,alpha=alpha,lw=lw,ls=ls,marker=msty,markersize=msz,zorder=zord);
                else: ax.plot(chunk[:,0],chunk[:,1],color=color,alpha=alpha,lw=lw,ls=ls,zorder=zord);


    def drawCurves(self,ax,poly=None,curves=[],typ='segmented'):
        if len(curves)==0: curves = list(self.CURVES);
        for i,curve in enumerate(curves):
            self.drawCurve(ax,curve,poly=poly,typ=typ);

    # def drawDiff(self,ax,tag1,tag2,k=None,k2=None,version='line',ins={}):
    def drawDiff(self,ax,tag1,tag2,poly=None,version='line',ins={}):

        typ = 'pts';
        if 'typ' in ins: typ = ins['typ'];
        # if k==None: k = 0;
        # if k2==None: k2 = k;

        if typ == 'pts':
            segs1 = [self.CURVES[tag1]['pts']]
            segs2 = [self.CURVES[tag2]['pts']]
            tsegs1 = [self.CURVES[tag1]['ts']]
            tsegs2 = [self.CURVES[tag2]['ts']]

        if False:  ## DOESN'T WORK YET
            if typ == 'filtered':
                segs = [self.CURVES[curve]['filtered'][poly]];
                tsegs = [self.CURVES[curve]['filtered_ts'][poly]];
        if False: ## PROBABLY USELESS
            if typ == 'segmented':
                segs = self.CURVES[curve]['segmented'][poly];
                tsegs = self.CURVES[curve]['segmented_ts'][poly];

        # for i,seg in enumerate(segs):
        #     ### version 1: defaults
        #     tseg = [];
        #     if len(tsegs)>0: tseg = tsegs[i]
        #     pts = self.scale*seg@self.AXES + self.loc;
        #     chunks = [pts];

        ####################


        # pts1 = self.pts[tag1][k];
        # pts2 = self.pts[tag2][k];

        pen = self.PENS['main']
        if version == 'arrow1to2' or version == 'line':
            if (tag1,tag2) in self.PENS: pen = self.PENS[(tag1,tag2)]
        if version == 'arrow2to1' or version == 'line':
            if (tag2,tag1) in self.PENS: pen = self.PENS[(tag2,tag1)]
        if 'pen' in ins: pen = ins['pen']

        fc = pen.fc; ec = pen.ec;
        fa = pen.fa; ea = pen.ea;
        elw = pen.elw; flw = pen.flw;
        ls = pen.ls; msty = pen.msty; msz = pen.msz
        params2 = {'width':1,'head_width':1.2}
        params2['linewidth'] = elw;
        params2['linestyle'] = ls;
        params2['facecolor'] = fc;
        params2['edgecolor'] = ec;        
        if 'width' in ins: params2['width'] = ins['width'];
        if 'wid' in ins: params2['width'] = ins['wid'];
        if 'headwidth' in ins: params2['headwidth'] = ins['headwidth'];
        if 'hwid' in ins: params2['headwidth'] = ins['hwid'];

        for i,seg1 in enumerate(segs1):
            seg2 = segs2[i];
            pts1 = seg1;
            pts2 = seg2;
            for t in range(len(pts1)):
                pt1 = pts1[t];
                pt2 = pts2[t];
                pt1 = self.scale*pt1@self.AXES + self.loc;
                pt2 = self.scale*pt2@self.AXES + self.loc;
                if version == 'line':
                    ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],color=ec,alpha=ea,linestyle=ls,linewidth=elw,marker=msty,markersize=msz);
                if version == 'arrow1to2':
                    draw_arrow(ax,pt1,pt2,np.eye(2),params=params2)
                if version == 'arrow2to1':
                    draw_arrow(ax,pt2,pt1,np.eye(2),params=params2)

    def drawDiffs(self,ax,tag1,tag2s=[],poly=None,version='line',ins={}):
        # if len(tags)==0: tags = list(self.pts);
        # typ = 'pts';
        # if 'typ' in ins: typ = ins['typ'];
        # if k1 == None: k1 = 0;
        # if len(ks)==0: ks = [0];
        # if len(k2s)==0: k2s = [0];
        for tag2 in tag2s:
            # for k2 in k2s:
            #     if k2 <= len(self.pts[tag2])-1:
            self.drawDiff(ax,tag1,tag2,poly=poly,version=version,ins=ins);


    # #################################################################################            

    # def plotTraj(self,ax,tag,k=None,ins={}):

    #     typ = 'pts';
    #     if 'typ' in ins: typ = ins['typ'];

    #     if k==None: k = 0;
    #     pts = self.pts[tag][k];
    #     pts = self.scale*pts@self.AXES + self.loc
    #     pen = self.PENS['main']
    #     if tag in self.PENS: pen = self.PENS[tag];
    #     fc = pen.fc; ec = pen.ec;
    #     fa = pen.fa; ea = pen.ea;
    #     elw = pen.elw; flw = pen.flw;
    #     ls = pen.ls; msty = pen.msty; msz = pen.msz
    #     ax.plot(pts[:,0],pts[:,1],color=ec,alpha=ea,linestyle=ls,linewidth=elw,marker=msty,markersize=msz);

    # def plotTrajs(self,ax,tags=[],ks=[],ins={}):
    #     typ = 'pts';
    #     if 'typ' in ins: typ = ins['typ'];
    #     if len(tags)==0: tags = list(self.pts);
    #     if len(ks)==0: ks = [0];
    #     for tag in tags:
    #         for k in ks:
    #             if k <= len(self.pts[tag])-1:
    #                 self.plotTraj(ax,tag,k=k,ins=ins);

    #################################################################################
    #################################################################################


    #################################################################################
    #################################################################################


    def addToLegend(self,ax,ins={}):
        leg = 'basic';
        if 'leg' in ins: leg = ins['leg'];
        tags = ['main'];
        if 'tags' in ins: tags = ins['tags'];
        lloc = np.zeros(2); lsize = 1.;
        if 'loc' in ins: lloc = ins['loc'];
        if 'sz' in ins: lsize = ins['sz'];
        if isinstance(lloc,list): lloc = np.array(lloc)

        info = ins['info'];
        ### defaults
        rloc = np.zeros(2);
        rsz = 1.;
        rann = 'tag';
        rann_loc = np.zeros(2);
        rann_sz = 1.;
        ### dicts
        rlocs = info['locs'];

        tags = list(rlocs);
        if 'tags' in info: tags = info['tags'];
        rsizs = info['sizes'];
        anns = info['anns'];
        alocs = info['alocs'];
        asizs = info['asizes'];

        typs = {tag:'hseg' for tag in tags};
        if 'typs' in info:
            for tag in tags:
                if tag in info['typs']: typs[tag] = info['typs'][tag]

        for tag in tags:
            typ = typs[tag]
            ####
            rloc = rlocs[tag];
            arloc = alocs[tag];
            if isinstance(rloc,list): rloc = np.array(rloc);
            if isinstance(arloc,list): arloc = np.array(arloc);
            loc = lloc + lsize * rloc;
            sz = lsize * rsizs[tag];

            if typ == 'hseg':
                params = {'pen':self.PENS['main']};
                shape = np.array([[-1.,0.],[1.,0.]])
                pts = sz*shape;
                draw_line(ax,pts[0],pts[1],np.eye(2),shift=loc,params=params);
            if typ == 'vseg':
                params = {'pen':self.PENS['main']};
                shape = np.array([[0.,-1.],[0.,1.]])
                pts = sz*shape;
                draw_line(ax,pts[0],pts[1],np.eye(2),shift=loc,params=params);
                
            #'scale':sz,'absloc':loc,'rad':rad});
            ####
            aloc = lloc + lsize * arloc;
            asz = lsize * asizs[tag];
            writeAnn(ax,text=anns[tag],loc=aloc,sz=asz)
        #'leg':self.name,'loc':self.loc,'scale':self.scale})                     


def ptInNHull(x,A,all_combs = False,eps=-0.0001):
    shp = np.shape(A);
    combs = [list(ind) for ind in itertools.combinations(list(range(shp[1])),shp[0]+1)]
    zs = [];
    cols = [];
    inside = False;
    for i,comb in enumerate(combs):
        out = ptInNTet(x,A[:,comb],eps=eps)
        if out['in'] == True:
            inside = True;
            cols.append(comb)
            zs.append(out['comb'])
            if not(all_combs):
                break
    zs = np.array(zs);
    return {'in':inside,'combs':zs,'inds':cols}

def ptsInNHull(pts,hullnD):
    outs = [];
    inds = [];
    for i,pt in enumerate(pts):
        out = ptInNHull(pt,hullnD); #,all_combs=True)
        if out['in']:
            outs.append(pt)
            inds.append(i);
    outs = np.array(outs)
    return outs,inds

########################################
########################################
### CAREFUL PROBLEMATIC FOR HIGH DIM
def ptInNTet(x,A,eps=-0.0001):
    # tries to solve x = Az and sum(z) = 1 for x in Rn and A in R^(n x n+1)
    shp = np.shape(A)
    m = shp[0]; n = shp[1]
    if not(len(x)==m):
        print('x:',x)
        print('A:',A)
        raise Exception("x and A not compatible dimensions.")
    if not(n==m+1): raise Exception("A doesn't have 1 more column than row.")
    xbar = np.hstack([x,1]); 
    M = np.vstack([A,np.ones(n)])
    inside = False
    try:
        z = mat.inv(M)@xbar
        if (z >= eps).all(): inside = True;
    except:
        z = mat.pinv(M)@xbar; inside = False; 
    return {'in':inside,'comb':z}

def filter_for_gaps(curve,eps=0.1):
    new_curves = [];
    current = [];
    new_inds = []; current_inds = []; ###########
    for i in range(len(curve)):
        current.append(curve[i])
        current_inds.append(i); ###########
        if i<(len(curve)-1):
            diff = curve[i+1]-curve[i];
            if mat.norm(diff)>eps:
                current = np.array(current)
                if len(current)>0:
                    new_curves.append(current)
                    new_inds.append(current_inds) ###########
                current = [];
                current_inds = []; ###########
    current = np.array(current);
    if len(current)>0:
        new_curves.append(current);
        new_inds.append(current_inds) ###########
    return new_curves,new_inds
    



parula_cm_data = [[0.2081, 0.1663, 0.5292], [0.2116238095, 0.1897809524, 0.5776761905], 
 [0.212252381, 0.2137714286, 0.6269714286], [0.2081, 0.2386, 0.6770857143], 
 [0.1959047619, 0.2644571429, 0.7279], [0.1707285714, 0.2919380952, 
  0.779247619], [0.1252714286, 0.3242428571, 0.8302714286], 
 [0.0591333333, 0.3598333333, 0.8683333333], [0.0116952381, 0.3875095238, 
  0.8819571429], [0.0059571429, 0.4086142857, 0.8828428571], 
 [0.0165142857, 0.4266, 0.8786333333], [0.032852381, 0.4430428571, 
  0.8719571429], [0.0498142857, 0.4585714286, 0.8640571429], 
 [0.0629333333, 0.4736904762, 0.8554380952], [0.0722666667, 0.4886666667, 
  0.8467], [0.0779428571, 0.5039857143, 0.8383714286], 
 [0.079347619, 0.5200238095, 0.8311809524], [0.0749428571, 0.5375428571, 
  0.8262714286], [0.0640571429, 0.5569857143, 0.8239571429], 
 [0.0487714286, 0.5772238095, 0.8228285714], [0.0343428571, 0.5965809524, 
  0.819852381], [0.0265, 0.6137, 0.8135], [0.0238904762, 0.6286619048, 
  0.8037619048], [0.0230904762, 0.6417857143, 0.7912666667], 
 [0.0227714286, 0.6534857143, 0.7767571429], [0.0266619048, 0.6641952381, 
  0.7607190476], [0.0383714286, 0.6742714286, 0.743552381], 
 [0.0589714286, 0.6837571429, 0.7253857143], 
 [0.0843, 0.6928333333, 0.7061666667], [0.1132952381, 0.7015, 0.6858571429], 
 [0.1452714286, 0.7097571429, 0.6646285714], [0.1801333333, 0.7176571429, 
  0.6424333333], [0.2178285714, 0.7250428571, 0.6192619048], 
 [0.2586428571, 0.7317142857, 0.5954285714], [0.3021714286, 0.7376047619, 
  0.5711857143], [0.3481666667, 0.7424333333, 0.5472666667], 
 [0.3952571429, 0.7459, 0.5244428571], [0.4420095238, 0.7480809524, 
  0.5033142857], [0.4871238095, 0.7490619048, 0.4839761905], 
 [0.5300285714, 0.7491142857, 0.4661142857], [0.5708571429, 0.7485190476, 
  0.4493904762], [0.609852381, 0.7473142857, 0.4336857143], 
 [0.6473, 0.7456, 0.4188], [0.6834190476, 0.7434761905, 0.4044333333], 
 [0.7184095238, 0.7411333333, 0.3904761905], 
 [0.7524857143, 0.7384, 0.3768142857], [0.7858428571, 0.7355666667, 
  0.3632714286], [0.8185047619, 0.7327333333, 0.3497904762], 
 [0.8506571429, 0.7299, 0.3360285714], [0.8824333333, 0.7274333333, 0.3217], 
 [0.9139333333, 0.7257857143, 0.3062761905], [0.9449571429, 0.7261142857, 
  0.2886428571], [0.9738952381, 0.7313952381, 0.266647619], 
 [0.9937714286, 0.7454571429, 0.240347619], [0.9990428571, 0.7653142857, 
  0.2164142857], [0.9955333333, 0.7860571429, 0.196652381], 
 [0.988, 0.8066, 0.1793666667], [0.9788571429, 0.8271428571, 0.1633142857], 
 [0.9697, 0.8481380952, 0.147452381], [0.9625857143, 0.8705142857, 0.1309], 
 [0.9588714286, 0.8949, 0.1132428571], [0.9598238095, 0.9218333333, 
  0.0948380952], [0.9661, 0.9514428571, 0.0755333333], 
 [0.9763, 0.9831, 0.0538]]

parula_colormap = LinearSegmentedColormap.from_list('parula', parula_cm_data)

