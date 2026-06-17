import numpy as np
import numpy.linalg as mat
import scipy as sp
import scipy.linalg as smat
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

import time
import warnings
warnings.filterwarnings('ignore')

from svo_classes import *
from svo_styles import *


#### PENS #### PENS #### PENS #### PENS #### PENS #### PENS #### PENS #### PENS #####
#### PENS #### PENS #### PENS #### PENS #### PENS #### PENS #### PENS #### PENS #####
#### PENS #### PENS #### PENS #### PENS #### PENS #### PENS #### PENS #### PENS #####



directions = [1,0,0]; #,0,0,0]; #,0]; #[2,0,0];
scale = 1.6; negscale= 0.4; #[2.2,1.6]; #0.8; #[0.2,0.2,0.2];
Ap = np.array([[1.25,0.35],[0.35,1.25]]).T; # A = np.hstack([Ap,Ap@U]) #CASE 3



show_version = 'LQ'
fig_version = 'regular'


filename = show_version + '_' + fig_version;



# example_version = '4D'

hcenter = np.array([np.pi/4,np.pi/4]); hrad = 0.6; # INTERSECTS
#hcenter = np.array([2.8*np.pi/8,2.8*np.pi/8]); hrad = 0.35; # BLOWUPS

nh = 80;
hangles = np.linspace(0,2*np.pi,nh);
hh1 = 0; hh2 = 80;
hangles = hangles[hh1:hh2];
th_highlights = [hcenter + hrad*np.array([np.cos(ang),np.sin(ang)]) for ang in hangles]




#figname = 'svocurves_NEtoSW_'+expans;

axisinfo = {}
axisinfo['axis1'] = {'nx':2,'AXES':DEFAULT_AXES,'loc':np.array([0,0]),'scale':1.0};
# axisinfo['axis2'] = {'nx':2,'AXES':DEFAULT_AXES,'loc':np.array([1.35,0.65]),'scale':0.35};
axisinfo['axis2'] = {'nx':2,'AXES':DEFAULT_AXES,'loc':np.array([0.,0.]),'scale':1.0};
axisinfo['axis3'] = {'nx':2,'AXES':DEFAULT_AXES,'loc':np.array([0.,0.]),'scale':0.4};
axisinfo['axis4'] = {'nx':2,'AXES':DEFAULT_AXES,'loc':np.array([0.,0.]),'scale':1.,'ascale':0.6,'negscale':0.4};


if show_version == 'LQ':

    quad_game_version = 'from_lq'


    ftags = {'u1a':'uframe2D','u1b':'uframe2D',
             'u2a':'uframe2D','u2b':'uframe2D',
             'ellipses':'uframe2D',
             'tframe':'tframe',
             'uframe':'uframelq'}

    # axisinfo  = {}
    axisinfo['axis1'] = {'nx':2,'AXES':DEFAULT_AXES,'loc':np.array([0,0]),'scale':1.0};
    # axisinfo['axis2'] = {'nx':2,'AXES':DEFAULT_AXES,'loc':np.array([1.6,0.7]),'scale':0.35};
    axisinfo['axis3'] = {'nx':2,'AXES':DEFAULT_AXES,'loc':np.array([0.0,0.0]),'scale':1.0};
    axisinfo['axis4'] = {'nx':2,'AXES':DEFAULT_AXES,'loc':np.array([0.0,0.0]),
                         'scale':1.0,'ascale':0.4,'negscale':0.4};

    nu1 = 2; nu2 = 2; nu = nu1 + nu2;
    nx1 = 2; nx2 = 2; nx = nx1 + nx2;


    use_ctrl = True;

    if True: #vers == 'basiclq':  ### PAPER VERSION 
        # ns = 10; nt = 50; dt = 0.04;
        ns = 10; nt = 50; dt = 0.04;
        # ns = 10; nt = 10; dt = 0.1;
        theta = np.pi/8; 
        Xsys = np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
        Xisys = mat.inv(Xsys)
        rott = 0.; stab = 0.;#2;
        # rott = 0.4; stab = 0.4;#2;
        Dsys = np.array([[stab,-rott],[rott,stab]]);
        Asys = Xsys@Dsys@Xisys;
        #Bsys = np.array([[1,1]]).T; #np.eye(2)[:,[0]];
        
        #nx = Asys.shape[0]; nu = Bsys.shape[0]; #nu1 = int(nu/2); nu2 = int(nu/2);



        if True: #### SVO2
            zscale1 = 1.2; zscale2 = 1.2;
            zscale1b = 0.0; zscale2b = 0.0;
            prox1 = -0.4; prox2 = -0.4;
            uscale1 = 0.1; uscale2 = 0.1;
            uscale1b = 0.000001; uscale2b = 0.000001;
            fwin = 5; ffac = 2.; use_ffac = False;
            

        paramsIndAgents = {'r1':uscale1*dt,'r1b':uscale1b*dt,'q1':zscale1*dt,'q1b':zscale1b*dt,'p1':-prox1*dt,
                           'r2':uscale2*dt,'r2b':uscale2b*dt,'q2':zscale2*dt,'q2b':zscale2b*dt,'p2':-prox2*dt}

        Q1base = np.eye(nx1); Q2base = np.eye(nx2);
        R1base = np.eye(nu1); R2base = np.eye(nu2);
        R12base = np.eye(nu2); R21base = np.eye(nu1); 
        Qprox_base = np.eye(nx1); #np.diag([1,1])

        Bsys = np.eye(nu1)
        Qsys = zscale1*np.eye(nx1)
        Rsys = uscale1*np.eye(nu1); #diag([10.,1.,10.,1.])
        Psys = smat.solve_continuous_are(Asys, Bsys, Qsys, Rsys);
        Ksys = -mat.inv(Rsys)@Bsys.T@Psys

        Qprox = np.block([[Qprox_base,-Qprox_base],
                          [-Qprox_base,Qprox_base]])

        Q1zero = smat.block_diag(zscale1*np.eye(nx1),zscale1b*np.eye(nx2)); #np.zeros([nx,nx]))
        Q2zero = smat.block_diag(zscale2b*np.eye(nx1),zscale2*np.eye(nx2))
        R1ctrl = uscale1*R1base
        R2ctrl = uscale2*R2base

        R12ctrl = uscale1b*R12base
        R21ctrl = uscale2b*R21base

        # x02player = 1.*np.array([-1.0,0.2342,0.3456,-1.0])
        # xf2player = 1.*np.array([ 1.0,0.1234,-0.3, 1.0])
        # x02player = 2.*np.random.rand(4)-1.;
        # xf2player = 2.*np.random.rand(4)-1.;


        if True:  # cross shape
            x02player = 1.*np.array([-1.0,0.0,0.0,-1.0])
            xf2player = 1.*np.array([ 1.0,0.0,0.0, 1.0])
        if False: # x-shape
            x02player = 1.*np.array([-1.0,-1.0,0.-1.0,1.0])
            xf2player = 1.*np.array([ 1.0, 1.0,1.0,-1.0])
        if False: # squat x-shape
            x02player = 1.*np.array([-1.0,-0.5,0.-1.0,0.5])
            xf2player = 1.*np.array([ 1.0, 0.5,1.0,-0.5])



        # Fs  = [np.eye(nx)+dt*smat.block_diag(Asys,Asys) for _ in range(nt)]
        if False: # basic 
            G1shape = np.eye(nu1);
            G2shape = np.eye(nu2);
        if False: # best separation
            G1shape = np.diag([1.15,0.85]);
            G2shape = G1shape;        
        if True: # keeping ordering
            G1shape = np.diag([1.01,0.99]);
            G2shape = G1shape;        

        # G1shape = np.diag([1.03,0.97])@np.array([[1,1],[-1,1]]);

        Fs  = [smat.expm(dt*smat.block_diag(Asys,Asys)) for _ in range(nt)]
        G1s = [dt*np.vstack([G1shape,np.zeros([nx2,nu1])]) for _ in range(nt)]
        G2s = [dt*np.vstack([np.zeros([nx1,nu2]),G2shape]) for _ in range(nt)]


        Q1s = [dt*(Q1zero) for _ in range(nt)];
        Q2s = [dt*(Q2zero) for _ in range(nt)];
        P1s = [dt*(prox1*Qprox) for _ in range(nt)];
        P2s = [dt*(prox2*Qprox) for _ in range(nt)];


        if use_ffac: 
            for t in range(fwin): #np.min([10,nt])):
                Q1s[-t] = Q1s[-t] + dt*ffac*Q1zero
                Q2s[-t] = Q2s[-t] + dt*ffac*Q2zero

        R1s = [dt*R1ctrl for _ in range(nt)]
        R2s = [dt*R2ctrl for _ in range(nt)]
        R12s = [dt*R12ctrl for _ in range(nt)]
        R21s = [dt*R21ctrl for _ in range(nt)]




    alphas2 = np.linspace(0,1,nt)
    alphas1 = 1. - alphas2;
    xdes_v2 = np.outer(alphas1,x02player) + np.outer(alphas2,xf2player)


    if True: 
        # x02player_p1 = 1.*np.array([-1.0,0.0,0.0,-1.0])
        # xf2player_p1 = 1.*np.array([ 1.0,0.0,0.0, 1.0])
        # x02player_p2 = 1.*np.array([-1.0,0.0,0.0,-1.0])
        # xf2player_p2 = 1.*np.array([ 1.0,0.0,0.0, 1.0])        
        x1des_v2 = xdes_v2; 
        x2des_v2 = xdes_v2; 
    

    udes_v2 = np.zeros([nt,nu]);

    # udes = udes.reshape(nt*nu)

    LTI2player_v2 = {'dt':dt,
                     'x0':x02player,'nt':nt,
                    'nx':nx,'nu1':nu1,'nu2':nu1,
                    'As':Fs,
                    'xdes':xdes_v2, #traj_data['xdes'],
                    'x1des':x1des_v2,
                    'x2des':x2des_v2,
                    'udes':udes_v2,
                    'Q1s':Q1s,'Q2s':Q2s,
                    'P1s':P1s,'P2s':P2s,
                    'R1s':R1s,'R2s':R2s,
                    'R12s':R12s,'R21s':R21s,
                    'B1s':G1s,'B2s':G2s,
                    'paramsIndAgents':paramsIndAgents};

    start_time = time.time()
    LQ1 = LQ2GAME(ins=LTI2player_v2);
    LQ1.buildRollMats();

    end_time = time.time()
    print('time to compute rollouts: ',end_time - start_time);





LEGS = {};
LEGS['leg1'] = {};
dyy = 0.2;
sxx = 0.; axx = 0.2;
sxx2 = 1.0; axx2 = sxx2 + 0.3;
sz2 = 0.15;

tag = 'pts1'
LEGS['leg1'][tag] = {};
LEGS['leg1'][tag]['anns'] = {'u1':'Opt P1','u2':'Opt P2','une':'Nash','una':'Altu. Nash','uso':'Social Opt'}
LEGS['leg1'][tag]['locs'] = {'u1':[sxx,-0.*dyy],'u2':[sxx,-1.*dyy],'une':[sxx,-2.*dyy],'una':[sxx,-3.*dyy],'uso':[sxx,-4.*dyy]} 
LEGS['leg1'][tag]['sizes'] = {'u1':1,'u2':1,'une':1,'una':1,'uso':1,'usamp':1};
LEGS['leg1'][tag]['alocs'] = {'u1':[axx,-0.*dyy],'u2':[axx,-1.*dyy],'une':[axx,-2.*dyy],'una':[axx,-3.*dyy],'uso':[axx,-4.*dyy]}
LEGS['leg1'][tag]['asizes'] = {'u1':1,'u2':1,'une':1,'una':1,'uso':1}
LEGS['leg1'][tag]['shows'] = {'u1':True,'u2':True,'une':True,'una':True,'uso':True}

tag = 'ptss'
LEGS['leg1'][tag] = {};
LEGS['leg1'][tag]['anns'] = {'usamp':'Sample'};
LEGS['leg1'][tag]['locs'] = {'usamp':[sxx,-5.*dyy]};
LEGS['leg1'][tag]['sizes'] = {'usamp':1};
LEGS['leg1'][tag]['alocs'] = {'usamp':[axx,-5.*dyy]};
LEGS['leg1'][tag]['asizes'] = {'usamp':1};
LEGS['leg1'][tag]['shows'] = {'usamp':True};


tag = 'u1a'
LEGS['leg1'][tag] = {};
LEGS['leg1'][tag]['typs'] = {ftags[tag]:'hseg'}
LEGS['leg1'][tag]['anns'] = {ftags[tag]:'dJ1/du1'}
LEGS['leg1'][tag]['locs'] = {ftags[tag]:[sxx2,-0.*dyy]}
LEGS['leg1'][tag]['sizes'] = {ftags[tag]:sz2}
LEGS['leg1'][tag]['alocs'] = {ftags[tag]:[axx2,-0.*dyy]}
LEGS['leg1'][tag]['asizes'] = {ftags[tag]:1}
LEGS['leg1'][tag]['shows'] = {ftags[tag]:True}

tag = 'u1b'
LEGS['leg1'][tag] = {};
LEGS['leg1'][tag]['typs'] = {ftags[tag]:'hseg'}
LEGS['leg1'][tag]['anns'] = {ftags[tag]:'dJ1/du2'}
LEGS['leg1'][tag]['locs'] = {ftags[tag]:[sxx2,-1.*dyy]}
LEGS['leg1'][tag]['sizes'] = {ftags[tag]:sz2}
LEGS['leg1'][tag]['alocs'] = {ftags[tag]:[axx2,-1.*dyy]}
LEGS['leg1'][tag]['asizes'] = {ftags[tag]:1}
LEGS['leg1'][tag]['shows'] = {ftags[tag]:True}

tag = 'u2a'
LEGS['leg1'][tag] = {};
LEGS['leg1'][tag]['typs'] = {ftags[tag]:'hseg'}
LEGS['leg1'][tag]['anns'] = {ftags[tag]:'dJ2/du2'}
LEGS['leg1'][tag]['locs'] = {ftags[tag]:[sxx2,-2.*dyy]}
LEGS['leg1'][tag]['sizes'] = {ftags[tag]:sz2}
LEGS['leg1'][tag]['alocs'] = {ftags[tag]:[axx2,-2.*dyy]}
LEGS['leg1'][tag]['asizes'] = {ftags[tag]:1}
LEGS['leg1'][tag]['shows'] = {ftags[tag]:True}

tag = 'u2b'
LEGS['leg1'][tag] = {};
LEGS['leg1'][tag]['typs'] = {ftags[tag]:'hseg'}
LEGS['leg1'][tag]['anns'] = {ftags[tag]:'dJ2/du1'}
LEGS['leg1'][tag]['locs'] = {ftags[tag]:[sxx2,-3.*dyy]}
LEGS['leg1'][tag]['sizes'] = {ftags[tag]:sz2}
LEGS['leg1'][tag]['alocs'] = {ftags[tag]:[axx2,-3.*dyy]}
LEGS['leg1'][tag]['asizes'] = {ftags[tag]:1}
LEGS['leg1'][tag]['shows'] = {ftags[tag]:True}

tag = 'highexp1'; ptag = 'tframe'
LEGS['leg1'][tag] = {};
LEGS['leg1'][tag]['typs'] = {ftags[ptag]:'hseg'}
LEGS['leg1'][tag]['anns'] = {ftags[ptag]:'Gam1'}
LEGS['leg1'][tag]['locs'] = {ftags[ptag]:[sxx2,-4.*dyy]}
LEGS['leg1'][tag]['sizes'] = {ftags[ptag]:sz2}
LEGS['leg1'][tag]['alocs'] = {ftags[ptag]:[axx2,-4.*dyy]}
LEGS['leg1'][tag]['asizes'] = {ftags[ptag]:1}
LEGS['leg1'][tag]['shows'] = {ftags[ptag]:True}

tag = 'highexp2'; ptag = 'tframe'
LEGS['leg1'][tag] = {};
LEGS['leg1'][tag]['typs'] = {ftags[ptag]:'hseg'}
LEGS['leg1'][tag]['anns'] = {ftags[ptag]:'Gam2'}
LEGS['leg1'][tag]['locs'] = {ftags[ptag]:[sxx2,-5.*dyy]}
LEGS['leg1'][tag]['sizes'] = {ftags[ptag]:sz2}
LEGS['leg1'][tag]['alocs'] = {ftags[ptag]:[axx2,-5.*dyy]}
LEGS['leg1'][tag]['asizes'] = {ftags[ptag]:1}
LEGS['leg1'][tag]['shows'] = {ftags[ptag]:True}



# EXP 1 - SIDEBAR    

# start_time = time.time()
# SVO1 = QUADGAME2(ins={'version':'from_lq','LQ':LQ1})
# SVO1.basicEquilibria();
# end_time = time.time()
# print('time to generate quad game: ',end_time - start_time)
# SVO1.blowupsSVO();    

start_time = time.time()

data = {'version': quad_game_version}
if quad_game_version == 'from_data':
    data['M1'] = M1; data['M2'] = M2;
    data['d1'] = d1; data['d2'] = d2;
    data['center1'] = center1; data['center2'] = center2;

if quad_game_version == 'from_lq':
    data['LQ'] = LQ1;


SVO1 = QUADGAME2(ins=data)
SVO1.basicEquilibria();

if show_version == 'LQ' and False:
    #thh1 = np.pi/3; thh2 = np.pi/7;r1
    thh1 = (np.pi/2)*np.random.rand(1)[0];
    thh2 = (np.pi/2)*np.random.rand(1)[0];
    SVO1.computeFitting(ins={'th1':thh1,'th2':thh2})

nph = 25; nt = 200; # don't change for some figs
if True: nph = 40; # UPDATED FOR COEFF PLOTS
minph = 0.; maxph = np.pi/2; phbuff = 0.001; mint = 0.; maxt = 400; tbuff = 0.04; ph1 = 0; ph2 = np.pi/2;
nth1 = 30; nth2 = nth1;
paramsb = {'bunch1':False,'bunch2':False,'maxpreth1':4.5,'maxpreth2':4.5}
# phs = edgebunch_linspace(ph1+phbuff,ph2-phbuff,nth1=nth1,nth2=nth2,mid=0.5,params = paramsb);
phs = diagtan_linspace(np.array([1,0]),np.array([0,1]),nph = nph,buff=0.001);

bunchparams0 = {'nt':200,'bunch1':False,'bunch2':False,'maxpreth1':2.6,'maxpreth2':2.6}

scale_evecs = 4;
maxpreth = 5.;
bunchparams =  {'nt':200,'bunch1':True,'bunch2':True,'maxpreth1':maxpreth,'maxpreth2':maxpreth}

blowupdata = {};
blowupdata['tbnds'] = (mint,maxt)
blowupdata['nt4'] = nt;
blowupdata['ts'] = atan_linspace(mint+tbuff,maxt-tbuff,nt);
blowupdata['phs'] = phs; #np.linspace(minph+phbuff,maxph-phbuff,nph);
blowupdata['bunchparams'] = bunchparams;
blowupdata['bunchparams0'] = bunchparams0;
blowupdata['scale_evecs'] = scale_evecs;



SVO1.blowupsSVO(ins=blowupdata);
SVO1.integrateCtrls()
end_time = time.time()
print('time to generate quad game: ',end_time - start_time)

### INTEGRATING CURVES...
# start_time = time.time()
# end_time = time.time()


### bould 

bndpolys = {}

scale1 = 1.;
xmax = 1.65; xmin = -0.4; ymax = 2.; ymin = -0.35;
xmax2 = 1.2;
uframe2D = scale1*np.array([[xmax,ymax],[xmin,ymax],[xmin,ymin],[xmax,ymin]]);
uframe2Db = scale1*np.array([[xmax2,ymax],[xmin,ymax],[xmin,ymin],[xmax2,ymin]]);

scalec = 1.5
xmaxc = 1.5; xminc = -1.; ymaxc = 1.8; yminc = -0.7;
uframe2Dc = scalec*np.array([[xmaxc,ymaxc],[xminc,ymaxc],[xminc,yminc],[xmaxc,yminc]]);
scale1 = 1.; xmax = 1.8; xmin = -0.4; ymax = 2.; ymin = -0.4;
tframe = scale1*np.array([[xmax,ymax],[xmin,ymax],[xmin,ymin],[xmax,ymin]]);

scale1 = 0.6;
# minc = scale1*np.array([-2.,-1.2,-1.]);
# maxc = scale1*np.array([2.5,1.2,1.]);
minc = scale1*np.array([-1.,-1.,-1.]);
maxc = scale1*np.array([1.,1.,1.]);


uframe3Dnodes = np.array([[maxc[0],maxc[1],maxc[2]],[maxc[0],maxc[1],minc[2]],[maxc[0],minc[1],minc[2]],[maxc[0],minc[1],maxc[2]],
                          [minc[0],maxc[1],maxc[2]],[minc[0],maxc[1],minc[2]],[minc[0],minc[1],minc[2]],[minc[0],minc[1],maxc[2]]]);



scalel = 3.
xminl = -1.; xmaxl = 1.; yminl = -1.; ymaxl = 1.;
uframelq = scalel*np.array([[xmaxl,ymaxl],[xminl,ymaxl],[xminl,yminl],[xmaxl,yminl]]);

scale_nl = 15.
xmin_nl = -1.; xmax_nl = 1.; ymin_nl = -1.; ymax_nl = 1.;
uframe_nl1 = scale_nl*np.array([[xmax_nl,ymax_nl],[xmin_nl,ymax_nl],[xmin_nl,ymin_nl],[xmax_nl,ymin_nl]]);
uframe_nl2 = scale_nl*np.array([[xmax_nl,ymax_nl],[xmin_nl,ymax_nl],[xmin_nl,ymin_nl],[xmax_nl,ymin_nl]]);
uframe_nl3 = scale_nl*np.array([[xmax_nl,ymax_nl],[xmin_nl,ymax_nl],[xmin_nl,ymin_nl],[xmax_nl,ymin_nl]]);

# minlq = scale1*np.array([-1.,-1.]); # maxlq = scale1*np.array([1.,1.]);

bndcenter1 = np.zeros(3);
bndcenter2 = np.zeros(3);

if hasattr(SVO1,'center1'): bndcenter1 = SVO1.center1;
if hasattr(SVO1,'center2'): bndcenter2 = SVO1.center2;


# bndcenter1 = 0.5 * (SVO1.center1 + SVO1.center2) + np.array([0.,0.3,0.]);
# bndcenter2 = bndcenter1;

if len(bndcenter1)==3:
    uframe3Dnodes1 = uframe3Dnodes + bndcenter1;
    uframe3Dnodes2 = uframe3Dnodes + bndcenter2;

xminface = np.array([[minc[0],maxc[1],maxc[2]],[minc[0],maxc[1],minc[2]],[minc[0],minc[1],minc[2]],[minc[0],minc[1],maxc[2]]])
xmaxface = np.array([[maxc[0],maxc[1],maxc[2]],[maxc[0],maxc[1],minc[2]],[maxc[0],minc[1],minc[2]],[maxc[0],minc[1],maxc[2]]])
yminface = np.array([[maxc[0],minc[1],maxc[2]],[maxc[0],minc[1],minc[2]],[minc[0],minc[1],maxc[2]],[minc[0],minc[1],maxc[2]]])
ymaxface = np.array([[maxc[0],maxc[1],maxc[2]],[maxc[0],maxc[1],minc[2]],[minc[0],maxc[1],maxc[2]],[minc[0],maxc[1],maxc[2]]])
zminface = np.array([[maxc[0],maxc[1],minc[2]],[maxc[0],minc[1],minc[2]],[minc[0],minc[1],minc[2]],[minc[0],maxc[1],minc[2]]])
zmaxface = np.array([[maxc[0],maxc[1],maxc[2]],[maxc[0],minc[1],maxc[2]],[minc[0],minc[1],maxc[2]],[minc[0],maxc[1],maxc[2]]])
uframe3Dfaces = [xminface,xmaxface,yminface,ymaxface,zminface,zmaxface];

if len(bndcenter1)==3:
    uframe3Dfaces1 = [face + bndcenter1 for face in uframe3Dfaces]
    uframe3Dfaces2 = [face + bndcenter2 for face in uframe3Dfaces]


if len(bndcenter1)==3:    
    bndpolys['uframe3Dp1'] = POLYTOPE(ins={'nodes':uframe3Dnodes1,'faces':uframe3Dfaces1})
    bndpolys['uframe3Dp2'] = POLYTOPE(ins={'nodes':uframe3Dnodes2,'faces':uframe3Dfaces2})

bndpolys['uframe2D'] = POLYTOPE(ins={'nodes':uframe2D})
bndpolys['uframe2Db'] = POLYTOPE(ins={'nodes':uframe2Db})
bndpolys['uframe2Dc'] = POLYTOPE(ins={'nodes':uframe2Dc})
bndpolys['uframelq'] = POLYTOPE(ins={'nodes':uframelq})
bndpolys['uframe_nl1'] = POLYTOPE(ins={'nodes':uframe_nl1})
bndpolys['uframe_nl2'] = POLYTOPE(ins={'nodes':uframe_nl2})
bndpolys['uframe_nl3'] = POLYTOPE(ins={'nodes':uframe_nl3})

bndpolys['tframe'] = POLYTOPE(ins={'nodes':tframe})
bndpolys['uframe3D'] = POLYTOPE(ins={'nodes':uframe3Dnodes,'faces':uframe3Dfaces})



### OVERWRITTEN IN FRAMES???? REMOVED

gridparams1 = {};
gridparams2 = {};
ndiv1 = 10; ndiv2 = 10; ndivt = 4000;  ## MADE HUGE AND NOT SURE IF IT'S A PROBLEM
xdivs = np.linspace(0,np.pi/2,ndiv1+1); ydivs = np.linspace(0,np.pi/2,ndiv2+1);
gridparams1['xdivs'] = xdivs;
gridparams1['ydivs'] = ydivs;
tdivs = atan_linspace(0,200,ndivt+1);
gridparams1['tdivs'] = tdivs;
gridparams1['version'] = 'polargrid' # colorgrid
gridparams1['divtype'] = 'polar'; # square
gridparams1['direction'] = 'tan'
gridparams2 = {**gridparams1}
gridparams2['direction'] = 'cot'; 

# gridparams1['alpha_manual'] = 0.;
# gridparams2['alpha_manual'] = 0.;


gridparams3 = {**gridparams1}
gridparams4 = {**gridparams1}

# gridparams3['divtype'] = 'theta'
# gridparams4['divtype'] = 'theta'
gridparams3['direction'] = 'theta1'
gridparams4['direction'] = 'theta2'
# gridparams3['version'] = 'polargrid'
# gridparams4['version'] = 'polargrid'



fins = {}
fins['QG'] = SVO1

exptoshow = ['exp1','exp2']

showtags = {};
showtags['use_time_grad'] = True;
showtags['show_version'] = show_version
showtags['fig_version'] = fig_version


fins['gridparams'] = {};
fins['gridparams']['exp1'] = gridparams1;
fins['gridparams']['exp2'] = gridparams2;
fins['gridparams']['exp3'] = gridparams3;
fins['gridparams']['exp4'] = gridparams4;


ph1s = [np.arctan2(np.tan(th[1]),np.tan(th[0])) for th in th_highlights]
ph2s = [np.arctan2(1./np.tan(th[1]),np.tan(th[0])) for th in th_highlights]
t1s = [np.tan(th[0])/np.cos(ph1s[i]) for i,th in enumerate(th_highlights)]
t2s = [np.tan(th[0])/np.cos(ph2s[i]) for i,th in enumerate(th_highlights)]

fins['th_highlights'] = th_highlights;
fins['ph_highlights'] = {}
fins['t_highlights'] = {}
fins['ph_highlights']['exp1'] = ph1s;
fins['ph_highlights']['exp2'] = ph2s;
fins['t_highlights']['exp1'] = t1s
fins['t_highlights']['exp2'] = t2s

# fins['phwin'] = [np.pi/12,5*np.pi/12];
# fins['twin'] = [1,2];
# fins['exptoshow'] = ['exp1','exp2']
fins['bndpolys'] = bndpolys;
fins['frametags'] = ftags
fins['axisinfo'] = axisinfo;


if show_version == 'NL':
    fins['spatials'] = spatials;
    fins['spatial_agentids'] = spatial_agentids;

fins['extras'] = {'rad1':0.015,'rad2':0.015}
fins['extras']['radsamplet'] = 0.01
fins['extras']['radsampleu'] = 0.01

#########################
fins['segment_eps'] = 5.;
fins['LEGS'] = LEGS;
fins['current_ind'] = 60


## cmap options: 'viridis','plasma','inferno','magma','cividis','spring','summer','autumn','winter','cool','hot'
## cmap good opts: 'plasma','inferno 
cmap = plt.get_cmap('inferno'); nc = 20; adj = 0;
# fins['curveColorGrad'] = {'nchunks':15,'lw':lw,'ls':'solid','msty':msty,'msz':msz,'colorlist':[cmap((j)/(nc+adj)) for j in range(nc)]}
maxlight = 0.8;
lw = 4.0; msty = 's'; msz = 2; maxlght = 0.8; bcolor = [0,0,1]; pname='blarg';# (1) NOTE: maxpreth = 5.;
colorlist = [[bcolor[0]+(j/nc)*(1-bcolor[0])*maxlight,
              bcolor[1]+(j/nc)*(1-bcolor[1])*maxlight,
              bcolor[2]+(j/nc)*(1-bcolor[2])*maxlight,1.] for j in range(nc)]
fins['curveColorGrad'] = {'nchunks':15,'lw':lw,'ls':'solid','msty':msty,'msz':msz,'colorlist':colorlist}
fins['coeff_cutoff'] = 0.3;
coeff_negindstoshow =  [0,1,2,3,4,5,6,7,8,9,10,11]; #,1,2,3,4,5,6,7,8,9,10]; #[]  ## has priority over coeff_cutoff, shows negative 
fins['coeff_negindstoshow'] = coeff_negindstoshow;

# fins['basis_colormap'] = plt.get_cmap('viridis'); #parula_colormap
# fins['basis_colormap'] = parula_colormap; #
# fins['basis_max_val'] = 10;
# fins['basis_num_offset'] = 0;
fins['use_basis_colormap'] = True;
cmap = parula_colormap
# basis_cmap_vals = [0.1,0.3,0.4,0.5,0.55,0.7,0.8,0.9,0.99]fsh
basis_cmap_vals = [0.1,0.1,0.35,0.35,0.5,0.5,0.6,0.6,0.85,0.85,0.99,0.99];
fins['basis_colors'] = [cmap(basis_cmap_vals[j]) for j in coeff_negindstoshow];
fins['basis_scaling'] = 10.0;
fins['coeff_ybounds'] = [-1,1.5]; 
fins['coeffs_fontsize'] = 30;
fins['coeff_aspect'] = 2

# showtags['use_specific_sing_pts'] = False;
fins['specific_sing_pts'] = [0,1,2,3,4,5,6,7,8,9,10,11,12]


coeffind = 0;
coeffmaxvals = [30,30,  30,30,  30,30,  30,30,  20,10,  10,10];
fins['coeffgridparams'] = {'basecolor':list(cmap(basis_cmap_vals[coeffind])),'maxval':coeffmaxvals[coeffind],'colorind':coeffind}

fins['showtags'] = showtags;
fins['PENS'] = PENS
