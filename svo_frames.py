import numpy as np
import numpy.linalg as mat
import scipy as sp
import scipy.linalg as smat
import copy
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import warnings
warnings.filterwarnings('ignore')

from svo_classes import *



directions = [1,0,0]; #,0,0,0]; #,0]; #[2,0,0];
scale = 1.6; negscale= 0.4; #[2.2,1.6]; #0.8; #[0.2,0.2,0.2];
Ap = np.array([[1.25,0.35],[0.35,1.25]]).T; # A = np.hstack([Ap,Ap@U]) #CASE 3



def frame_SVO(ins={},
              show_opts = [],
              svo_opts = None,
              basis_opts1 = [],
              basis_opts2 = [],
              scale_opt = None,
              theta2_opt = []):

    # Step 2: Reload the module to apply update



    phstart = 2.8*(np.pi/12); phdel = 0.4*(np.pi/12); 
    phwin = [phstart,phstart+phdel]
    if svo_opts == 'svo1': tstart = 0.6; tdel = 0.1;
    if svo_opts == 'svo2': tstart = 2.6; tdel = 0.1;
    ins['sampwin'] = {'phwin':phwin,'twin':[tstart,tstart+tdel]};
        


    basis_opts = basis_opts1 + basis_opts2;
    basis_opts_even = []
    if (0 in basis_opts) or (1 in basis_opts): basis_opts_even.append(0)
    if (2 in basis_opts) or (3 in basis_opts): basis_opts_even.append(2)
    if (4 in basis_opts) or (5 in basis_opts): basis_opts_even.append(4)
    if (6 in basis_opts) or (7 in basis_opts): basis_opts_even.append(6)
    if (8 in basis_opts) or (9 in basis_opts): basis_opts_even.append(8)
    if (10 in basis_opts) or (11 in basis_opts): basis_opts_even.append(10)



    
    figsize = (12,10);
    if 'figsize' in ins: figsize = ins['figsize']
    
    fig = plt.figure(figsize=figsize)
    axs = [];
    axs_specs = {};
    dx1 = 0.55; dy1 = 0.2;
    axs_specs[0] = [0.05,0.35, dx1,dx1]
    axs_specs[1] = [0.05,0.05, dx1,dy1];
    axs_specs[2] = [0.65,0.5,  0.35,0.45];
    axs_specs[3] = [0.65,0.05,  0.35,0.35];
    ind = 0; axs.append(fig.add_axes(axs_specs[ind]))
    ind = 1; axs.append(fig.add_axes(axs_specs[ind]))
    ind = 2; axs.append(fig.add_axes(axs_specs[ind]))
    ind = 3; axs.append(fig.add_axes(axs_specs[ind]))

    ax = axs[0];
    ax1 = axs[1];
    ax2 = axs[2];
    ax3 = axs[3];

    ax1.set_aspect(dy1/dx1)



    loc = np.array([0,0])


    showtags = {}
    show_version = '2D'
    if 'showtags' in ins: showtags = ins['showtags'];
    if 'show_version' in showtags: show_version = showtags['show_version'];
    if 'fig_version' in showtags: fig_version = showtags['fig_version'];

    showtags['lq'] = True;
    showtags['ellipses'] = False;
    showtags['basic_struct'] = True;
    showtags['exp1'] = True;
    showtags['exp2'] = False;
    showtags['exp3'] = False;
    showtags['exp4'] = False;


    showtags['intersect'] = False;
    showtags['intersect1'] = False;
    showtags['intersect2'] = False;
    showtags['evecs'] = True;
    showtags['figs'] = [True,True,False,False];
    showtags['axes'] = [False,True,False,True];
    showtags['blowups'] = True;
    showtags['xblowups'] = False;
    showtags['high'] = False;
    showtags['hintersect'] = False;
    showtags['legend'] = False;
    showtags['trajs'] = True;
    showtags['upts'] = {'u1':True,'u2':True,'une':True,'una':True,'uso':True}
    showtags['tpts'] = {'u1':True,'u2':True,'une':True,'una':True,'uso':True}

    filename_add = '_exp2';
    filename = show_version + '_' + fig_version;
    filename = filename + filename_add

    # if showtags['exp1']: showtags['intersect1'] = True;
    # if showtags['exp2']: showtags['intersect2'] = True;


    showtags['coeff_plots'] = True;
    temp = True
    showtags['theta_space'] = True;
    showtags['theta_cutoff_square'] = False;
    showtags['theta_curves'] = False;#False; #stemp; 
    showtags['samples'] = False; #False; #temp;
    showtags['xdes'] = True;
    showtags['xne'] = True;
    showtags['xna'] = False;
    showtags['exact'] = True;
    showtags['approx'] = False;
    showtags['basis'] = True;
    showtags['basis_scaled'] = False;
    showtags['use_coeff_colors'] = True;












    isgif = False;
    if 'isgif' in ins: isgif = ins['isgif'];
    current_ind = 0;
    if 'current_ind' in ins: current_ind = ins['current_ind'];
    if isgif:
        current_ind = ins['current_ind'];
        all_inds = ins['all_inds'];

    current_inds = [current_ind]
    if 'current_inds' in ins: current_inds = ins['current_inds']


    #ins['current_ind'];

    segeps = 0.1;
    if 'segment_eps' in ins: segeps = ins['segment_eps']


    spatials = [[0,1],[2,3]];
    spatial_agentids = [0,1];
    if 'spatials' in ins: spatials = ins['spatials']
    if 'spatial_agentids' in ins: spatial_agentids = ins['spatial_agentids']



    ccolors = {'1a':[1,0,0.5],'1b':[0.5,0.,1.],'2a':[1.,0.,0.8],'2b':[0.8,0,1.]}
    # from svo_inputs import *
    ticksAlpha = 0.4;
    directions = [2,0,0];
    GAME1 = ins['QG'];
    # phs_to_highlight = data['phs_to_highlight'];
    # draw_ellipse_bnds = data['draw_ellipse_bnds'];
    LEGS = {};
    if 'LEGS' in ins: LEGS = ins['LEGS']


    # axs = plt.subplots(1,2); #add_subplot(121,aspect='equal');

    minx = 0.2; maxx = 1.2; miny = 0.2; maxy = 1.2;
    minx1 = 0.2; maxx1 = 1.2; minx2 = 0.2; maxx2 = 1.2;

    AXES3 = np.array([[1,0],[0,1],[-0.2,-0.7]]); #AXES3 = normalize(AXES3);
    AXES2 = np.array([[1,0],[0,1]]);
    AXES = AXES2;
    # draw_axes(ax,loc,AXES,normalize =False,head_width=0.05,scale=scale,**paramsAxes)
    # draw_axes(ax,loc,-AXES,normalize = False,head_width=0,scale=negscale,**paramsAxes)
    PENS = {};
    ph_highlights = {};
    th_highlights = {};
    t_highlights = {};
    # exptoshow = ['exp1','exp2'];
    ftags = {};
    if 'PENS' in ins: PENS = {**PENS,**ins['PENS']}

    frgbune = [0.5,0.,1.];
    PENS['exact'] = PEN({'frgba':[0,0,1,1],'ergba':frgbune + [0.7],'ls':'-','lw':4.,'zord':10})
    PENS['basis_scaled'] = PEN({'frgba':[0,0,0,1],'ergba':[0,0,1,0.6],'ls':'-','lw':4.,'zord':10})
    PENS['basis'] = PEN({'frgba':[0,0,0,1],'ergba':[0,0,0,1.],'ls':'-','lw':4.,'zord':11})



    if 'ph_highlights' in ins: ph_highlights = ins['ph_highlights']
    if 'th_highlights' in ins: th_highlights = ins['th_highlights']
    if 't_highlights' in ins: t_highlights = ins['t_highlights']
    # if 'exptoshow' in ins: exptoshow = ins['exptoshow'];
    if 'frametags' in ins: ftags = ins['frametags']

    wintag = 'exp1'
    if 'wintag' in ins: wintag = ins['wintag'];

    sampwin = {};
    sampwinb = {};
    if 'sampwin' in ins: sampwin = ins['sampwin']
    if 'sampwinb' in ins: sampwinb = ins['sampwinb']

    bndpolys = {};
    if 'bndpolys' in ins: bndpolys = ins['bndpolys']


    d1 = GAME1.d1; d2 = GAME1.d2;
    dd = GAME1.d1 + GAME1.d2;

    if show_version == 'LQ' or show_version == 'NL':
        nt = GAME1.LINQUAD.nt;
        nx = GAME1.LINQUAD.nx;
        nx1 = GAME1.LINQUAD.nx1;
        nx2 = GAME1.LINQUAD.nx2;
        nu = GAME1.LINQUAD.nu;
        nu1 = GAME1.LINQUAD.nu1;
        nu2 = GAME1.LINQUAD.nu2;

    if show_version == 'NL':
        xnl = GAME1.LINQUAD.xbase.reshape([nt,nx])
        unl1 = GAME1.LINQUAD.u1base;
        unl2 = GAME1.LINQUAD.u2base;#.reshape([nt,nu2])

        


    axisinfo  = {}
    axisinfo['axis1'] = {'nx':dd,'AXES':DEFAULT_AXES,'loc':np.array([0,0]),'scale':1.};
    axisinfo['axis2'] = {'nx':2,'AXES':DEFAULT_AXES,'loc':np.array([0,0]),'scale':1};
    axisinfo['axis3'] = {'nx':2,'AXES':DEFAULT_AXES,'loc':np.array([0,0.]),'scale':1.};
    axisinfo['axis4'] = {'nx':d1,'AXES':DEFAULT_AXES,'loc':np.array([0,0]),'scale':1.,'negscale':1.};
    
    if 'axisinfo' in ins: axisinfo = ins['axisinfo']
    axis1 = axisinfo['axis1']
    axis2 = axisinfo['axis2']
    axis3 = axisinfo['axis3']
    axis4 = axisinfo['axis4']


    params1 = {'pen_main':PENS['axes'],'ascale':1., 'headwid':0.03,**axis1}
    params2 = {'pen_main':PENS['axes'],'ascale':0.7,'headwid':0.03,**axis2}
    params3 = {'pen_main':PENS['axes'],'ascale':0.3,'headwid':0.03,**axis3}
    params4 = {'pen_main':PENS['axes'],'ascale':0.3,'headwid':0.03,**axis4}    

    # print(params1)
    # print(params2)

    AXIS1 = AXIS(ins=params1);
    AXIS2 = AXIS(ins=params2);
    AXIS3 = AXIS(ins=params3);
    AXIS4 = AXIS(ins=params4);
    if showtags['axes'][0]: AXIS1.draw(ax)
    if showtags['axes'][1] and showtags['theta_space']: AXIS2.draw(ax2)
    if showtags['axes'][2]: AXIS3.draw(ax)
    if showtags['axes'][3]: AXIS4.draw(ax)

    if showtags['legend']:
        lparams = {'name':'leg1',**axis4};
        LEG1 = LEGEND(ins=lparams);

    
    PEN3 = PEN({'frgba':[0,0,0,0.8],'ergba':[0,0,0,0.8]})

    # ######## WORKING
    extras = {};
    rad1 = 0.03; rad2 = 0.03;
    radsample = 0.01;
    radsamplet = 0.01; radsampleu = 0.01;
    if 'extras' in ins: extras = ins['extras']
    if 'rad1' in extras: rad1 = extras['rad1']
    if 'rad2' in extras: rad2 = extras['rad2']
    if 'radsample' in extras: radsample = extras['radsample']
    if 'radsamplet' in extras: radsamplet = extras['radsamplet']
    if 'radsampleu' in extras: radsampleu = extras['radsampleu']


    params1 = {'P':GAME1.M1,'r':GAME1.c1,'main':PENS['sph1'],'pairs':PENS['pairs1'],**axis1}
    params2 = {'P':GAME1.M2,'r':GAME1.c2,'main':PENS['sph2'],'pairs':PENS['pairs2'],**axis1}


    exptoshow = [];
    if 'exp1' in showtags:
        if showtags['exp1']: exptoshow.append('exp1');
    if 'exp2' in showtags:
        if showtags['exp2']: exptoshow.append('exp2');
    if 'exp3' in showtags:
        if showtags['exp3']: exptoshow.append('exp3');
    if 'exp4' in showtags:
        if showtags['exp4']: exptoshow.append('exp4');


    gridparams = {}; 
    if 'gridparams' in ins: gridparams = ins['gridparams']
    coeffgridparams = {};
    if 'coeffgridparams' in ins: coeffgridparams = ins['coeffgridparams'];


    if showtags['figs'][0] and (show_version == '2D' or show_version == '3D'):
        params1['rads'] = np.array([0.1,0.2,0.3,0.4]);
        params2['rads'] = 0.7*np.array([0.1,0.2,0.3,0.4]);
        QUAD1 = QUADFORM(ins=params1)
        QUAD2 = QUADFORM(ins=params2)
        QUAD1.drawLevelSets(ax)
        QUAD2.drawLevelSets(ax)

        inds1 = list(range(GAME1.d1))
        inds2 = [i + GAME1.d1 for i in range(GAME1.d2)]

        params1a = {'A':GAME1.M1[:,inds1].T,'b':-GAME1.c1[inds1],'main':PENS['opt1a'],'bnds':bndpolys,**axis1}
        params1b = {'A':GAME1.M1[:,inds2].T,'b':-GAME1.c1[inds2],'main':PENS['opt1b'],'bnds':bndpolys,**axis1}
        params2a = {'A':GAME1.M2[:,inds2].T,'b':-GAME1.c2[inds2],'main':PENS['opt2a'],'bnds':bndpolys,**axis1} 
        params2b = {'A':GAME1.M2[:,inds1].T,'b':-GAME1.c2[inds1],'main':PENS['opt2b'],'bnds':bndpolys,**axis1}

        OPT1a = AFFINE(ins=params1a)
        OPT1b = AFFINE(ins=params1b)
        OPT2a = AFFINE(ins=params2a)
        OPT2b = AFFINE(ins=params2b)

        OPT1a.computeHullBnds(ftags['u1a']); OPT1a.drawInHull(ax,ftags['u1a']);
        OPT1b.computeHullBnds(ftags['u1b']); OPT1b.drawInHull(ax,ftags['u1b']);
        OPT2a.computeHullBnds(ftags['u2a']); OPT2a.drawInHull(ax,ftags['u2a']);
        OPT2b.computeHullBnds(ftags['u2b'],eps=-0.07); OPT2b.drawInHull(ax,ftags['u2b']);


    # print(ftags['u1a'])
    # print(ftags['u1b'])addObj
    # print(ftags['u2a'])
    # print(ftags['u2b'])

        if showtags['legend']:
            ltag = 'leg1';
            tag = 'u1a'; LEG1.addObj(ins={'tag':tag,'obj':OPT1a,'tags':[ftags[tag]],'show':True,'info':LEGS[ltag][tag]})
            tag = 'u1b'; LEG1.addObj(ins={'tag':tag,'obj':OPT1b,'tags':[ftags[tag]],'show':True,'info':LEGS[ltag][tag]})
            tag = 'u2a'; LEG1.addObj(ins={'tag':tag,'obj':OPT2a,'tags':[ftags[tag]],'show':True,'info':LEGS[ltag][tag]})
            tag = 'u2b'; LEG1.addObj(ins={'tag':tag,'obj':OPT2b,'tags':[ftags[tag]],'show':True,'info':LEGS[ltag][tag]})

        # LEG1.draw(ax);
        # OPT2b.computeHullBnds(ftags['u2b']); OPT2b.draw(ax);

        #params1a = {'A':GAME1.M1[:,[0]].T,'b':-GAME1.c1[[0]],'main':PENS['opt1a'],'bnds':bndpolys}   
        #params = {};
        # params = {};        
        # ndiv1 = 10; ndiv2 = 10; ndivt = 10;
        # xdivs = np.linspace(0,np.pi/2,ndiv1+1); ydivs = np.linspace(0,np.pi/2,ndiv2+1);
        # params['xdivs'] = xdivs; params['ydivs'] = ydivs;
        # tdivs = atan_linspace(0,100,ndivt+1);    
        # params['tdivs'] = tdivs;
        # params['version'] = 'polargrid' # colorgrid
        # params['divtype'] = 'polar'; # square    



        # uCURVES1 = {};
        # tCURVES1 = {};

        hind = current_ind;
        if isgif:
            if len(ph_highlights['exp1'])==len(all_inds): hinds = current_ind

        ph1 = ph_highlights['exp1'][hind];
        ph2 = ph_highlights['exp2'][hind];
        t1 = t_highlights['exp1'][hind];
        t2 = t_highlights['exp2'][hind];
        EXP1 = GAME1.EXP['exp1']['computed_curves_full']
        EXP2 = GAME1.EXP['exp2']['computed_curves_full']
        tint1 = list(EXP1['tcurves'])[0][1]
        tint2 = list(EXP2['tcurves'])[0][1]
        ph1snaps = np.array([temp[0] for temp in list(EXP1['tcurves'])]);
        ph2snaps = np.array([temp[0] for temp in list(EXP2['tcurves'])]);
        ind1 = np.argmin(np.abs(ph1snaps-ph1));
        ind2 = np.argmin(np.abs(ph2snaps-ph2));
        ph1 = ph1snaps[ind1];
        ph2 = ph2snaps[ind2];
        tsnaps1 = EXP1['tcurves'][(ph1,tint1)]
        tsnaps2 = EXP2['tcurves'][(ph2,tint2)]

        ind1 = np.argmin(np.abs(tsnaps1-t1));
        ind2 = np.argmin(np.abs(tsnaps2-t2));
        t1high = tsnaps1[ind1];
        t2high = tsnaps2[ind2];
        usamp = EXP1['ucurves'][(ph1,tint1)][ind1];
        # u2 = EXP2['ucurves'][(ph2,tint2)][ind2];
        thsamp = EXP1['thcurves'][(ph1,tint1)][ind1];
        # th2 = EXP1['thcurves'][(ph2,tint2)][ind2];
        # print('u1: ',u1)
        # print('u2: ',u2)
        # print('th expected: ',(np.pi/3,np.pi/3))
        # print('th1: ',th1)
        # print('th2: ',th2)
        # asdasdfa
        ph1high = ph1;
        ph2high = ph2;


    for tag in exptoshow:
        ### MAIN CURVES 

        
        EXP = GAME1.EXP[tag]['computed_curves_full']
        ucurves = EXP['ucurves'];
        thcurves = EXP['thcurves'];
        tcurves = EXP['tcurves'];


        uparams1 = {'curves':ucurves,'tcurves':tcurves,'bnds':bndpolys}
        thparams1 = {'curves':thcurves,'tcurves':tcurves,'bnds':bndpolys}

        uparams1 = {**uparams1,**axis1}
        thparams1 = {**thparams1,**axis2}

        uparams1['curveColorInputs'] = gridparams[tag];
        uparams1['curveColorFunc'] = curveColor2DGrid1;
        uparams1['pen_main'] = PENS['tphcurvs']


        thparams1['curveColorInputs'] = gridparams[tag];
        thparams1['curveColorFunc'] = curveColor2DGrid1 

        if showtags['use_coeff_colors']:
            # coeffgridparams[tag]['data'] = GAME1.EXP[tag];
            gridparams[tag] = {**gridparams[tag],**coeffgridparams}
            gridparams[tag]['data'] = GAME1.EXP[tag]['coeffs_color_data']
            thparams1['curveColorInputs'] = gridparams[tag];
            thparams1['curveColorFunc'] = curveColorByCoeffs

            ### BLARG ### BLARG ### BLARG ### BLARG ### BLARG ### BLARG ### BLARG ### BLARG 
            ### BLARG ### BLARG ### BLARG ### BLARG ### BLARG ### BLARG ### BLARG ### BLARG 
            ### BLARG ### BLARG ### BLARG ### BLARG ### BLARG ### BLARG ### BLARG ### BLARG


        thparams1['pen_main'] = PENS['tphcurvs']


        if showtags['figs'][0] and (show_version == '2D' or show_version == '3D'):
            uCURVES1 = CURVES(ins=uparams1);
            uCURVES1.filterAndSegment(ftags['uframe'],eps=segeps);
            uCURVES1.drawCurves(ax,ftags['uframe'],typ='segmented');


        thCURVES1 = CURVES(ins=thparams1);
        thCURVES1.filterAndSegment(ftags['tframe'],eps=segeps);
        #     # tCURVES.filterAndSegment('tframe')
        if showtags['theta_curves']:            
            thCURVES1.drawCurves(ax2,ftags['tframe'],typ='segmented');   #### DRAWING MAIN THETA CURVES

        ### HIGHLIGHT CURVES
        TRAJS = GAME1.EXP[tag]['trajs_full']
        phsnaps = np.array(list(TRAJS));
        # if False: 

        uCURVESh = {}
        thCURVESh = {}


        # print(ph_highlights['exp1'])
        # print(ph_highlights['exp2'])
        # asdfasdfas

        if tag in ph_highlights:
            phs = ph_highlights[tag];

            for curr_ind in current_inds:
                phs0 = [phs[curr_ind]]
                for k,ph0 in enumerate(phs0):
                    ind = np.argmin(np.abs(phsnaps-ph0));
                    ph = phsnaps[ind];
                    intervals = list(TRAJS[ph])
                    ucurvesh = {(ph,interv):TRAJS[ph][interv]['ctrls'] for interv in intervals}
                    thcurvesh = {(ph,interv):TRAJS[ph][interv]['ths'] for interv in intervals}
                    tcurvesh = {(ph,interv):TRAJS[ph][interv]['ts'] for interv in intervals}

                    uparamsh = {'curves':ucurvesh,'tcurves':tcurvesh,'bnds':bndpolys}
                    thparamsh = {'curves':thcurvesh,'tcurves':tcurves,'bnds':bndpolys}
                    uparamsh = {**uparamsh,**axis1}
                    thparamsh = {**thparamsh,**axis2}

                    # uparamsh['curveColorInputs'] = gridparams[tag];
                    # uparamsh['curveColorFunc'] = curveColor2DGrid1 
                    uparamsh['pen_main'] = PENS['high'+tag]
                    # thparamsh['curveColorInputs'] = gridparams[tag];
                    # thparamsh['curveColorFunc'] = curveColor2DGrid1 
                    thparamsh['pen_main'] = PENS['high'+tag]
                    uparamsh['use_shadow'] = True;
                    thparamsh['use_shadow'] = True;


                    if showtags['figs'][0] and (show_version == '2D' or show_version == '3D'):
                        uCURVESh[(tag,ph0)] = CURVES(ins=uparamsh);
                        uCURVESh[(tag,ph0)].filterAndSegment(ftags['uframe'],eps=segeps);


                    thCURVESh[(tag,ph0)] = CURVES(ins=thparamsh);
                    thCURVESh[(tag,ph0)].filterAndSegment(ftags['tframe'],eps=segeps);
                    # tCURVES.filterAndSegment('tframe')
                    if showtags['high']:
                        if showtags['figs'][0] and (show_version == '2D' or show_version == '3D'):
                            uCURVESh[(tag,ph0)].drawCurves(ax,ftags['uframe'],typ='segmented');
                        if showtags['theta_space']:                            
                            thCURVESh[(tag,ph0)].drawCurves(ax2,ftags['tframe'],typ='segmented');

                        if k==0:
                            if showtags['legend']:
                                otag = 'high' + tag;
                                LEG1.addObj(ins={'tag':otag,'obj':thCURVESh[(tag,ph0)],'tags':[ftags['tframe']],'show':True,'info':LEGS[ltag][otag]})

        if showtags['blowups']:
            bparams2 = {'bnds':bndpolys,'rad':rad2,**axis2,'pen_main':PENS['tblow'],'id':'blowups','dtyp':'marker'}
            cparams2 = {'bnds':bndpolys,'rad':rad2,**axis2,'pen_main':PENS['cutoff'],'id':'cutoffs','dtyp':'marker'}
            # wparams2 = {'bnds':bndpolys,'rad':rad2,**axis2,'pen_main':PENS['wblow'],'drawtype':'marker'}
            # tparamsb['curveColorInputs'] = gridparams[wintag];
            # tparamsb['curveColorFunc'] = curveColor2DGrid1            
            bPTS = PTS(ins=bparams2);
            cPTS = PTS(ins=cparams2);
            # wPTS = PTS(ins=wparams2);
            PROBPTS = GAME1.EXP[tag]['probpts'];
            PROBVIS = GAME1.EXP[tag]['probVis'];
            CUTPTS = GAME1.EXP[tag]['pdcutoffs']


            if False: 
                for ph in PROBPTS:
                    for kk,ths in enumerate(PROBPTS[ph]):
                        probtag = (ph,kk)
                        # if tag == 'exp1': th1,th2 = np.arctan2(tt*np.cos(ph),1),np.arctan2(tt*np.sin(ph),1)
                        # if tag == 'exp2': th1,th2 = np.arctan2(1,tt*np.cos(ph)),np.arctan2(tt*np.sin(ph),1)
                        if ths[0] <= GAME1.maxtheta1 and ths[1] <= GAME1.maxtheta2:
                            if (len(ins['specific_sing_pts'])==0) or (kk in ins['specific_sing_pts']):
                                bPTS.addPts(pts = {probtag:PROBPTS[ph][kk]},ins={'pens':{probtag:PENS['tblow1']}})
                            else: bPTS.addPts(pts = {probtag:PROBPTS[ph][kk]},ins={'pens':{probtag:PENS['tblow3']}})
                        else:
                            if (len(ins['specific_sing_pts'])==0) or (kk in ins['specific_sing_pts']):
                                bPTS.addPts(pts = {probtag:PROBPTS[ph][kk]},ins={'pens':{probtag:PENS['tblow2']}})
                            else: bPTS.addPts(pts = {probtag:PROBPTS[ph][kk]},ins={'pens':{probtag:PENS['tblow3']}})


            for ph in CUTPTS:
                cPTS.addPts(pts = {ph:CUTPTS[ph]})

                    # wPTS.addPts(pts = {tag:thstemp[j]})
            bPTS.filterPts(ftags['tframe']); #,params={'filterstyle':'boundary','nbndpts':20});
            # cPTS.filterPts(ftags['tframe']); #,params={'filterstyle':'boundary','nbndpts':20});
            # wPTS.filterPts(ftags['tframe']);
            if showtags['theta_space']:
                bPTS.drawPts(ax2,ftags['tframe'])
                cPTS.drawPts(ax2,ftags['tframe'],typ='pts')
            # wPTS.drawPts(ax)



    if showtags['theta_space'] and showtags['theta_cutoff_square']:
        mth1 = GAME1.maxtheta1 
        mth2 = GAME1.maxtheta2 
        if mth1 >= np.pi/2: mth1 = np.pi/2
        if mth2 >= np.pi/2: mth2 = np.pi/2
        temp1 = axis2['scale']*np.array([[0.,mth2],[mth1,mth2]]) + axis2['loc'];
        temp2 = axis2['scale']*np.array([[mth1,0.],[mth1,mth2]]) + axis2['loc'];
        ax2.plot(temp1[:,0],temp1[:,1],color='black',lw=2)
        ax2.plot(temp2[:,0],temp2[:,1],color='black',lw=2)
        temp1 = axis2['scale']*np.array([[0.,0.],[mth1,0.]]) + axis2['loc'];
        temp2 = axis2['scale']*np.array([[0.,0.],[0.,mth2]]) + axis2['loc'];
        ax2.plot(temp1[:,0],temp1[:,1],color='black',lw=2)
        ax2.plot(temp2[:,0],temp2[:,1],color='black',lw=2)


    
    # temp1 = np.array([[0.,[0.,GAME1.lammin1])
    # temp2 = 


    phs0 = list(GAME1.EXP[wintag]['trajs'])
    phtsinwin = {};
    phwin = [-1,-1]; twin = [-1,-1];
    if 'twin' in sampwin:
        sampleboxtyp = 'box'
        twin = sampwin['twin'];
        phwin = sampwin['phwin'];
    elif 'thcenter' in sampwin:
        sampleboxtyp = 'circle'
        thcenter = sampwin['thcenter'];
        thrad = sampwin['thrad'];
        if wintag == 'exp1': phwinver = 'tan'
        if wintag == 'exp2': phwinver = 'cot'
        phwin = phwin_fromcircle(thcenter,thrad,version=phwinver);

    for ph0 in phs0:
        if (phwin[0] <= ph0) & (ph0 <= phwin[1]):
            phtsinwin[ph0] = twin; #phsinwin.append(ph0)


    phs0 = list(GAME1.EXP[wintag]['trajs'])
    phtsinwinb = {};
    phwinb = [-1,-1]; twinb = [-1,-1];
    if 'twin' in sampwinb:
        sampleboxtypb = 'box'
        twinb = sampwinb['twin'];
        phwinb = sampwinb['phwin'];
    elif 'thcenter' in sampwinb:
        sampleboxtypb = 'circle'
        thcenterb = sampwinb['thcenter'];
        thradb = sampwinb['thrad'];
        if wintag == 'exp1': phwinverb = 'tan'
        if wintag == 'exp2': phwinverb = 'cot'
        phwinb = phwin_fromcircle(thcenterb,thradb,version=phwinverb);
    for ph0 in phs0:
        if (phwinb[0] <= ph0) & (ph0 <= phwinb[1]):
            phtsinwinb[ph0] = twinb; #phsinwin.append(ph0)


    #phwin2 = [2.2*np.pi/12,3.8*np.pi/12];
    twin2 = [0,10];
    PROBTS = GAME1.EXP[wintag]['probts'];
    PROBPTS = GAME1.EXP[wintag]['probpts'];
    PROBVIS = GAME1.EXP[wintag]['probVis'];
    PROBXVIS = GAME1.EXP[wintag]['probxVis'];
    # print(list(PROBTS))
    # print(list(PROBPTS))
    # print(list(PROBVIS))
    # print(list(PROBXVIS))

    if showtags['xblowups'] and (show_version == 'LQ' or show_version == 'NL'):
        bparams2 = {'bnds':bndpolys,'rad':rad2,**axis2,'pen_main':PENS['tblow'],'id':'blowups','dtyp':'marker'}
        bPTSh = PTS(ins=bparams2);
        xblowups = {}; #thblowups = {};
        x1btags = [];
        x2btags = [];
        xpens = {}; thpens = {};
        # for ph in PROBPTS:
        for ph in phtsinwinb:
            # if (phwin[0] < ph) & (ph < phwin[1]):
            #twin2 = phtsinwin[ph]
            ##########################################
            for kk,tt in enumerate(PROBTS[ph]):
                tinside = False;
                if sampleboxtypb == 'box':
                    if twinb[0] <= tt and tt <= twinb[1]:
                        tinside = True;
                if sampleboxtypb == 'circle':
                    th1 = np.arctan2(tt*np.cos(ph),1)
                    if wintag == 'exp1': th2 = np.arctan2(tt*np.sin(ph),1)
                    if wintag == 'exp2': th2 = np.arctan2(1,tt*np.sin(ph))
                    if mat.norm(np.array([th1,th2]) - thcenterb) < thradb: tinside = True;

                if tinside: #twintemp[0] <= tt and tt<= twintemp[1]:
                #     probtag = (phh,tt)
                # if (twin2[0] < tt) & (tt < twin2[1]):
                    probtag = (ph,kk)
                    ths = PROBPTS[ph][kk];
                    # thblowups[probtag] = PROBPTS[ph][kk];
                    xtemps = PROBXVIS[ph][tt]

                    for ll,xtemp in enumerate(xtemps):
                        xtemp = xtemp.reshape([nt,nx]);
                        xtemp1 = xtemp[:,:nx1]
                        xtemp2 = xtemp[:,nx1:]

                        if True: 
                            for jj,spat in enumerate(spatials):
                                xtempj = xtemp[:,spat] + xnl[:,spat];
                                if ths[0] <= GAME1.maxtheta1 and ths[1] <= GAME1.maxtheta2:
                                    aid = spatial_agentids[jj]
                                    bPTSh.addPts(pts = {probtag:PROBPTS[ph][kk]},ins={'pens':{probtag:PENS['tblow1']}})
                                    xblowups[(ph,kk,aid,ll,jj)] = xtempj; #np.array([xtemp1]);
                                    # xblowups[(ph,kk,1,ll,jj)] = xtemp2; #np.array([xtemp2]);
                                    xpens[(ph,kk,aid,ll,jj)] = PENS['xblowups'];
                                    # xpens[(ph,kk,1,ll)] = PENS['xblowups'];
                                    if aid == 0: x1btags.append((ph,kk,aid,ll,jj))
                                    if aid == 1: x2btags.append((ph,kk,aid,ll,jj))


        if showtags['theta_space']:
            bPTSh.drawPts(ax2,ftags['tframe'])
        if wintag == 'exp1': xbase = GAME1.xne
        if wintag == 'exp2': xbase = GAME1.xopt1;
        xbase = xbase.reshape([nt,nx]);

        xblowups['xbase1'] = xbase[:,:nx1]; #np.array([xbase[:,:nx1]])
        xblowups['xbase2'] = xbase[:,nx1:]; #np.array([xbase[:,nx1:]]);
        if show_version == 'NL':
            xblowups['xbase1'] = xbase[:,spatials[0]] + xnl[:,spatials[0]];
            xblowups['xbase2'] = xbase[:,spatials[1]] + xnl[:,spatials[1]];

        xpens['xbase1'] = PENS['xbase'];
        xpens['xbase2'] = PENS['xbase'];


        #TRAJSb = TRAJECTORIES(pts=xblowups,ins={'pen_main':PENS['ablowups'],'pens':xpens,**axis4})
        TRAJSb = CURVES(ins={'curves':xblowups,'bnds':bndpolys,
                             'pen_main':PENS['ablowups'],'pens':xpens,**axis4})
        TRAJSb.filterAndSegment(ftags['uframe'],eps=segeps);
        if showtags['trajs']:
            TRAJSb.drawCurves(ax,curves=['xbase1','xbase2'],poly=ftags['uframe'],typ='segmented');
        # TRAJSb.drawCurves(ax,poly=ftags['uframe'],typ='segmented');

        temptags = x1btags + x2btags + ['xbase1','xbase2']
        awid = 0.001
        if showtags['evecs']: 
            # TRAJSb.plotTrajs(ax,tags=temptags);
            TRAJSb.drawCurves(ax,curves=temptags,typ='pts');

        temptags = x1btags + x2btags;
        awid = 0.001
        if showtags['evecs']: 
            TRAJSb.drawDiffs(ax,'xbase1',x1btags,version='line');#arrow1to2',ins={'wid':awid,'hwid':1.2*awid})
            TRAJSb.drawDiffs(ax,'xbase2',x2btags,version='line');#arrow1to2',ins={'wid':awid,'hwid':1.2*awid})


    ############################################################################################################
    ############################################################################################################
    ############################################################################################################
    ############################################################################################################

    # if showtags['basic_struct']:
    # if showtags['basic_struct']:
    # if showtags['basic_struct']:
    # if showtags['basic_struct']:
    # if showtags['basic_struct']:
    # if showtags['basic_struct']:
    # if showtags['basic_struct']:

    #### ONLY FOR SVPO

    if showtags['basic_struct']:

        tag = 'exp1';
        EXP = GAME1.EXP[tag];
        xxi = GAME1.xxi_exp1; # self.una - self.une;
        eigs = EXP['eigs'];
        VVs = EXP['Vs'];
        VVis = EXP['Vis']
        probts = EXP['probts']


        ### DEFAULTS
        phs = list(eigs);
        ph = phs[int(len(phs)/2)]; # for ph in [phs[10]]:
        tt = 0.2;


        if 'phwin' in sampwin:
            ph0 = 0.5*(sampwin['phwin'][-1]+sampwin['phwin'][0]);
            ind = np.argmin(np.abs(ph0 - np.array(phs)))
            ph = phs[ind];

        if 'twin' in sampwin: tt = 0.5*(sampwin['twin'][-1]+sampwin['twin'][0]);

        eig = eigs[ph]
        probt = probts[ph]
        eig_real = np.real(eig);
        VV = VVs[ph]
        VVi = VVis[ph]

        sortinds = np.argsort(eig_real)
        eig_sorted = eig[sortinds];
        VV_sort = VV[:, sortinds]
        VVi_sort = VVi[sortinds,:];

        U2 = (1./np.sqrt(2))*np.array([[1.,-1j],[1.,1j]]); temp = []
        i = 0; 
        while i < len(eig_sorted):
            # if not(np.isclose(np.imag(eig[i]),0)): temp.append(U2); i = i + 2;
            if not(np.imag(eig_sorted[i])==0): temp.append(U2); i = i + 2;
            else: temp.append(1); i = i + 1;

        TRANS = smat.block_diag(*temp);
        V = VV_sort@TRANS
        LAM = np.real(TRANS.conj().T@np.diag(eig_sorted)@TRANS)
        Vi = TRANS.conj().T@VVi_sort;
        deig = np.diag(LAM)
        
        coeff_cutoff = 0.1;
        if 'coeff_cutoff' in ins: coeff_cutoff = ins['coeff_cutoff'];
        negindstoshow = [];
        if 'coeff_negindstoshow' in ins: negindstoshow = ins['coeff_negindstoshow'];
        basis_scaling = 1.;
        if 'basis_scaling' in ins: basis_scaling = ins['basis_scaling'];
        
        geig = 1./(1+deig/tt);
        coeffs2 = xxi@V@np.diag(geig)

        if len(negindstoshow) > 0: inds = np.where(deig < 0.)[0][::-1]; inds = inds[negindstoshow];
        else: mask = np.abs(coeffs2) > coeff_cutoff; inds = np.where(mask)[0];
        # inds = np.where(mask)[0];




        deig_trimmed = deig[inds];
        coeffs2_trimmed = coeffs2[inds]
        toshow_Vi = basis_scaling*Vi[inds]
        toshow_coeffs2 = coeffs2[inds]*(1/basis_scaling)


        sgns = np.array([1,-1,1,1,1,1,1,1,1,1,1,1])
        # sgns = np.sign(toshow_coeffs2)
        toshow_coeffs2 = sgns*toshow_coeffs2
        toshow_Vi = np.diag(sgns)@toshow_Vi


        if showtags['coeff_plots']:
            #fig,axb = plt.subplots(1,1,figsize=(7,7))
            axb = axs[1]
            # axs = [axs]
            # axs[0].plot(coeffs,color='grey')
            # axs[0].plot(geig,color='blue')
            positions = np.array(list(range(len(inds))))+1
            for jj,_ in enumerate(positions): #pos,_,coeff in enumerate(zip(positions,np.ones(len(inds)),toshow_coeffs2)): #coeffs2_trimmed):
                pos = positions[jj]; coeff = toshow_coeffs2[jj]; clr = ins['basis_colors'][jj]
            # for pos,coeff,clr in zip(positions,toshow_coeffs2,ins['basis_colors']): #coeffs2_trimmed):
                # pos = positions[jj]; coeff = toshow_coeffs2[jj]                
                axb.bar(pos, np.abs(coeff), width=0.8, align='center', edgecolor='black', color=clr,alpha=0.8)
            # axs[0].plot(coeffs2_trimmed,color='red',marker='o')f
            ymax = np.max(np.abs(toshow_coeffs2));

            if not(ins['coeff_aspect']==None): axb.set_aspect(ins['coeff_aspect']);
            elif len(ins['coeff_ybounds'])>0: axb.set_ylim(ins['coeff_ybounds']);
            else: axb.set_ylim([-1.25*ymax,1.25*ymax])
            coeffs_fontsize = 20;
            if 'coeffs_fontsize' in ins: coeffs_fontsize = ins['coeffs_fontsize'];
            # axb.set_xticklabels(r'$j$',fontsize=coeffs_fontsize)
            # axb.set_yticklabels(fontsize=coeffs_fontsize)
            axb.tick_params(axis='both', which='major', labelsize=coeffs_fontsize)  # Major ticks
            axb.tick_params(axis='both', which='minor', labelsize=coeffs_fontsize)  # Minor ticks (if any)            
            axb.set_xlabel(r'$\ell$',fontsize=coeffs_fontsize)
            axb.set_xticks(positions)
            axb.set_xticklabels(positions)
            
            axb.set_ylim()
            if False:
                # Example data
                x = [1, 2.5, 4, 6.5]  # Center positions of the bars
                widths = [1, 1.5, 2, 1]  # Widths of the bars
                heights = [3, 7, 5, 9]  # Heights of the bars


                # Loop through the data to plot each bar with variable width
                for xi, width, height in zip(x, widths, heights):
                    ax.bar(xi, height, width=width, align='center', edgecolor='black', color='skyblue')

                # # Add labels and title
                # ax.set_xlabel('X-axis', fontsize=12)
                # ax.set_ylabel('Y-axis', fontsize=12)
                # ax.set_title('Bar Plot with Variable Width Bars', fontsize=14)

                # Customize tick marks
                ax.set_xticks(x)
                ax.set_xticklabels([f'Bar {i+1}' for i in range(len(x))], fontsize=10)

        # # Adjust layout for better spacing
        # plt.tight_layout()

        # # Show the plot
        # plt.show()        

        ################################################
        ################################################
        # uparams2 = {'bnds':bndpolys,'rad':radsampleu,**axis1,'pen_main':PENS['usamp']}
        # tparams2 = {'bnds':bndpolys,'rad':radsamplet,**axis2,'pen_main':PENS['tsamp']}
        # uparams2['curveColorInputs'] = gridparams[wintag];
        # tparams2['curveColorInputs'] = gridparams[wintag];
        # uparams2['curveColorFunc'] = curveColor2DGrid1
        # tparams2['curveColorFunc'] = curveColor2DGrid1
        # uPTS = PTS(ins=uparams2);
        # tPTS = PTS(ins=tparams2);

        AA  = GAME1.LINQUAD.AA;
        BB1 = GAME1.LINQUAD.BB1;
        BB2 = GAME1.LINQUAD.BB2;
        BB = np.hstack([BB1,BB2]);
        HH = GAME1.LINQUAD.HH;
        x0 = GAME1.LINQUAD.x0;
        xdes = GAME1.LINQUAD.xdes.reshape([nt,nx])
        une = GAME1.une;
        una = GAME1.una;

        xtrajsb = {}; xpensb = {};
        xtrajs1 = {}; xpens1 = {};
        xtrajs2 = {}; xpens2 = {};

        if showtags['xdes']:
            xtag = 'xdes'
            ind = 0; spat = spatials[ind]; xtrajsb[(xtag,ind)] = xdes[:,spat];
            ind = 1; spat = spatials[ind]; xtrajsb[(xtag,ind)] = xdes[:,spat];
            ind = 0; xpensb[(xtag,ind)] = PENS['xdes'];#PEN({'ergba':[0,0,0,1],'lw':templw,'ls':templs,'msty':tempmsty,'msz':tempmsz})
            ind = 1; xpensb[(xtag,ind)] = PENS['xdes'];#PEN({'ergba':[0,0,0,1],'lw':templw,'ls':templs,'msty':tempmsty,'msz':tempmsz})            
            if not(PENS['xne1']==None): ind = 0; xpensb[(xtag,ind)] = PENS['xdes1'];
            if not(PENS['xne2']==None): ind = 1; xpensb[(xtag,ind)] = PENS['xdes2']; 


        if showtags['xne'] and 'Nash' in show_opts:
            xtag = 'xne'; 
            templw = 3; templs = '--'; tempmsty = ''; tempmsz = 1.;
            temp = AA@BB@une + HH@x0; temp = temp.reshape([nt,nx]); temp = np.vstack([x0,temp])
            ind = 0; spat = spatials[ind]; xtrajsb[(xtag,ind)] = temp[:,spat];
            ind = 1; spat = spatials[ind]; xtrajsb[(xtag,ind)] = temp[:,spat];
            ind = 0; xpensb[(xtag,ind)] = PENS['xne'];#PEN({'ergba':[0,0,0,1],'lw':templw,'ls':templs,'msty':tempmsty,'msz':tempmsz})
            ind = 1; xpensb[(xtag,ind)] = PENS['xne'];#PEN({'ergba':[0,0,0,1],'lw':templw,'ls':templs,'msty':tempmsty,'msz':tempmsz})
            if not(PENS['xne1']==None): ind = 0; xpensb[(xtag,ind)] = PENS['xne1'];
            if not(PENS['xne2']==None): ind = 1; xpensb[(xtag,ind)] = PENS['xne2']; 

        if showtags['xna']: 
            xtag = 'xna';
            ### JUST A DOT... PRETTY WEIRD, but great
            temp = AA@BB@una + HH@x0; temp = temp.reshape([nt,nx]); 
            ind = 0; spat = spatials[ind]; xtrajsb[(xtag,ind)] = temp[:,spat];
            ind = 1; spat = spatials[ind]; xtrajsb[(xtag,ind)] = temp[:,spat];
            ind = 0; xpensb[(xtag,ind)] = PENS['xna'] #PEN({'ergba':[0,0,0,1],'lw':templw,'ls':templs,'msty':tempmsty,'msz':tempmsz})
            ind = 1; xpensb[(xtag,ind)] = PENS['xna'] #PEN({'ergba':[0,0,0,1],'lw':templw,'ls':templs,'msty':tempmsty,'msz':tempmsz})


        if showtags['approx']:
            xtag = 'approx'
            uu = une + toshow_coeffs2@toshow_Vi;
            temp = AA@BB@uu + HH@x0; temp = temp.reshape([nt,nx]); temp = np.vstack([x0,temp])
            templw = 2; templs = '--'; tempmsty = ''; tempmsz = 1.;
            ind = 0; spat = spatials[ind]; xtrajsb[(xtag,ind)] = temp[:,spat];
            ind = 1; spat = spatials[ind]; xtrajsb[(xtag,ind)] = temp[:,spat];
            ind = 0; xpensb[(xtag,ind)] = PENS['approx'] #PEN({'ergba':[1,0,0,1],'lw':templw,'ls':templs,'msty':tempmsty,'msz':tempmsz})
            ind = 1; xpensb[(xtag,ind)] = PENS['approx'] #PEN({'ergba':[0,0,1,1],'lw':templw,'ls':templs,'msty':tempmsty,'msz':tempmsz})
            if not(PENS['approx1']==None): ind = 0; xpensb[(xtag,ind)] = PENS['approx1']; #PEN({'ergba':[1,0,0,1],'lw':templw,'ls':templs,'msty':tempmsty,'msz':tempmsz})
            if not(PENS['approx2']==None): ind = 1; xpensb[(xtag,ind)] = PENS['approx2']; #PEN({'ergba':[0,0,1,1],'lw':templw,'ls':templs,'msty':tempmsty,'msz':tempmsz})


        if showtags['exact'] and 'SVO' in show_opts:
            xtag = 'exact'
            uu = une + coeffs2@Vi;
            temp = AA@BB@uu + HH@x0; temp = temp.reshape([nt,nx]); temp = np.vstack([x0,temp])
            templw = 4; templs = '-'; tempmsty = ''; tempmsz = 1.;
            ind = 0; spat = spatials[ind]; xtrajsb[(xtag,ind)] = temp[:,spat];
            ind = 1; spat = spatials[ind]; xtrajsb[(xtag,ind)] = temp[:,spat];
            ind = 0; xpensb[(xtag,ind)] = PENS['exact']; #PEN({'ergba':[1,0,0,1],'lw':templw,'ls':templs,'msty':tempmsty,'msz':tempmsz})
            ind = 1; xpensb[(xtag,ind)] = PENS['exact']; #PEN({'ergba':[0,0,1,1],'lw':templw,'ls':templs,'msty':tempmsty,'msz':tempmsz})
            if not(PENS['exact1']==None): ind = 0; xpensb[(xtag,ind)] = PENS['exact1']; #PEN({'ergba':[1,0,0,1],'lw':templw,'ls':templs,'msty':tempmsty,'msz':tempmsz})
            if not(PENS['exact2']==None): ind = 1; xpensb[(xtag,ind)] = PENS['exact2']; #PEN({'ergba':[0,0,1,1],'lw':templw,'ls':templs,'msty':tempmsty,'msz':tempmsz})




        for jj,toshow_vec in enumerate(toshow_Vi):

            if jj in basis_opts:
                uu = une + toshow_vec
                temp = AA@BB@uu + HH@x0; temp = temp.reshape([nt,nx]); temp = np.vstack([x0,temp])
                ind = 0; spat = spatials[ind]; xtrajs1[(jj,ind)] = temp[:,spat];
                ind = 1; spat = spatials[ind]; xtrajs1[(jj,ind)] = temp[:,spat];
                ind = 0; xpens1[(jj,ind)] = copy.deepcopy(PENS['basis']); #PEN({'ergba':[0,0,0,1],'lw':templw,'ls':templs,'msty':tempmsty,'msz':tempmsz})
                ind = 1; xpens1[(jj,ind)] = copy.deepcopy(PENS['basis']); #PEN({'ergba':[0,0,0,1],'lw':templw,'ls':templs,'msty':tempmsty,'msz':tempmsz})
                if not(PENS['basis1']==None): ind = 0; xpens1[(jj,ind)] = copy.deepcopy(PENS['basis1']);
                if not(PENS['basis2']==None): ind = 1; xpens1[(jj,ind)] = copy.deepcopy(PENS['basis2']); #PEN({'ergba':[0,0,0,1],'lw':templw,'ls':templs,'msty':tempmsty,'msz':tempmsz})
                if ins['use_basis_colormap']:
                    newcolor = ins['basis_colors'][jj]
                    # cmap = ins['basis_colormap']
                    # newcolor = cmap(jj+ins['basis_num_offset']/ins['basis_max_val'])
                    # if len(ins['basis_cmap_vals'])>0: newcolor = cmap(ins['basis_cmap_vals'][jj])
                    ind = 0; xpens1[(jj,ind)].update({'ec':newcolor[:3]})
                    ind = 1; xpens1[(jj,ind)].update({'ec':newcolor[:3]})



                uu = une + toshow_vec*toshow_coeffs2[jj]
                temp = AA@BB@uu + HH@x0; temp = temp.reshape([nt,nx]); temp = np.vstack([x0,temp])
                ind = 0; spat = spatials[ind]; xtrajs2[(jj,ind)] = temp[:,spat];
                ind = 1; spat = spatials[ind]; xtrajs2[(jj,ind)] = temp[:,spat];
                ind = 0; xpens2[(jj,ind)] = copy.deepcopy(PENS['basis_scaled']); #PEN({'ergba':[0,0,0,1],'lw':templw,'ls':templs,'msty':tempmsty,'msz':tempmsz})
                ind = 1; xpens2[(jj,ind)] = copy.deepcopy(PENS['basis_scaled']); #PEN({'ergba':[0,0,0,1],'lw':templw,'ls':templs,'msty':tempmsty,'msz':tempmsz})
                if ins['use_basis_colormap']:
                    newcolor = ins['basis_colors'][jj]
                    # cmap = ins['basis_colormap']
                    # newcolor = cmap(jj+ins['basis_num_offset']/ins['basis_max_val'])
                    # if len(ins['basis_cmap_vals'])>0: newcolor = cmap(ins['basis_cmap_vals'][jj])
                    ind = 0; xpens2[(jj,ind)].update({'ec':newcolor[:3]})
                    ind = 1; xpens2[(jj,ind)].update({'ec':newcolor[:3]})


            # for jj,spat in enumerate(spatials):
            #     xtempj = xtemp[:,spat]
            #     # if show_version == 'NL': xtempj = xtempj + xnl[:,spat];
            #     aid = spatial_agentids[jj]
            #     xtempj = np.vstack([x0[spat],xtempj])
            #     xtrajs[(phh,tt,aid,jj)] = xtempj; #np.array([x1temp])
            #     _,_,cc = curveColor2DGrid1(tag,ins=gridparams[wintag]) #,ins={}):
            #     temp = {'frgba':cc,'ergba':cc,'ls':'-','lw':2.,'zord':10}
            #     xpens[(phh,tt,aid,jj)] = PEN(temp)
        if False:         
            if showtags['use_time_grad'] == True:
                temp['useColorGrad'] = True; temp['curveColorGrad'] = ins['curveColorGrad'];
                templw = 4; templs = '--'; tempmsty = ''; tempmsz = 1.;
                if 'lw' in ins['curveColorGrad']: templw = ins['curveColorGrad']['lw'];
                if 'ls' in ins['curveColorGrad']: templs = ins['curveColorGrad']['ls'];
                if 'msty' in ins['curveColorGrad']: tempmsty = ins['curveColorGrad']['msty'];
                if 'msz' in ins['curveColorGrad']: tempmsz = ins['curveColorGrad']['msz'];
                temp['pen_main'] = PEN({'lw':templw,'ls':templs,'msty':tempmsty,'msz':tempmsz})

        temp = {'curves':xtrajsb,'bnds':bndpolys,'pens':xpensb,**axis4}
        STRUCTXsb = CURVES(ins=temp);
        STRUCTXsb.filterAndSegment(ftags['uframe'],eps=segeps);
        STRUCTXsb.drawCurves(ax,poly=ftags['uframe'],typ='segmented');

        if 'basis' in show_opts:
            if scale_opt == 'unscaled': #showtags['basis']:
                temp = {'curves':xtrajs1,'bnds':bndpolys,'pens':xpens2,**axis4}
                STRUCTXs1 = CURVES(ins=temp);
                STRUCTXs1.filterAndSegment(ftags['uframe'],eps=segeps);
                STRUCTXs1.drawCurves(ax,poly=ftags['uframe'],typ='segmented');
            if scale_opt == 'scaled': #showtags['basis_scaled']:
                temp = {'curves':xtrajs2,'bnds':bndpolys,'pens':xpens2,**axis4};
                STRUCTXs2 = CURVES(ins=temp);
                STRUCTXs2.filterAndSegment(ftags['uframe'],eps=segeps);
                STRUCTXs2.drawCurves(ax,poly=ftags['uframe'],typ='segmented');


    ###########################################################################################
    ###########################################################################################
    ###########################################################################################
    ###########################################################################################
    ###########################################################################################
    ###########################################################################################

    # ax3.set_aspect('equal')
    #########

    if showtags['trajs']:
        AXIS4 = AXIS(ins=params4);
        if showtags['axes'][3]: AXIS4.draw(ax)



    ###########################################################################################
    ###########################################################################################
    ###########################################################################################
    ###########################################################################################
    ###########################################################################################
    ###########################################################################################





    ###########################################################################################
    ###########################################################################################
    ###########################################################################################
    ###########################################################################################
    ###########################################################################################
    ###########################################################################################


    SVO1 = GAME1; dt = SVO1.LINQUAD.dt    # print(dt0)
    paramsIndAgents = GAME1.LINQUAD.paramsIndAgents;


    #siginds = [0,2,2,3,4,5,6,7,8,9,10,11]
    siginds = [0,2,4,6,8,10]
    # siginds = [1,3,5,7,9,11];

    basis_opt_max = np.max(basis_opts_even)

    slices = {};
    for sigind in siginds:
        # sigind = 8;
        # print(paramsIndAgents)
        nx = 2;
        x1des = SVO1.LINQUAD.x1des[:,:nx];
        x2des = SVO1.LINQUAD.x2des[:,nx:];
        x10 = x1des[0];
        x20 = x2des[0];
        nt = SVO1.LINQUAD.nt;

        x1des_flat = x1des.reshape(nt*nx);
        x2des_flat = x2des.reshape(nt*nx);

        Q12 = np.eye(nt*2); HH = np.vstack([np.eye(2) for _ in range(nt)]);
        Fk = np.eye(2); Gk = dt*np.diag([1.01,0.99])
        Fbar = np.tril(np.ones(nt))
        FF = np.kron(Fbar,Fk); GG = np.kron(np.eye(nt),Gk);
        X,sigs,YT = mat.svd(Q12@FF@GG);
        Xl = X[:,sigind]; Yl = YT[sigind]; sigl = sigs[sigind];

        r1 = paramsIndAgents['r1']; r2 = paramsIndAgents['r2']
        q1 = paramsIndAgents['q1']; q2 = paramsIndAgents['q2']
        r1b = paramsIndAgents['r1b']; r2b = paramsIndAgents['r2b']
        q1b = paramsIndAgents['q1b']; q2b = paramsIndAgents['q2b']
        p1 = paramsIndAgents['p1']; p2 = paramsIndAgents['p2']

        zetap1 = p1/(p1+p2); zetap2 = p2/(p1+p2);
        zetar1 = r1/(p1+p2); zetar2 = r2/(p1+p2);
        zetaq1 = q1/(p1+p2); zetaq2 = q2/(p1+p2);

        alpha1l = r1 + (q1 - p1)*sigl*sigl;
        alpha2l = r2 + (q2 - p2)*sigl*sigl;
        beta1l = p1*sigl*sigl;
        beta2l = p2*sigl*sigl;
        delta1l = r1b + (q1b - p1)*sigl*sigl;
        delta2l = r2b + (q2b - p2)*sigl*sigl;

        DelMl = alpha1l*alpha2l - beta1l*beta2l;
        DelNl = delta1l*delta2l - beta1l*beta2l;
        Del1l = alpha1l*delta1l - beta1l*beta2l;
        Del2l = alpha2l*delta2l - beta1l*beta2l;

        gam1l = np.abs(Del1l)/np.abs(DelNl);
        gam2l = np.abs(Del2l)/np.abs(DelNl);
        gamMl = np.abs(DelMl)/np.abs(DelNl);

        barr1l = r1/(r1+q1*sigl*sigl); barr2l = r2/(r2+q2*sigl*sigl);
        barq1l = q1/(r1+q1*sigl*sigl); barq2l = q2/(r2+q2*sigl*sigl);
        barp1l = p1/(r1+q1*sigl*sigl); barp2l = p2/(r2+q2*sigl*sigl);

        barr1lb = r1b/(r1b+q1b*sigl*sigl); barr2lb = r2b/(r2b+q2b*sigl*sigl);
        barq1lb = q1b/(r1b+q1b*sigl*sigl); barq2lb = q2b/(r2b+q2b*sigl*sigl);
        barp1lb = p1/(r1b+q1b*sigl*sigl); barp2lb = p2/(r2b+q2b*sigl*sigl);

        barDelMl = 1 - barp1l - barp2l;
        barDelNl = 1 - barp1lb - barp2lb;

        barxi1l = Xl@Q12@x1des_flat;
        barxi2l = Xl@Q12@x2des_flat;
        xi10l = Xl@Q12@FF@HH@x10;
        xi20l = Xl@Q12@FF@HH@x20;

        h1l = ((1-barp2l)/(sigl*barDelMl))*barq1l*(barxi1l - xi10l)
        h2l = ((1-barp1l)/(sigl*barDelMl))*barq2l*(barxi2l - xi20l)
        h1l = h1l - (barp1l/(sigl*barDelMl))*(barq2l*barxi2l+barr2l*xi20l-xi10l)
        h2l = h2l - (barp2l/(sigl*barDelMl))*(barq1l*barxi1l+barr1l*xi10l-xi20l)

        h1lb = ((1-barp2lb)/(sigl*barDelNl))*barq1lb*(barxi1l - xi10l)
        h2lb = ((1-barp1lb)/(sigl*barDelNl))*barq2lb*(barxi2l - xi20l)
        h1lb = h1lb - (barp1lb/(sigl*barDelNl))*(barq2lb*barxi2l+barr2lb*xi20l-xi10l)
        h2lb = h2lb - (barp2lb/(sigl*barDelNl))*(barq1lb*barxi1l+barr1lb*xi10l-xi20l)


        # def g1calc(th1,th2):
        #     return 1 + (1./np.tan(th1))*zetap1 + (1./np.tan(th2)) * (zetap2 - (1/(sigl*sigl))*zetar2 - zetaq2)
        # def g2calc(th1,th2):
        #     return 1 + (1./np.tan(th2))*zetap2 + (1./np.tan(th1)) * (zetap1 - (1/(sigl*sigl))*zetar1 - zetaq1)


        def g1lcalc(phi,sgn=1):
            term1 = Del1l/(2*np.cos(phi)) - Del2l/(2*np.sin(phi))
            temp2 = Del1l/(2*np.cos(phi)) + Del2l/(2*np.sin(phi))
            term2 = np.sqrt(temp2**2 - DelMl*DelNl/(np.cos(phi)*np.sin(phi)))
            term3 = (delta1l + beta1l)*beta1l/np.cos(phi) - (alpha2l+beta2l)*beta2l/np.sin(phi)
            return term1 + sgn*term2 + term3

        def g2lcalc(phi,sgn=1):
            term1 = Del2l/(2*np.sin(phi)) - Del1l/(2*np.cos(phi))
            temp2 = Del2l/(2*np.sin(phi)) + Del1l/(2*np.cos(phi))
            term2 = np.sqrt(temp2**2 - DelMl*DelNl/(np.cos(phi)*np.sin(phi)))
            term3 = (delta2l + beta2l)*beta2l/np.sin(phi) - (alpha1l+beta1l)*beta1l/np.cos(phi)
            return term1 + sgn*term2 + term3

        def lamphil(phi,sgn=1):
            term1 = Del1l/(2*np.cos(phi)) + Del2l/(2*np.sin(phi))
            term2 = np.sqrt(term1**2 - DelMl*DelNl/(np.cos(phi)*np.sin(phi)))
            return (1/DelNl)*(term1 + sgn*term2)

        def Thetainv(theta):    
            th1 = theta[0]; th2 = theta[1];
            phi = np.arctan2(np.tan(th2),np.tan(th1))
            s = np.tan(th2)/np.sin(phi)
            return phi,s

        def flarg(theta,maxval=100,useabs=True,eigsonly=True):
            phi,s = Thetainv(theta)
            g1p = g1lcalc(phi,sgn= 1); g2p = g2lcalc(phi,sgn= 1); 
            g1m = g1lcalc(phi,sgn=-1); g2m = g2lcalc(phi,sgn=-1);
            normplus = mat.norm([g1p,g2p]);
            if False: norm = mat.norm([g1p*g2m,g2p*g1m]);
            if False: norm = 1; 
            if True: norm = 10*(g1p*g2m - g2p*g1m)/normplus
            lamphi = lamphil(phi)
            term2 = 1./(1+lamphi/s)
            term1 = (g2m/norm)*(h1lb - h1l)-(g1m/norm)*(h2lb - h2l)
            if eigsonly: out = term2;        
            else: out = term2*term1;
            if useabs: out = np.abs(out)
            if out > maxval: out = maxval
            # print(out)
            return out
        
        plotbuff = 0.001
        maxval = 50+plotbuff; buff = 0.0001; theta2 = theta2_opt; eigsonly = True;
        theta1s = np.linspace(buff,np.pi/2-buff,1000)
        slices[sigind] = np.array([flarg((theta1,theta2),maxval=maxval,eigsonly=eigsonly) for theta1 in theta1s])

        if sigind == basis_opt_max:
            nth1 = 200; nth2 = 200; buff = 0.001; maxval = 100.
            th1s = np.linspace(buff, np.pi/2-buff, nth1)
            th2s = np.linspace(buff, np.pi/2-buff, nth2)
            TH1,TH2 = np.meshgrid(th1s,th2s)
            ZZ = np.zeros([nth1,nth2]);
            IM = np.zeros([nth1,nth2]);
            for i,th1 in enumerate(th1s):
                for j,th2 in enumerate(th2s):
                    val = flarg([th1,th2],maxval=maxval,useabs=True,eigsonly=eigsonly);
                    IM[i,j] = val;
                    ZZ[i,j] = val;
            if False: IM = np.sqrt(np.sqrt(IM))
            IM = IM.T[::-1]
            cmap = parula_colormap
            basis_cmap_vals = [0.1,0.1,0.35,0.35,0.5,0.5,0.6,0.6,0.85,0.85,0.99,0.99];
            maincolor = list(cmap(basis_cmap_vals[sigind]));
            ntemp = 100;
            custom_colors = [maincolor[:3] + [i/ntemp] for i in range(ntemp)]
            imcmap = ListedColormap(custom_colors)
            levels = np.linspace(1.1,maxval,int(maxval*1)); #3.*np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]);
            #levels = np.linspace(1.1,1000,20); #int(maxval*1)); #3.*np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]);
            levels = np.hstack([(1./levels)[::-10],levels]);
            ZZ0 = ZZ/np.max(ZZ);
            ZZ0[np.where(ZZ0>0.99)] = 1.
            im = ax2.imshow(IM, cmap=imcmap,extent=[th1s[0], th1s[-1], th2s[0], th2s[-1]]); #, origin='lower', cmap='viridis', aspect='auto')
            CS = ax2.contour(TH1, TH2, ZZ, levels=levels, colors=[maincolor[:3]+[1.]], linewidths=1)
            CS = ax2.contour(TH1, TH2, ZZ0, levels=[0.999], colors=[[0,0,0,1.]], linewidths=2,zorder=300)    




    cmap = parula_colormap
    basis_cmap_vals = [0.1,0.1,0.35,0.35,0.5,0.5,0.6,0.6,0.85,0.85,0.99,0.99];
    basis_colors = [cmap(basis_cmap_vals[j]) for j in siginds];
    for j,sigind in enumerate(slices):
        if sigind in basis_opts_even:
            ax3.plot(theta1s,slices[sigind],color=basis_colors[j],lw=4)
        else: ax3.plot(theta1s,slices[sigind],color=basis_colors[j],alpha=0.5,lw=2)

    # ax.plot(theta1s,np.ones(len(theta1s)),color=[0,0,0,0.8],linestyle='--',lw=2)
    ax3.tick_params(axis='both', labelsize=32)  # Set both x and y tick label fontsize to 20
    ax3.set_xlim([buff,np.pi/2-buff])
    # ax3.set_xticks([])
    ax3.set_ylim([0,50]);#maxval-2*plotbuff])





    ###########################################################################################
    ###########################################################################################
    ###########################################################################################
    ###########################################################################################
    ###########################################################################################
    ###########################################################################################





    ax.set_aspect('equal')
    # ax.axis('off')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim([-2.5,3.5])
    ax.set_ylim([-2.5,3.5])
    # ax.set_xlim([-2.5,2.5])
    # ax.set_ylim([-2.5,2.5])
    # ax.set_xlim=


    fntsz1 = 22;
    ax2.set_aspect('equal')
    
    ax2.set_ylabel(r'$\theta_2$',fontsize=fntsz1+6,rotation=0,labelpad=-20)
    ax2.set_xlabel(r'$\theta_1$',fontsize=fntsz1+6,labelpad=-20)
    ax2.set_xticks([0,np.pi/2],['0',r'$\pi/2$'],fontsize=fntsz1)
    theta2_str = str(np.round(theta2_opt,2))
    ax2.set_yticks([0,theta2_opt,np.pi/2],['0',theta2_str,r'$\pi/2$'],fontsize=fntsz1)
    ax2.set_xlim([0,np.pi/2])
    ax2.set_ylim([0,np.pi/2])
    ax2.axhline(y=theta2_opt, color='black', linestyle='--',linewidth=1)    


    ax3.set_xlabel(r'$\theta_1$',fontsize=fntsz1+6)

    title3 = r'$|w_{1\ell}(\phi)/(1+\lambda_{\ell+}^\phi/s)| \ \ \  $  $\theta_2 = \ $' + theta2_str
    ax3.set_title(title3,fontsize=fntsz1,pad=10)


    ylabel1 = r'$|\frac{w_{1\ell}(\phi)}{(1+\lambda_{\ell+}^\phi/s)}|$'


    ax1.set_ylim([0,2.5]);    
    ax1.annotate(ylabel1, xy=(9, 1.2),fontsize=fntsz1 + 10)

    # ax1.set_ylabel(ylabel1,fontsize=fntsz1+6,rotation=0,labelpad=20)



