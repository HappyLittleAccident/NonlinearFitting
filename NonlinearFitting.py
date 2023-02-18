# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.optimize import curve_fit
import os
import lmfit

#%% 
def PhaseChange(xs,ys,phase):
    """
    Change phase of complex data according to
    
    xs_new + 1j*ys_new = (xs + 1j*ys)*np.exp(1j*phase)
    
    Parameters
    
    ---------
    
    xs, ys : input data, real
    
    phase : angle to shift phase by, in radians
    
    Returns
    
    -------
    
    xs_new, ys_new : data with changed phase
    """
    zs = np.copy(xs + 1j*ys)
    xs_new = np.real(zs*np.exp(1j*phase))
    ys_new = np.imag(zs*np.exp(1j*phase))
    return xs_new, ys_new

#%%
def BgndRoughCorr(xs,ys,show:bool=False):
    """
    Corrects constant background noise
    """
    xshift = xs[-1]
    yshift = (ys[0]+ys[-1])/2
    xs_new = xs - xshift
    ys_new = ys - yshift
    if show==True:
        print('bgnd was changed by\n'
              'Xshift = ', xshift, '\nYshift = ', yshift, '\n')
    return xs_new, ys_new, xshift, yshift
    
#%%
def PhaseRoughCorr(xs,ys,show:bool=False):
    """
    Rough phase correction for reasonance curve, based on projection in x-y plane
    """
    phase = np.angle(ys[0]-ys[-1] + 1j*(xs[0]-xs[-1])) 
    xs_new, ys_new = PhaseChange(xs,ys,phase)
    if show==True:
        print('Phase was changed by\n'
              'Phi_shift = ', phase, '\n')
    return xs_new, ys_new, phase

#%%
def GetWeightsDefault(xs,fs,f0,width):
    """
    Gives artifficial errors to data according, emphasizing peak data and 
    deprecating edge data
    
    for fiting purposes
    """
    weights = 1/(((fs-f0)/width)**2 +1) + 1
    return weights

#%%
def GetWeightsBackground(xs,treshold):
    """
    Creates weights passed to EdgeFit procedure, eliminating peak data above "treshold" from fit
    
    --------------
    
    xs : array
            absorption peak data
            
    
    treshold : real from interval [0,1]
                determine relative height above which datapoints are excluded from fit
    """
    if treshold > 1 or treshold <0:
        raise Exception('treshold variable must be within (0,1)')
        

    weights = np.where(xs/np.max(xs) > treshold,np.full_like(xs,np.inf),np.ones_like(xs))
    return weights

#%%
def ParameterEstimate(fs,xs,ys):
    """
    Initial parameters for linear resonance curve with linear background.
    """

    In_f0 = fs[np.argmax(xs)]       #resonance freq
    
    halfmax = np.max(xs)/2          #height at half maximum
    test_xs = xs[3:np.argmax(xs)]  
    test_ws = fs[3:np.argmax(xs)]
    w1 = np.interp(halfmax, test_xs, test_ws)
    
    test_xs = xs[-3:np.argmax(xs):-1]
    test_ws = fs[-3:np.argmax(xs):-1]
    w2 = np.interp(halfmax, test_xs, test_ws)
    
    In_w = np.abs(w2-w1)/np.sqrt(3)
    
    In_amp = np.max(xs)
    In_lor_amp = In_amp*In_w
    
    return [In_f0,In_w,0,0,0,0,In_lor_amp]

#%%
def ResonanceLin(fs,f0,w,a,b,c,d,A):
    """
    real and imaginary part for velocity linear resonance curve, with linear background given by:
        zs = 1j*fs/(f0**2 - fs**2 + 1j*fs*w)*np.exp(1j*Phi)
        Real = A*np.real(zs)+a+c*(fs-f0)
        Imaginary = A*np.imag(zs)+b+d*(fs-f0)
        
    Parameters
    
    ----------
    
    fs : frequency
    
    f0 : reasonance frequency
    
    w : width
    
    Phi : complex phase between real and imaginary part
    
    a : constant background in real part
    
    b : constant background in imaginary part
    
    c : linear background in real part
    
    d : linear background in imaginary part
    
    A : amplitude of curve
    
    Returns
    
    -------
    
    np.array([ Real , Imaginary ]) where length of Real and Imaginary matches length of fs
    """
    zs = 1j*fs/(f0**2 - fs**2 + 1j*fs*w)     
    return np.ravel([A*np.real(zs)+a+c*(fs-f0), A*np.imag(zs)+b+d*(fs-f0)])


#%% EdgeFit
def EdgeFit(f,x,y,treshold,show=False):
    """
    Takes in phase corrected resonance data, 
    estimates best linear reasonance
    curve according to edges of said curve
    and corrects absorption x and dispersion y
    by estimated linear background
    
    Returns
    
    -------
    
    x, y : corrected absorption and dispersion respectively
    
    params : [f0, width, Absorbtion const bgnd, Dispersion const bgnd, Abs lin, Disp lin, Amplitude]
    """

    x,y,xshift,yshift = BgndRoughCorr(x, y)
    if show:
        fig,ax = plt.subplots()
        ax.plot(f,x,'k+')
        ax.plot(f,y,'b+')

    weights = GetWeightsBackground(x, treshold)

    init_par = ParameterEstimate(f, x, y)

    if show:
        ax.plot(f,np.split(ResonanceLin(f,*init_par),2)[0],'r--',label='init guess absorption')
        ax.plot(f,np.split(ResonanceLin(f,*init_par),2)[1],'g--',label='init guess dispersion')
    params,_ = curve_fit(ResonanceLin,f, np.ravel([x,y]), init_par, sigma=np.ravel([weights,weights]))

    bgnd = params[2:6]
    
    x-= bgnd[2]*(f-params[0]) + bgnd[0]
    y-= bgnd[3]*(f-params[0]) + bgnd[1]
    
    if show:
        ax.plot(f,np.split(ResonanceLin(f,*params),2)[0])
        ax.plot(f,np.split(ResonanceLin(f,*params),2)[1])

    
    return x,y,params


#%% Correct data for phase, lin background
def Corrections(fq,x,y):
    """
    Corrects absorption x, dispersion y using "RoughPhaseCorr" and "EdgeFit"
    
    Returns
    
    -------
    
    x, y : corrected absorption and dispersion
    
    z : corrected overall amplitude z = np.sqrt(x**2 + y**2)
    
    params_edge : fit parameters from "EdgeFit" in list of form \n
    [f0, g1, x_const_bgnd, y_const_bgnd, x_lin_bgnd, y_lin_bgnd, amplitude]
    """
    x,y,phase = PhaseRoughCorr(x,y)
    x,y,params_edge = EdgeFit(fq,x,y,0.1,show=False)
    z = np.sqrt(x**2 + y**2)
    return x,y,z,params_edge

#%%

class ResonanceFit:
    """
    Class with minimizing method, and data output methods. Does not work alone, 
    new child class containing "model" and "params_estimate" methods needs to be created.
    
    -------------------
    
    fqs : list of numpy arrays \n
        each array corresponds to measured frequencies from one peak
    
    xs, ys : same type as fqs \n
        corresponds to absorption and dispersion
    
    drives : list of numbers or numpy array \n
        contains drive info for each peak
        """
    run_fit = False
    def __init__(self,fqs=None,xs=None,ys=None,drives=None,filenames = None):
        """
        fqs : list of numpy arrays \n
            each array corresponds to measured frequencies from one peak
        
        xs, ys : same type as fqs \n
            corresponds to absorption and dispersion
        
        drives : list of numbers or numpy array \n
            contains drive info for each peak
        """
        if type(self).__name__ == 'ResonanceFit':
            raise Exception('Create new class with "model" and "params_estimate" functions defined!')
        
        if xs == None or ys == None or fqs == None:
            raise Exception('Required arguments fqs, xs or ys not supplied!')
        
        if np.array_equal(drives,None):
            drives = [None]*len(xs)
        
        if np.array_equal(filenames,None):
            filenames = [None]*len(xs)
            
        if not (len(drives) == len(filenames) == len(fqs) == len(xs) == len(ys)):
            raise Exception('All suplied variables fqs, xs, ys, drives and filenames must have the same length!')
        
        if filenames[0] == None:
            filenames == None
            
        if drives[0] == None:
            drives = None
        
        if np.array_equal(drives,None) and np.array_equal(filenames,None):
            drives,fqs,xs,ys = map(zip(*sorted(zip(drives,fqs,xs,ys))) )
            
        elif np.array_equal(drives,None):
            drives,fqs,xs,ys,filenames = map(list,zip(*sorted(zip(drives,fqs,xs,ys,filenames))) )
        
        
        self.fqs = fqs
        self.xs = xs
        self.ys = ys        
        self.zs = None
        self.drives = drives
        self.forwards = None
        self.GetWeights = None
        self.results = None
        self.init_pars = None
        self.filenames = filenames
        

            

        

    def model(self, params, fq, z):
        return None
    
    def model_weighted(self,params, fq, z, weights):
        return np.abs(self.model(params, fq, z))*weights
    
    def params_estimate(self):
        return None
    
    def remove_transient(self,show=False):
        """
        Removes measured point too far away from other measurements, 
        useful for data exhibiting duffing nonlinearities.
        
        show : optional parameter \n
            creates graph for each peak, showing removed datapoints
        """
        for i,x in enumerate(self.xs):
            x = self.xs[i]
            i_delete = []
            for j in range(len(x)-4):
                dist1 = np.abs(x[j+1]-x[j+2])
                dist2 = np.abs(x[j+3]-x[j+2])
                dist_min = np.min([dist1,dist2])
                dist_max = np.max([dist1,dist2])
                
                min_condition=dist_min > (np.max(x)-np.min(x))*0.05
                max_condition=dist_max > (np.max(x)-np.min(x))/5
                
                if min_condition and max_condition and np.max(x) != x[j+2]:
                    i_delete.append(j+2)
            if show:
                plt.figure()
                plt.plot(self.fqs[i],x,'r+')
                
            self.xs[i] = np.delete(x,i_delete)
            self.ys[i] = np.delete(self.ys[i],i_delete) 
            self.fqs[i] = np.delete(self.fqs[i],i_delete) 
            if self.zs!= None:
                self.zs[i] = np.delete(self.zs[i],i_delete) 
            
            if show:
                plt.plot(self.fqs[i],self.xs[i],'k+')

        
    def fit(self,vary_pars:dict={},init_values:dict={},bounds:dict={}):
        """
        Correct supplied data xs and ys, and find optimal parameters for all supplied peaks
        for given "model" and "params_estimate".
        
        Data xs and ys are replaced by their corrections, fit results are stored in class property "results",
        which is a list of lmfit.MinimizerResult objects corresponding to fitted peaks.
        
        Initial parameters are stored in properties "init_pars", which is a list of lmfit.params objects.
        
        Optional parameters
        
        -------
        
        vary_pars : dictionary {'name': True/False, .. },\n
                     given parameters will be forcibly varied/fixed
        
        init_values : dictionary {'name': init_value, .. },\n
                    override estimate for given parameters
        
        bounds : dictionary {'name': (low,high), .. },\n
                     override initial bounds for given parameters
        """


            
        self.forwards = np.ones(len(self.xs),dtype=bool)
        self.zs = [None for x in self.xs]
        self.results = [None for x in self.xs]
        self.init_pars = [None for x in self.xs]
        
        for i in range(len(self.fqs)):
            fq = self.fqs[i]*2*np.pi
            x = self.xs[i]
            y = self.ys[i]
            
            if fq[0] > fq[-1]:
                forward = False
                fq = fq[::-1]
                x = x[::-1]
                y = y[::-1]
            else:
                forward = True
            self.forwards[i] = forward
            
            self.fqs[i]=fq
            try:
                x,y,z,params_edge = Corrections(fq,x,y)
        
        
                self.xs[i] = x
                self.ys[i] = y
                self.zs[i] = z
                
                fq0 = params_edge[0]
                g1_edge = params_edge[1]
                
                init_params = self.params_estimate(params_edge,i)
                
                params = lmfit.Parameters()
                for param in init_params:
                    if type(param) == tuple:
                        params.add(param[0], value = param[1], vary = param[2], min = param[3], max = param[4])
                    else:
                        raise Exception('"init_params" must be list of tuples')
                
                if type(vary_pars) != dict:
                    if type(vary_pars) == bool:
                        vary_pars = {par_name:vary_pars for (par_name, _ ) in params.items()}
                        
                    else:
                        raise Exception('"vary_pars" must be dictionary or boolean')
                
                if type(init_values) != dict:
                    raise Exception('"init_values" must be dictionary')
                    
                if type(bounds) != dict:
                    raise Exception('"bounds" must be dictionary')
                    
                for key in vary_pars:
                    params[key].vary = vary_pars[key]
                
                for key in init_values:
                    params[key].value = init_values[key]
                
                for key in bounds:
                    params[key].min = bounds[key][0]
                    params[key].min = bounds[key][1]
                
                self.init_pars[i]=params
                
                if self.GetWeights == None:
                    weights = GetWeightsDefault(z,fq,fq0,g1_edge)
                else:
                    weights = self.GetWeights(z,fq,init_params)
               
                self.results[i] = lmfit.minimize(self.model_weighted, params, args=(fq, z, weights),gtol=1e-20)
            except:
                print('fit number {} did not converge, or numpy.interpolate crashed'.format(i+1))
                self.zs[i] = np.sqrt(x**2 + y**2)
            
            self.fqs[i] = fq/ 2/np.pi
        self.run_fit = True
            
    def figures(self,xlabel:str='Driving frequency (Hz)',ylabel:str='Speed amplitude (m/s)'
                    ,extension:str='.png',names='',show:bool=True,path=None,save=False,plot_estimate=True):
        """
        Draw or save figures with measured data and fit curves. 
        Allways draws a graph containing all measured peaks.
        
        Parameters
        
        ---------
        
        xlabel, ylabel : string, axis labels
        
        extension : string, which extension to save figure with, example '.png, ...'
        
        names : string or list of strings. If string is given, figure names will be name1, name2, ...\n
                If list is given, each list element corresponds to figure name with same order as fqs, xs are given.
                
        show : bool or list of bools\n
                If True, show all figures. \n
                If list is given, it must be the same length as fqs. Each element states if figure is to be shown or not,
                in same order as fqs are given.
                
        path : if None, figures are not saved. If string is given, figures are saved to given path.
        
        """
        if not self.run_fit:
            raise Exception('"fit" method has not been executed')
            
        
        plt.ioff()
        
        if type(show) == bool:
            shows = [show for x in self.xs]
        else:
            shows = show
        
        if type(save) == bool:
            saves = [save for x in self.xs]
        else:
            saves = save
            
        if type(plot_estimate) == bool:
            plot_estimates = [plot_estimate for x in self.xs]
        else:
            plot_estimates = plot_estimate
        
        for i,result in enumerate(self.results):
            if result == None:
                fq_plot = np.linspace(self.fqs[i][0],self.fqs[i][-1],800)
                x_plot = np.linspace(0,np.max(self.zs[i])+0.6*np.max(self.zs[i]),800)
                FQ,X = np.meshgrid(fq_plot,x_plot)
                
                fig,ax = plt.subplots()
                datapoints = ax.plot(self.fqs[i],self.zs[i],'k+',label=r'Amplitude $z = \sqrt{x^2 + y^2}$')

                if plot_estimates[i]:
                    ax.contour(FQ,X,self.model(self.init_pars[i],FQ*2*np.pi,X),levels = [0] ,colors = 'r',linestyles='dashed')
                    red_line = Line2D([],[],color='red',linestyle='dashed', label='Initial estimate')
                    ax.legend(handles = [red_line,datapoints[0]])
            elif shows[i] or (saves[i] and type(path)==str):
                #get values for countour plot
                fq_plot = np.linspace(self.fqs[i][0],self.fqs[i][-1],800)
                x_plot = np.linspace(0,np.max(self.zs[i])+0.6*np.max(self.zs[i]),800)
                FQ,X = np.meshgrid(fq_plot,x_plot)
                
                #plot datapoints
                fig,ax = plt.subplots()
                datapoints = ax.plot(self.fqs[i],self.zs[i],'k+',label=r'Amplitude $z = \sqrt{x^2 + y^2}$')
                
                #plot fit curve
                ax.contour(FQ,X,self.model(result.params,FQ*2*np.pi,X),levels = [0] ,colors = 'b')
                blue_line = Line2D([],[],color='blue', label='Fit curve')
                #plot estimate
                if plot_estimates[i]:
                    ax.contour(FQ,X,self.model(self.init_pars[i],FQ*2*np.pi,X),levels = [0] ,colors = 'r',linestyles='dashed')
                    red_line = Line2D([],[],color='red',linestyle='dashed', label='Initial estimate')
                    ax.legend(handles = [red_line,blue_line,datapoints[0]])
                else:
                    ax.legend(handles = [blue_line,datapoints[0]])

                fig.canvas.manager.set_window_title('Number: {:d}'.format(i+1))
                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel)
                

                if type(names) == str:
                    fig.savefig(os.path.join(path,names+str(i+1)+extension))
                elif len(names) != len(self.results):
                    raise Exception('Len of "names" does not match with len of "results"')
                else:
                    fig.savefig(os.path.join(path,names[i]+extension))
                

                if shows[i]:
                    fig.show()
                else:
                    plt.close(fig)

        
        fig_all,ax_all = plt.subplots()
            #show all measurements in one graph
        for fq,z in zip(self.fqs,self.zs):
            ax_all.plot(fq,z,'--')
        ax_all.set_xlabel(xlabel)
        ax_all.set_ylabel(ylabel)
        fig_all.show()
            
        plt.ion()
            
    def save_params(self,folder:str='',name:str='params',which_peaks=None):
        """
        Saves fit parameters and other important data into .txt file, each line corresponding to one peak.\n
        
        Saved data have the format:
            
        ---------------
            
        #forward, drive, f_max, z_max, parameters, parameter errors, redchi
        
        forward : 1 if sweep was forward, 0 if sweep was backward
        
        drive : supplied driving parameter
        
        f_max : frequency for maximal amplitude in Hz
        
        z_max : maximal amplitude
        
        parameters : model parameters in order supplied
        
        parameter errors : in same order as parameters
        
        redchi : reduced chi squared (sum of residuals) of the fit
        
        Parameters
        
        -------
        
        folder : string \n
            to which folder is a file saved
        
        name : name of the file\n
            without extension
        
        which_peaks : optional, array of boolean \n
            give specific peaks to save
        """
        
        if self.drives == None:
            self.drives = [0 for x in self.fqs]
        
        header = 'forward\tdrive\tfq_max\tz_max'    #make header

        for key in self.results[0].params:  #add params
            header+='\t'+key
        for key in self.results[0].params:  #add errors
            header+='\t'+key+'_err'
        header+='\tredchi'                  #second last is reduced chi squared
        header+='\tfname'                   #last is the name of file
        
        if np.array_equal(which_peaks,None):
            which_peaks = [True for x in self.xs]
        
        if np.array_equal(self.filenames,None):
            self.filenames = [None for x in self.xs]
        
        lines_to_save = [] #prepare output datafile
        j=0
        for result,drive,fq,z,forward,this_peak,fname in zip(self.results,self.drives, self.fqs,
                                                             self.zs,self.forwards,which_peaks,self.filenames):
            if result == None:
                lines_to_save.append(['Did not converge'])
            elif this_peak:
                num_of_params = len(result.params)
                line = ['']*(3 + 2*num_of_params + 2)  #each line corresponds to one fit
                
                line[0] = str(int(forward))
                line[1] = '{:e}'.format(drive)
                line[2] = '{:e}'.format(fq[np.argmax(z)])
                line[3] = '{:e}'.format(np.max(z))
                for i,key in enumerate(result.params):
                    line[4+i]='{:e}'.format(result.params[key].value)
                    if result.params[key].stderr == None:
                        line[4+ num_of_params +i] = '{:f}'.format(np.nan)
                    else:
                        line[4+ num_of_params +i] = '{:f}'.format(result.params[key].stderr)
                line[-2]='{:e}'.format(result.redchi)
                line[-1]='{}'.format(fname)
                
                lines_to_save.append(line)
                j+=1
            
        data_to_save = np.stack(lines_to_save) #create numpy 2D array
            
        location = os.path.join(folder,name+'.txt') 
        np.savetxt(fname=location,X=data_to_save,header=header,fmt = '%s')  #save file

    def save_reports(self,folder:str='',name:str='report',which_peaks=None):
        """
        Saves fit reports from lmfit, each peak goes to one .txt file.
        
        Parameters
        
        ------
        
        folder : string \n
                path to folder
        
        name : string or list of strings\n
                if string, results are numbered name1, name2, ...
                
        which_peaks : optional, array of bool\n
                each element states if given peak is saved or not
        """
        if np.array_equal(which_peaks,None):
            which_peaks = [True for x in self.xs]
        for i,result in enumerate(self.results):
            if result == None:
                result = 'The fit did not converge'
            if which_peaks[i]:
                if type(name) == str:
                    path = os.path.join(folder,name+str(i+1)+'.txt')
                elif type(name) == list:
                    path = os.path.join(folder,name[i]+'.txt')
                else:
                    raise Exception('"name" must be either string or list of strings') 
                file = open(path,'w')
                file.write(lmfit.fit_report(result))
                file.close()
                

                        

class HelsinkyFit(ResonanceFit):
    def model(self,params,fq, z):
        """
        Implicit equation for speed reasonance curve of nonlinear oscillator, F(fq,z) = 0
        
        Parameters
        ----------
        params : result.params object from lmfit package \n
        fq : frequency \n
        z : amplitude \n
    
        Output
        -----
        Numpy compatible output with regards to fq, z. Implicit equation for reasonance curve, F(fq,z) = 0
        
        Equation specification
        ----------------------
        params : [a1,a2,a3,g1,g2,g3,f]
        
        **Restoring force**
        
        a1 : alpha1 \n
        a3 : alpha2 \n
        
        **Damping**
        
        g1 : gamma1 \n
        g2 : gamma2 \n
        g3 : gamma2 \n
        
        **Forcing**
        
        f : f (forcing amplitude)
        
        **Independent variables**
        
        fq : omega (forcing frequency in rad/s) \n
        z : z (amplitude)
        
        """
        
        a1, a3, g1, g2, g3, f = params.values()
        
        return ((z/fq)**2 *( (a1  + 0.75*a3*(z/fq)**2 - fq**2 )**2
                + fq**2 * (g1 + ( 8/(3*np.pi) )*g2*z + 0.75*g3*z**2)**2 ) - f**2)

    
    def params_estimate(self,par_edge,i:int):
        z = self.zs[i]
        fq = self.fqs[i]
        fq0=par_edge[0]
        fq_max = fq[np.argmax(z)]
        z_max = np.max(z)
        z_half = z_max/2
        if not np.array_equal(self.drives, None):
            drive=self.drives[i]
        #    #name #init #vary #min #max
        a1 = ('a1',fq0**2,True, 0,   np.inf)

        a3_guess = 4/3*(fq_max**2 - fq0**2)*fq_max**2/z_max**2
        a3= ('a3',a3_guess,True, -np.inf,np.inf)
        
        z_left = z[3:np.argmax(z)]
        fq_left = fq[3:np.argmax(z)]
        fq_l = np.interp(z_half, z_left, fq_left)
        
        z_right = z[-3:np.argmax(z):-1]
        fq_right = fq[-3:np.argmax(z):-1]
        fq_r = np.interp(z_half, z_right, fq_right)
        
        width = np.abs(fq_r-fq_l)/np.sqrt(3)
        if i==0:
            g1_guess = width
        else:
            g1_guess = self.init_pars[0]['g1'].value
        g1 = ('g1',g1_guess,True,0,np.inf)
        
        
        if self.drives != None:
            HWD = np.max(self.zs[0])*g1_guess/self.drives[0]
        else:
            HWD = None
        
        if HWD == None:
            g2_guess = 0.5*(width - g1_guess)/z_max
            g3_guess = 0.5*(width - g1_guess)/z_max**2
        else:
            g2_guess = 3*np.pi*(HWD*drive/z_max - g1_guess)/z_max/16
            g3_guess = 2/3*(HWD*drive/z_max - g1_guess)/z_max**2
        g2 = ('g2',g2_guess,True,-np.inf,np.inf)
        g3 = ('g3',g3_guess,True,0,np.inf)

        
        if i==0 or HWD == None:
            f_guess = z_max*width
        else:
            f_guess = drive*HWD
            
        f = ('f',f_guess,True,0,np.inf)
        return [a1,a3,g1,g2,g3,f]
    
    def __init__(self,fqs=None,xs=None,ys=None,drives=None,filenames=None):
        ResonanceFit.__init__(self,fqs=fqs,xs=xs,ys=ys,drives=drives,filenames = filenames)
    
class NbTiFit(ResonanceFit):
    def model(self,params,fq, z):
        """
        Implicit equation for speed reasonance curve of nonlinear oscillator, F(fq,z) = 0
        
        Parameters
        ----------
        params : result.params object from lmfit package \n
        fq : frequency \n
        z : amplitude \n
    
        Output
        -----
        Numpy compatible output with regards to fq, z. Implicit equation for reasonance curve, F(fq,z) = 0
        
        Equation specification
        ----------------------
        params : [a1,a2,a3,g1,g2,g3,f]
        
        **Restoring force**
        
        a1 : alpha1 \n
        a3 : alpha2 \n
        
        **Damping**
        
        g1 : gamma1 \n
        g2 : gamma2 \n
        g3 : gamma2 \n
        
        **Forcing**
        
        f : f (forcing amplitude)
        
        **Independent variables**
        
        fq : omega (forcing frequency in rad/s) \n
        z : z (amplitude)
        
        """
        
        a1, a3, g1, g2, g3, f = params.values()
        
        return ((z/fq)**2 *( (a1  + 0.75*a3*(z/fq)**2 - fq**2 )**2
                + fq**2 * (g1 + ( 8/(3*np.pi) )*g2*z + 0.75*g3*z**2)**2 ) - f**2)

    
    def params_estimate(self,par_edge,i:int):
        z = self.zs[i]
        fq = self.fqs[i]
        fq0=par_edge[0]
        fq_max = fq[np.argmax(z)]
        z_max = np.max(z)
        z_half = z_max/2
        if self.drives != None:
            drive=self.drives[i]
        
        #    #name #init #vary #min #max
        a1 = ('a1',fq0**2,True, 0,   np.inf)

        a3_guess = 4/3*(fq_max**2 - fq0**2)*fq_max**2/z_max**2
        a3= ('a3',a3_guess,True, -np.inf,np.inf)
        
        z_left = z[3:np.argmax(z)]
        fq_left = fq[3:np.argmax(z)]
        fq_l = np.interp(z_half, z_left, fq_left)
        
        z_right = z[-3:np.argmax(z):-1]
        fq_right = fq[-3:np.argmax(z):-1]
        fq_r = np.interp(z_half, z_right, fq_right)
        
        width = np.abs(fq_r-fq_l)/np.sqrt(3)
        if i==0:
            g1_guess = width
        else:
            g1_guess = self.init_pars[0]['g1'].value
        g1 = ('g1',g1_guess,True,0,np.inf)
        
        
        if self.drives != None:
            HWD = np.max(self.zs[0])*g1_guess/self.drives[0]
        else:
            HWD = None
        
        if HWD == None:
            g2_guess = 0.5*(width - g1_guess)/z_max
            g3_guess = 0.5*(width - g1_guess)/z_max**2
        else:
            g2_guess = 3*np.pi*(HWD*drive/z_max - g1_guess)/z_max/16
            g3_guess = 2/3*(HWD*drive/z_max - g1_guess)/z_max**2
        g2 = ('g2',g2_guess,True,0,np.inf)
        g3 = ('g3',g3_guess,True,0,np.inf)

        
        if i==0 or HWD == None:
            f_guess = z_max*width
        else:
            f_guess = drive*HWD
            
        f = ('f',f_guess,True,0,np.inf)
        return [a1,a3,g1,g2,g3,f]
    
    def __init__(self,fqs=None,xs=None,ys=None,drives=None,filenames=None):
        ResonanceFit.__init__(self,fqs=fqs,xs=xs,ys=ys,drives=drives,filenames = filenames)
    
    
"""
Not finished - better g3 estimate needed
"""    
class GoalpostFit(ResonanceFit):
     def model(self,params,fq, z):
         """
         Implicit equation for speed reasonance curve of nonlinear oscillator, F(fq,z) = 0
         
         Parameters
         ----------
         params : result.params object from lmfit package \n
         fq : frequency \n
         z : amplitude \n
     
         Output
         -----
         Numpy compatible output with regards to fq, z. Implicit equation for reasonance curve, F(fq,z) = 0
         
         Equation specification
         ----------------------
         params : [a1,a2,a3,g1,g2,g3,f]
         
         **Restoring force**
         
         a1 : alpha1 \n
         a3 : alpha2 \n
         
         **Damping**
         
         g1 : gamma1 \n
         g2 : gamma2 \n
         g3 : gamma2 \n
         
         **Forcing**
         
         f : f (forcing amplitude)
         
         **Independent variables**
         
         fq : omega (forcing frequency in rad/s) \n
         z : z (amplitude)
         
         """
         
         a1, a3, g1, g2, g3, f = params.values()
         
         return ((z/fq)**2 *( (a1  + 0.75*a3*(z/fq)**2 - fq**2 )**2
                 + fq**2 * (g1 + ( 8/(3*np.pi) )*g2*z + 0.75*g3*z**2)**2 ) - f**2)

     
     def params_estimate(self,par_edge,i:int):
         z = self.zs[i]
         fq = self.fqs[i]
         fq0=par_edge[0]
         fq_max = fq[np.argmax(z)]
         z_max = np.max(z)
         z_half = z_max/2
         if self.drives != None:
             drive=self.drives[i]
         #    #name #init #vary #min #max
         a1 = ('a1',fq0**2,True, 0,   np.inf)

         a3_guess = 4/3*(fq_max**2 - fq0**2)*fq_max**2/z_max**2
         a3= ('a3',a3_guess,True, -np.inf,np.inf)
         
         z_left = z[3:np.argmax(z)]
         fq_left = fq[3:np.argmax(z)]
         fq_l = np.interp(z_half, z_left, fq_left)
         
         z_right = z[-3:np.argmax(z):-1]
         fq_right = fq[-3:np.argmax(z):-1]
         fq_r = np.interp(z_half, z_right, fq_right)
         
         width = np.abs(fq_r-fq_l)/np.sqrt(3)
         if i==0:
             g1_guess = width
         else:
             g1_guess = self.init_pars[0]['g1'].value
         g1 = ('g1',g1_guess,True,0,np.inf)
         
         
         if self.drives != None:
             HWD = np.max(self.zs[0])*g1_guess/self.drives[0]
         else:
             HWD = None
         
         if HWD == None:
             g3_guess = 0.5*(width - g1_guess)/z_max**2
         else:
             g3_guess = 2/3*(HWD*drive/z_max - g1_guess)/z_max**2
         g2 = ('g2',0,False,-np.inf,np.inf)
         g3 = ('g3',g3_guess,True,0,np.inf)

         
         if i==0 or HWD == None:
             f_guess = z_max*width
         else:
             f_guess = drive*HWD
             
         f = ('f',f_guess,True,0,np.inf)
         return [a1,a3,g1,g2,g3,f]
     
     def __init__(self,fqs=None,xs=None,ys=None,drives=None,filenames=None):
         ResonanceFit.__init__(self,fqs=fqs,xs=xs,ys=ys,drives=drives,filenames = filenames)
       
        
        


