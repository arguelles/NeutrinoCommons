def CalculateParameterScan(datapath):
    """Given a set of data files for several configurations will calculate the ratio between standard and
    the given configuration. Then will save the data.
    
    @type   datapath:   string
    @param  datapath:   path to zip files.
    
    @rtype          :   plot
    @return         :   generates the countour plot
    """
    
    sparam = PC.PhysicsConstants()
    param  = PC.PhysicsConstants()
    
    channels = ['bb','WW']
    
    DM_mass = 1000.0*sparam.GeV
    DMsig_soft   = 1.8e-41*param.cm**2
    DMsig_hard   = 7.2e-43*param.cm**2
    
    ice_mu_flux_lim_hard = 3.6e2
    ice_mu_flux_lim_soft = 1.3e3
    
    ratio_flux_soft = []
    ratio_flux_hard = []
    mu_flux = []    
    mu_flux_std = []
    
    # calculate STD flux values
    print "Calculate STD flux"
    for ch in channels:
        if ch == 'bb':
            sig = DMsig_soft
        elif ch == 'WW':
            sig = DMsig_hard
            
        mu_inter = DM.DMNeuFluxDetNoInt(ch,DM_mass,sig,sparam,onlyosc = False,datapath = datapath)
        # integrating
        int_flux = integrate.quad(lambda E: mu_inter(E)*sparam.km**2*sparam.year,1.0*param.GeV,DM_mass,epsrel=1.0e-10,limit = 250)[0]
        # saving
        mu_flux_std.append(int_flux)
        #print int_flux*sparam.km**2*sparam.year
        
    # calculating non std fluxes
    
    # reading filenames
    filenames = []
    for filename in os.listdir(datapath):
        if fnmatch.fnmatch(filename,"output*.zip"):
            filenames.append(filename)
            
    for i,filename in enumerate(filenames):
        
        # open zip file #
        zipfile = zp.ZipFile(datapath+filename, mode = "r")
        zipfile.extractall(path = datapath)
        
        #begin get param data
        num = re.search('(?<=_)\w+',filename).group(0)
        tmp_filename = "DataNeuOscProb_RK_neutrino_Emin_1.0_GeV_Emax_1000000.0_GeV_ineu_0_param_CONF"+num+".dat"
        data = []
        file = open(datapath+tmp_filename,'r')
        gt.hreadfilev4(file,data,param)
        file.close()
        
        #print param.name,param.th14,param.dm41sq
        th = param.th14
        dmsq = param.dm41sq
        #end get param data
        print "Calculating for CONF"+str(num)
        
        flux = []
        
        for ch in channels:
            if ch == 'bb':
                sig = DMsig_soft
            elif ch == 'WW':
                sig = DMsig_hard
                
            mu_inter = DM.DMNeuFluxDetNoInt(ch,DM_mass,sig,param,onlyosc = False,datapath = datapath)
            # integrating
            #int_flux = integrate.quad(mu_inter,1.0*param.GeV,DM_mass)[0]
            # testing integration precision
            #int_flux = integrate.quad(mu_inter,1.0*param.GeV,DM_mass,epsrel=1.0e-20,limit=100)
            #
            #print th,dmsq
            #print int_flux[0],int_flux[1]
            #print int_flux[0]*sparam.km**2*sparam.year
            
            iflux = integrate.quad(lambda E : mu_inter(E)*sparam.km**2*sparam.year,1.0*param.GeV,DM_mass,epsrel=1.0e-10,limit = 250)
            int_flux  = iflux[0]
            #int_error = iflux[1]
            ##int_flux = integrate.quadrature(lambda E : mu_inter(E)*sparam.km**2*sparam.year,1.0*param.GeV,DM_mass,tol=1.0e-8,maxiter = 500, vec_func = False)[0]
            #
            #print int_flux,int_error
            #
            #if ch == 'bb':
            #    ratio = mu_flux_std[0]/int_flux
            #elif ch == 'WW':
            #    ratio = mu_flux_std[1]/int_flux
            #    
            #print ratio
            #
            #if ch == 'bb':
            #    iratio = integrate.quad(lambda E : mu_inter(E)*sparam.km**2*sparam.year/mu_flux_std[0],1.0*param.GeV,DM_mass,epsrel=1.0e-10,limit = 250)[0]
            #elif ch == 'WW':
            #    iratio = integrate.quad(lambda E : mu_inter(E)*sparam.km**2*sparam.year/mu_flux_std[1],1.0*param.GeV,DM_mass,epsrel=1.0e-10,limit = 250)[0]
            #
            #print ratio-iratio**-1
            
            
            # saving
            #flux.append(int_flux*sparam.km**2*sparam.year)
            flux.append(int_flux)
            #print int_flux*sparam.km**2*sparam.year
            
        mu_flux.append([th,dmsq,flux])
        # clean up
        print "cleaning up"
        for filename in os.listdir(datapath):
            if fnmatch.fnmatch(filename,"DataNeuOscProb_RK_*_Emin_*_GeV_Emax_*_GeV_ineu_*_param_CONF"+num+".dat"):
                os.remove(datapath+filename)
        
        #if i > 10:
        #    quit()
        #    break
    
    for flux_conf in mu_flux:
        th,dmsq,flux = flux_conf
        
#        R_soft = flux[0]/mu_flux_std[0]
#        R_hard = flux[1]/mu_flux_std[0]
        
        R_soft = mu_flux_std[0]/flux[0]
        R_hard = mu_flux_std[1]/flux[1]        
        
        #ratio_flux_hard.append([th,dmsq,R_hard**-1])
        #ratio_flux_soft.append([th,dmsq,R_soft**-1])
        
        ratio_flux_hard.append([th,dmsq,R_hard])
        ratio_flux_soft.append([th,dmsq,R_soft])        
    
    #begin save this data #
    filename_soft = "3+3_soft_ratio.dat"
    file = open(datapath+filename_soft,'w')
    gt.quickprint(file,ratio_flux_soft)
    file.close()
    filename_hard = "3+3_hard_ratio.dat"
    file = open(datapath+filename_hard,'w')
    gt.quickprint(file,ratio_flux_hard)
    file.close()
    #end save this data #


def GenerateConvolutedDMFlux(param):
    DMm     = [250.0*param.GeV,500.0*param.GeV,1000.0*param.GeV,3000.0*param.GeV,5000.0*param.GeV]
    DMsig_soft   = [0.0,2.5e-41*param.cm**2,1.8e-41*param.cm**2,5.3e-41*param.cm**2,1.1e-40*param.cm**2]
    DMsig_hard   = [3.70e-43*param.cm**2,2.9e-43*param.cm**2,7.2e-43*param.cm**2,7.4e-42*param.cm**2,2.6e-41*param.cm**2]
    
    channels = ['WW','bb']
    
    for ch in channels:
        
        if ch == 'bb':
            DMsig = DMsig_soft
        elif ch == 'WW':
            DMsig = DMsig_hard
        else :
            print "Missing cross sections for this channel."
            quit()
            
        dat_pair = [[DMm[i],DMsig[i]] for i in range(len(DMm))]        
        
        for dm,sig in dat_pair:
            print "BEGIN :",ch,dm,sig
            
            param.neutype = "neutrino"
            datapath = "../data/myMC/"+param.name+"/"+param.neutype+"/"
            DM.DMFNeuFluxMCDetv2(ch,dm,sig,param,use_old_data = False,datapath = datapath,crosscheck = False,binnum=100)
            
            param.neutype = "antineutrino"
            datapath = "../data/myMC/"+param.name+"/"+param.neutype+"/"
            DM.DMFNeuFluxMCDetv2(ch,dm,sig,param,use_old_data = False,datapath = datapath,crosscheck = False,binnum=100)
            
            print "END :",ch,dm/param.GeV,sig/param.cm**2
            raw_input("Press Enter to continue...")

