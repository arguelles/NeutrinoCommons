"""
General tools for writing/reading data,simple numerics, creating logscales, etc.
"""

import fnmatch
cimport numpy as np
cimport cython

# external class definitions #
#cdef extern from "Python.h":
#        #ctypedef struct FILE
#        ctypedef file FILE

#ctypedef np.ndarray ARRAY

#################################################################################
################################ GENERAL TOOLS ##################################
#################################################################################

############################ FILE INPUT / OUTPUT ################################

def quickprint(file,array):
        """ Write array into file without a header.
    
        @type  file     : file
        @param file     : file open in write/append mode.
        @type  array    : array
        @param array    : array to save.
    
        @rtype          : void
        @return         : void
        """
        # WRITE ARRAY TO FILE
        file.write(str(array).replace('[','').replace(']','\n').replace(',',''))
	
def quickread(file):
        """ Read array from file without a header.
    
        @type  file     : file
        @param file     : file open in read mode.
        @type  array    : array
        @param array    : array to save.
    
        @rtype          : void
        @return         : void
        """
        # READ DATA FILE TO ARRAY
        data = []
        while 1 :
                line = file.readline()
                if not line : break
                if line.isspace() : break
                data.append([float(s) for s in line.split()])
        return data

def hwritefile(file,array,param,header = True):
        """WRITE DATA INTO FILE USING HEADER INFORMATION.

        USED WITH hreadfilev4

        @type   file     :      file  
        @param  file     :      file open in write/append mode.
        @type   array    :      list or numpy array
        @param  array    :      array to save.
        @type   param    :      header
        @param  param    :      parameters included in header.
	@type	overwrite:	bool
	@param	overwrite:	if true will overwrite the file, else will append data	
        
        @rtype           :      void
        @return          :      void
        """
        # BEGIN WRITE HEADER
        p_dict = []
#        for slot in dir(param):
#                if slot != 'Refresh' and slot != 'MemoryLoadObjects':
#                        p_dict.append([slot,getattr(param,slot)])
        saveprop = ['name','linestyle','markerstyle','colorstyle','savefilename',
                    'numneu','numneumax','neutype','th12','th13','th23',
                    'th14','th24','th34','th15','th25','th35','th45','delta21',
                    'delta32','delta31','dm21sq','dm31sq','dm32sq','dm41sq','dm51sq']
        for slot in saveprop:
                p_dict.append([slot,
                               type(getattr(param,slot)).__name__,
                               getattr(param,slot)
                               ])
        if header : 
                file.write('######################################################\n')
                for dic in p_dict:
                        file.write('# '+dic[0]+' '+dic[1]+' '+str(dic[2])+'\n')
                file.write('######################################################\n')
        # END WRITE HEADER
        # WRITE DATA
        file.write(str(array).replace('[','').replace(']','\n').replace(',',''))
        #file.close()

def hreadfilev2(file):
        """ Reads file header and data omiting lines that start with '#'.
        @type  file  : file
        @param file  : file open in read mode.

        @rtype       : list
        @return      : header, data
        """
        data = []
        header = []
        while 1 :
                line = file.readline()
                if not line : break
                elif line.isspace() : break
                elif fnmatch.fnmatch(line,'#*') :
                        header.append([s for s in line.split()][-1])
                else :
                        data.append([float(s) for s in line.split()])
        return header,data

def hreadfilev3(file):
        """ Reads file header and data omiting lines that start with '#'
        separetes chunks of data by blank lines.

        WARNING : Small bug when two blank lines are together. 
        
        @type  file  : file
        @param file  : file open in read mode.
        
        @rtype       : list
        @return      : header, data
        """
        tdat = []
        data = []
        header = []
        while 1 :
                line = file.readline()
                if not line : break
                elif line.isspace() : 
                        tdat.append(data)
                        data = []
                elif fnmatch.fnmatch(line,'#*') :
                        header.append([s for s in line.split()][-1])
                else :
                        data.append([float(s) for s in line.split()])
        # WARNING SMALL BUG WHEN TWO BLANCS TOGETHER # MUST FIX
        return header,tdat

def hreadfilev4(file,data,param,header_read = True):
        """READ FILE HEADER AND DATA OMITING LINES THAT START WITH #
           SEPARATING CHUNKS OF DATA BY SPACE AND CONSIDER HEADER TYPES

        USE WITH hwritefile

        @type   file     :      file  
        @param  file     :      file open in read mode
        @type   data     :      list or numpy array
        @param  data     :      array to read on.
        @type   param    :      header
        @param  param    :      parameters included in header.

        
        @rtype           :      void
        @return          :      void
        """
        tdata   = []
        header  = []
        while 1 :
                line = file.readline()
                if not line :
                        data.append(tdata)
                        break
                elif line.isspace(): 
                        data.append(tdata)
                        tdata = []
                elif fnmatch.fnmatch(line,'#*') :
                        header.append([s for s in line.split()])
                else :
                        tdata.append([float(s) for s in line.split()])
        ## UPDATING PARAM WITH HEADER INFORMATION
        if header_read:
                for dic in header :
                        if len(dic) == 4:
                                if dic[2] == 'str':
                                        setattr(param,dic[1],str(dic[3]))
                                elif dic[2] == 'float' or dic[2] == 'float64':
                                        setattr(param,dic[1],float(dic[3]))
                                elif dic[2] == 'int':
                                        setattr(param,dic[1],int(dic[3]))
                                else:
                                        print "ERROR:LEGION:SOLARNEU:GENERALTOOLS: Wrong type in header."
                                        break

        try :
                param.Refresh()
        except :   
                pass


def Compareparams(param1,param2):
        """Compare if neutrino oscillations parameters are the same in the
        two parameter set.

        @type   param1    :      header
        @param  param1    :      set of parameters 1
        @type   param2    :      header
        @param  param2    :      set of parameters 1

        
        @rtype           :      bool
        @return          :      True if param1 == param2.
        """
        compprop = ['name','linestyle','markerstyle','colorstyle','savefilename',
                    'numneu','numneumax','neutype','th12','th13','th23',
                    'th14','th24','th34','th15','th25','th35','th45','delta21',
                    'delta32','delta31','dm21sq','dm31sq','dm32sq','dm41sq','dm51sq']

        p_dict1 = []
        p_dict2 = []

        for slot in compprop:
                p_dict1.append([slot,
                               type(getattr(param1,slot)).__name__,
                               getattr(param1,slot)
                               ])
                p_dict2.append([slot,
                               type(getattr(param2,slot)).__name__,
                               getattr(param2,slot)
                               ])

        comp = True
        for i in range(len(p_dict1)):
                if p_dict1[i][1] == 'str':
                        if p_dict1[i][2] != p_dict2[i][2]:
                                comp = False
                elif p_dict1[i][1] == 'float' or p_dict1[i][1] == 'float64':
                        if np.abs(p_dict1[i][2]-p_dict2[i][2])>1.0e-6:
                                comp = False
                elif p_dict1[i][1] == 'int':
                        if p_dict1[i][2] != p_dict2[i][2]:
                                comp = False
        return comp


################################### NUMERICS ##########################################
def norm(cnum):
        """ Calculate norm of a complex number.
        @type  cnum     : complex
        @param cnum     : complex number (z)

        @rtype          : float
        @return         : norm(z)
        """
        return np.sqrt(cnum.real**2+cnum.imag**2)

def eigenvectors(M):
        """ Calculates the eigenvectors and eigenvalues ordered by eigenvalue size
        
        @type  M     : matrix
        @param M     : matrix M
        
        @rtype       : list
        @return      : [eigenvalues list, eigenvector list]
        """
        D,V = np.linalg.eig(M)
        DV = []
        VT = V.T
        for i,eigenvalue in enumerate(D):
                DV.append([eigenvalue,VT[i]]) 
        
        DV = sorted(DV,key = lambda x : x[0].real)#np.abs(x[0].real))
        
        V2 = []
        for e in DV:
                V2.append(e[1])
        return D,V2

def Log10SpaceEnergies(Emin,Emax):
        """Returns log10 space energies.

        @type   Emin     :      float
        @param  Emin     :      min energy
        @type   Emax     :      float
        @param  Emax     :      max energy
        
        @rtype           :      list
        @return          :      log spaced energies
        """
        E_list = []
        ordmag = int(np.log10(Emax)-np.log10(Emin))
        for n in range(1,ordmag+1):
                for EE in np.arange(Emin,Emin*10,Emin):
                        E_list.append(EE)
                Emin = Emin*10
        E_list.append(Emin)
        
        return E_list

def LogSpaceEnergies(Emin,Emax,binnum = 30):
        """Returns log space energies.

        @type   Emin     :      float
        @param  Emin     :      min energy
        @type   Emax     :      float
        @param  Emax     :      max energy
        
        @rtype           :      list
        @return          :      log spaced energies
        """
        E_list = []
        Emin_log, Emax_log = np.log(Emin),np.log(Emax)
        if Emin_log == -float('inf') :
                Emin_log = 0.0
        Estep_log = (Emax_log-Emin_log)/float(binnum)
        Elog = np.arange(Emin_log,Emax_log+Estep_log,Estep_log)
        for logE in Elog:
                E_list.append(np.exp(logE))

        return E_list

def MidPoint(bin_list):
        """Returns the mid point list of a given list of bins
        @type   bin_list :      array
        @param  bin_list :      bin array
        
        @rtype           :      list
        @return          :      mid energies
        """
        mid_point = []
        for i,E in enumerate(bin_list):
                if i != len(bin_list)-1:
                        mid_point.append((bin_list[i]+bin_list[i+1])/2.0)
        return mid_point

if __name__ == '__main__':
        pass
