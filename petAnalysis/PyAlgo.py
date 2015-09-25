
from Centella.AAlgo import AAlgo
from Centella.physical_constants import *

from array import *

from string import upper

class PyAlgo(AAlgo):

    def __init__(self,param=False,level = 1,label="",**kargs):

        """
        """
        
        self.name='PyAlgo'
        AAlgo.__init__(self,param,level,self.name,0,label,kargs)
        
        
        
    def initialize(self):

        """        
        """
    
        self.m.log(1,'+++Init method of PyAlgo algorithm+++')
        self.hman.h1("MyPyHisto","MyPyHisto",100,0,10)

        self.tman.book("petalo","petalo")

        self.eventID = array('i',[0])
        self.x = array('d',[0])
        self.y = array('d',[0])
        self.z = array('d',[0])
#        self.plane0 = array('d',[0])
        self.plane0 = array('d',(0,)*64)

        #test
        type=upper(self.plane0.typecode)
        print type
        print "plane0"+'['+str(64)+']/'+type

        self.tman.addBranch("petalo","id",self.eventID)
        self.tman.addBranch("petalo","plane0",self.plane0,64)
        self.tman.addBranch("petalo","x",self.x)
        self.tman.addBranch("petalo","y",self.y)
        self.tman.addBranch("petalo","z",self.z)

        return

    def execute(self,event=""):

        """
  
        """
        
        cppalgo =  self.logman["CNTJob"].instances["petAnalysis"]
        #sipm = cppalgo.fetch_ivstore("sipm")
        #print sipm[0]

        #var = cppalgo.fetch_istore("eventID")
        #self.m.log(1,'This is EventID from from C++Algo:',var)

        self.eventID[0] = cppalgo.fetch_istore("eventID")
        self.x[0] = cppalgo.fetch_dstore("x")
        self.y[0] = cppalgo.fetch_dstore("y")
        self.z[0] = cppalgo.fetch_dstore("z")
        self.plane0 = cppalgo.fetch_dvstore("plane0")

        for i in range(0,64):
                print "%d, " %self.plane0[i]

        self.tman.fill("petalo");
        
        return True

    def finalize(self):

        self.m.log(1,'+++End method of PyAlgo algorithm+++')

        return

    
