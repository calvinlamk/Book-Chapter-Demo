from __future__ import division
from PySteppables import *       
import math                      
import numpy                   
import CompuCell
import sys
import random
RNG=random.SystemRandom()

#Global Parameters
#Motility Parameters
CtoM=52  
BASAL=100 
SCF=0.5

#EndTime Parameters
ENDMCS=20000 

#Cell Size and Division Parameters
RADAVG=3 
RADDEV=.5 
MTFORCEMIN=-3*10**(-3.88)
MTFORCEMAX=4*10**(-3.88)

#Signaling Parameters

#Constitutive Ligand Parameters
CONEXPSCF=10000 
THETA=0 
XI=1000 
#YG Signaling Parameters
ALPHAYG=1 
BETAYG=1750 
EPSILONYG=1000
KAPPAYG=25000
THRESHOLDUPYG=5263
THRESHOLDDOYG=5263
#BR Signaling Parameters
ALPHABR=1
BETABR=921.181
EPSILONBR=526.389
KAPPABR=25000
THRESHOLDUPBR=5301
THRESHOLDDOBR=5301

#Sampling and Computational Threads Parameters
RESOL=100
USEDNODES=8


class ELUGMSteppable(SteppableBasePy):

    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):
#Graph the Quantification Code
        self.pW1 = self.addNewPlotWindow(
            _title='Types',
            _xAxisTitle='MonteCarlo Step (MCS)',
            _yAxisTitle='Count',
            _xScaleType='linear',
            _yScaleType='linear',
            _grid=True)                
        self.pW1.addPlot('Y', _style='Dots', _color='gray', _size=3)
        self.pW1.addPlot('G', _style='Dots', _color='green', _size=3)
        self.pW1.addPlot('B', _style='Dots', _color='blue', _size=3)   
        self.pW1.addPlot('R', _style='Dots', _color='red', _size=3)
        
        self.pW2 = self.addNewPlotWindow(
            _title='Point System',
            _xAxisTitle='MonteCarlo Step (MCS)',
            _yAxisTitle='Count',
            _xScaleType='linear',
            _yScaleType='linear',
            _grid=True)                
        self.pW2.addPlot('YG', _style='Dots', _color='green', _size=3)
        self.pW2.addPlot('BR', _style='Dots', _color='blue', _size=3)

#Calling Adhesion Values From the XML File Code             
        global YtoY,YtoG,GtoY,YtoB,BtoY,YtoR,RtoY,GtoG,GtoB,BtoG,GtoR,RtoG,BtoB,BtoR,RtoB,RtoR
        YtoY=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','Y','Type2','Y']))
        YtoG=GtoY=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','Y','Type2','G']))
        YtoB=BtoY=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','Y','Type2','B']))
        YtoR=RtoY=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','Y','Type2','R']))
        GtoG=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','G','Type2','G']))
        GtoB=BtoG=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','G','Type2','B']))
        GtoR=RtoG=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','G','Type2','R']))
        BtoB=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','B','Type2','B']))
        BtoR=RtoB=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','B','Type2','R']))
        RtoR=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','R','Type2','R']))
        
#Initialization of Cells Parameters
        for cell in self.cellList:
            cell.dict["RDM"]=RNG.gauss(RADAVG,RADDEV)
            cell.lambdaSurface=2.5                    
            cell.targetSurface=4*math.pi*cell.dict["RDM"]**2  
            cell.lambdaVolume=2.5                     
            cell.targetVolume=(4/3)*math.pi*cell.dict["RDM"]**3 
            cell.dict["PTS"]=[0]

    def step(self,mcs):                 
#Defining Per Step Variables Required for Quantifications Code
        NUMTY=0
        NUMTG=0
        NUMTB=0
        NUMTR=0
        
        YGPTS=0
        BRPTS=0     
        
        if mcs==1:
            self.changeNumberOfWorkNodes(USEDNODES)             

        for cell in self.cellList:
#Defining Per Step Per Cell Variables for Signaling and Motility Code
            CSAY=0
            CSAG=0
            CSAB=0
            CSAR=0
            CSAM=0
            
            PTSY=0
            PTSG=0
            PTSB=0
            PTSR=0
            DTRES=0
#Iterating Over Neighbors Code
            for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell):
                if neighbor is None:
                    continue
                if neighbor.type==1:
                    CSAY+=commonSurfaceArea
                    PTSY+=0
                if neighbor.type==2:
                    CSAG+=commonSurfaceArea
                    PTSG+=commonSurfaceArea*neighbor.dict["PTS"][0]/(neighbor.surface)
                if neighbor.type==3:
                    CSAB+=commonSurfaceArea                 
                    PTSB+=commonSurfaceArea*(CONEXPSCF/(1+math.exp(-(mcs-THETA)/XI)))/neighbor.surface
                if neighbor.type==4:
                    CSAR+=commonSurfaceArea
                    PTSR+=commonSurfaceArea*(CONEXPSCF/(1+math.exp(-(mcs-THETA)/XI)))/neighbor.surface

            CSAM=cell.surface-(CSAY+CSAG+CSAB+CSAR)               
#Changing Reporter as a Result of Signaling Changes           
            if (cell.type==1 or cell.type==2):
                DTRES=(1/(ALPHAYG+math.exp(-((PTSR+PTSB)-BETAYG)/EPSILONYG)))-(1/KAPPAYG)*cell.dict["PTS"][0]
                cell.dict["PTS"][0]+=DTRES
            if (cell.type==3 or cell.type==4):
                DTRES=(1/(ALPHABR+math.exp(-((PTSG)-BETABR)/EPSILONBR)))-(1/KAPPABR)*cell.dict["PTS"][0]
                cell.dict["PTS"][0]+=DTRES
#Changing State as a Result of Reporter Changes               
            if cell.type==1:
                if cell.dict["PTS"][0]>=THRESHOLDUPYG:
                    cell.type=2
            if cell.type==2:
                if cell.dict["PTS"][0]<THRESHOLDDOYG:
                    cell.type=1
            if cell.type==3:
                if cell.dict["PTS"][0]>=THRESHOLDUPBR:
                    cell.type=4
            if cell.type==4:
                if cell.dict["PTS"][0]<THRESHOLDDOBR:
                    cell.type=3
#Defining and Calculating Cell Physical Properties Code/Parameters
            if cell.type==1:                             
                cell.lambdaSurface=2.2          
                cell.lambdaVolume=2.2           
                NUMTY+=1                     
                cell.fluctAmpl=BASAL+SCF*(CtoM*CSAM+YtoY*CSAY+YtoG*CSAG+YtoB*CSAB+YtoR*CSAR)/cell.surface
                YGPTS+=cell.dict["PTS"][0]

            if cell.type==2:
                cell.lambdaSurface=1.0            
                cell.lambdaVolume=1.0             
                NUMTG+=1                       
                cell.fluctAmpl=BASAL+SCF*(CtoM*CSAM+GtoY*CSAY+GtoG*CSAG+GtoB*CSAB+GtoR*CSAR)/cell.surface
                YGPTS+=cell.dict["PTS"][0]
             
            if cell.type==3:                           
                cell.lambdaSurface=2.2          
                cell.lambdaVolume=2.2               
                NUMTB+=1                     
                cell.fluctAmpl=BASAL+SCF*(CtoM*CSAM+BtoY*CSAY+BtoG*CSAG+BtoB*CSAB+BtoR*CSAR)/cell.surface
                BRPTS+=cell.dict["PTS"][0]
                
            if cell.type==4:                    
                cell.lambdaSurface=1.0          
                cell.lambdaVolume=1.0                
                NUMTR+=1                       
                cell.fluctAmpl=BASAL+SCF*(CtoM*CSAM+RtoY*CSAY+RtoG*CSAG+RtoB*CSAB+RtoR*CSAR)/cell.surface
                BRPTS+=cell.dict["PTS"][0]
#Adding Data Points to Graphs Code
        if mcs%RESOL==0:
            
            self.pW1.addDataPoint("Y", mcs, NUMTY)
            self.pW1.addDataPoint("G", mcs, NUMTG) 
            self.pW1.addDataPoint("B", mcs, NUMTB)   
            self.pW1.addDataPoint("R", mcs, NUMTR)

            self.pW2.addDataPoint("YG", mcs, YGPTS/(NUMTY+NUMTG))
            self.pW2.addDataPoint("BR", mcs, BRPTS/(NUMTB+NUMTR))
#Exporting the Graphs as Data in Text Files Code         
            if mcs==ENDMCS:
                fileName = "FOU" + str(mcs) + ".txt"
                self.pW1.savePlotAsData(fileName)
                fileName = "SIG" + str(mcs) + ".txt"
                self.pW2.savePlotAsData(fileName)                
                self.stopSimulation()                 
    def finish(self):
        pass
#Mitosis Codes
from PySteppablesExamples import MitosisSteppableBase

class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,_simulator,_frequency=1):
        MitosisSteppableBase.__init__(self,_simulator, _frequency)
        self.setParentChildPositionFlag(0)
    def step(self,mcs):        
        cells_to_divide=[]         
        for cell in self.cellList:
            cell.dict["RDM"]+=RNG.uniform(MTFORCEMIN,MTFORCEMAX)
            cell.targetSurface=4*math.pi*cell.dict["RDM"]**2
            cell.targetVolume=(4/3)*math.pi*cell.dict["RDM"]**3 
            if cell.volume>2*(4/3)*math.pi*RADAVG**3:            
                cells_to_divide.append(cell)
                
        for cell in cells_to_divide:
            self.divideCellRandomOrientation(cell) 

    def updateAttributes(self):
        self.parentCell.dict["RDM"]=RNG.gauss(RADAVG,RADDEV)
        self.parentCell.targetVolume=(4/3)*math.pi*self.parentCell.dict["RDM"]**3 
        self.parentCell.targetSurface=4*math.pi*self.parentCell.dict["RDM"]**2 
        self.cloneParent2Child()
        