from Gaudi.Configuration import *
from PhysSelPython.Wrappers import Selection, SelectionSequence, DataOnDemand, \
    AutomaticData
from Configurables import DecayTreeTuple, FitDecayTrees, TupleToolRecoStats, \
    TupleToolTrigger, TupleToolTISTOS, CondDB, SelDSTWriter, FilterDesktop
from DecayTreeTuple.Configuration import *
from Configurables import TupleToolApplyIsolation

from Configurables import EventNodeKiller
eventNodeKiller = EventNodeKiller('Stripkiller')
eventNodeKiller.Nodes = [ '/Event/AllStreams', '/Event/Strip' ]

# Standard stripping20 
name = 'inclmumu'

"""
Options for building Stripping20,
with tight track chi2 cut (<3)
"""

from Configurables import EventNodeKiller
eventNodeKiller = EventNodeKiller('Stripkiller')
eventNodeKiller.Nodes = ['/Event/AllStreams','/Event/Strip']

from Gaudi.Configuration import *
MessageSvc().Format = "% F%30W%S%7W%R%T %0W%M"

# Tighten Trk Chi2 to <3
from CommonParticles.Utils import DefaultTrackingCuts
DefaultTrackingCuts().Cuts  = { "Chi2Cut" : [ 0, 4 ],
                                "CloneDistCut" : [5000, 9e+99 ] }

#
# Build the streams and stripping object
#
from StrippingArchive.Stripping21.StrippingB2XMuMuInclusive \
    import B2XMuMuInclusiveConf as builder
from StrippingArchive.Stripping21.StrippingB2XMuMuInclusive \
    import defaultConfig as config
from StrippingConf.Configuration import StrippingConf, StrippingStream

stripping='stripping20'
lb = builder('B2XMuMuIncl',config)
print config
AllStreams = StrippingStream("MyStream")

for line in lb.lines():
    if line.name() == 'StrippingB2XMuMuIncl_InclDiMuLowQ2Line':
        AllStreams.appendLines([line])


sc = StrippingConf( Streams = [ AllStreams ],
                    MaxCandidates = 2000
                    )

stripsel = AutomaticData(Location = "Phys/B2XMuMuIncl_InclDiMuLowQ2Line/Particles")

stripfilter = FilterDesktop("stripfilter",
                             Preambulo = ["from LoKiPhysMC.decorators import *",
                                          "from LoKiPhysMC.functions import mcMatch"],
                             Code = "ALL")

inclmumu = Selection ("Sel"+name,
                     Algorithm = stripfilter,
                     RequiredSelections = [stripsel])
seq = SelectionSequence("seq",
                      TopSelection = inclmumu)


tuple = DecayTreeTuple("Incl_Tuple")

tuple.Inputs = [stripsel.outputLocation()]
tuple.ToolList =  [
      "TupleToolKinematic"
    , "TupleToolEventInfo"
    , "TupleToolRecoStats"
    , "TupleToolMCBackgroundInfo"
    , "TupleToolGeometry"
    , "TupleToolPid"
    , "TupleToolPropertime"
    , "TupleToolPrimaries"
    , "TupleToolTrackInfo"
]
MCtruth = tuple.addTupleTool("TupleToolMCTruth")
MCtruth.addTupleTool("LoKi::Hybrid::MCTupleTool/LoKi_Photos").Variables = {
    "nPhotons"  : "MCNINTREE(('gamma' == MCABSID))",
    "photons_TRUEP_X" : "MCSUMTREE(MCPX,('gamma' == MCABSID))", 
    "photons_TRUEP_Y" : "MCSUMTREE(MCPY,('gamma' == MCABSID))",
    "photons_TRUEP_Z" : "MCSUMTREE(MCPZ,('gamma' == MCABSID))", 
    "photons_TRUEP_E"  : "MCSUMTREE(MCE,('gamma'==MCABSID))" 
    }
      
tuple.addBranches ({         
      "muplus" :  "[B0 -> ^mu+ mu-]CC",
      "muminus" :  "[B0 -> mu+ ^mu-]CC",
      "B0" : "[B0 -> mu+ mu-]CC",
})

LoKi_All=tuple.addTupleTool("LoKi::Hybrid::TupleTool/LoKi_All")
LoKi_All.Variables = {
        'MINIPCHI2' : "MIPCHI2DV(PRIMARY)", 
        'MINIP' : "MIPDV(PRIMARY)",
        'IPCHI2_OWNPV' : "BPVIPCHI2()", 
        'IP_OWNPV' : "BPVIP()"
}

LoKi_muplus=tuple.muplus.addTupleTool("LoKi::Hybrid::TupleTool/LoKi_muplus")
LoKi_muplus.Variables = {
       'PIDmu' : "PIDmu",
       'ghost' : "TRGHP",
       'TRACK_CHI2' : "TRCHI2DOF",
       'NNK' : "PPINFO(PROBNNK)",
       'NNpi' : "PPINFO(PROBNNpi)",
       'NNmu' : "PPINFO(PROBNNmu)"
}

LoKi_muminus=tuple.muminus.addTupleTool("LoKi::Hybrid::TupleTool/LoKi_muminus")
LoKi_muminus.Variables = {
       'PIDmu' : "PIDmu",
       'ghost' : "TRGHP",
       'TRACK_CHI2' : "TRCHI2DOF",
       'NNK' : "PPINFO(PROBNNK)",
       'NNpi' : "PPINFO(PROBNNpi)",
       'NNmu' : "PPINFO(PROBNNmu)"
}

list = [
      "L0DiMuonDecision"
    , "L0MuonDecision"
    , "Hlt1TrackAllL0Decision"
    , "Hlt1TrackMuonDecision"
    , "Hlt1DiMuonLowMassDecision"
    , "Hlt1DiMuonHighMassDecision"
    , "Hlt1SingleMuonHighPTDecision"
    , "Hlt2TopoMu2BodyBBDTDecision"
    , "Hlt2TopoMu3BodyBBDTDecision"
    , "Hlt2Topo2BodyBBDTDecision"
    , "Hlt2Topo3BodyBBDTDecision"
    , "Hlt2DiMuonDetachedDecision"
    , "Hlt2SingleMuonDecision"
    , "Hlt2DiMuonDetachedHeavyDecision"
]


# Isolation tools
tuple.B0.ToolList += [ "TupleToolApplyIsolation" ]
tuple.B0.addTool(TupleToolApplyIsolation, name="TupleToolApplyIsolation")

# TISTOS tools
tuple.B0.ToolList += [ "TupleToolTISTOS" ]
tuple.B0.addTool( TupleToolTISTOS, name = "TupleToolTISTOS" )
tuple.B0.TupleToolTISTOS.Verbose = True
tuple.B0.TupleToolTISTOS.TriggerList = list
tuple.muplus.ToolList += [ "TupleToolTISTOS" ]
tuple.muplus.addTool( TupleToolTISTOS, name = "TupleToolTISTOS" )
tuple.muplus.TupleToolTISTOS.Verbose = True
tuple.muplus.TupleToolTISTOS.TriggerList = list
tuple.muminus.ToolList += [ "TupleToolTISTOS" ]
tuple.muminus.addTool( TupleToolTISTOS, name = "TupleToolTISTOS" )
tuple.muminus.TupleToolTISTOS.Verbose = True
tuple.muminus.TupleToolTISTOS.TriggerList = list
tuple.Decay = "[B0 -> ^mu+ ^mu-]CC"

from Configurables import DaVinci
DaVinci().TupleFile = "TupleToolBGMC.root"

DaVinci().EvtMax = -1
DaVinci().DataType = '2012'
DaVinci().Simulation   = True
DaVinci().Lumi = False
#CondDB().UseOracle = True
#DaVinci().DDDBtag  = "dddb-20120831"
#DaVinci().CondDBtag = "sim-20121025-vc-md100"
myseq = GaudiSequencer("myseq")
myseq.Members += [ eventNodeKiller ]
myseq.Members += [ sc.sequence() ] 
myseq.Members += [ seq.sequence() ] 
myseq.Members += [tuple]
##DaVinci().UserAlgorithms = [smear,scaler,_myseq]
DaVinci().UserAlgorithms = [myseq]
DaVinci().MainOptions  = ""
