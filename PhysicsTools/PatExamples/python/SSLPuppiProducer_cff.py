import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatExamples.SSLPuppiProducer_cfi import SSLPuppiProducer as _SSLPuppiProducer

SSLPuppiProducer = _SSLPuppiProducer.clone(
    Client = dict(
        timeout = 300,
        mode = "Async",
        modelName = "PuppiGNN",
        modelConfigPath = "/afs/cern.ch/user/y/yujil/CMSSW_13_3_1/src/TritonDemo/models/PuppiGNN/config.pbtxt",
        # version "1" is the resolutionTune
        # version "2" is the responeTune
        modelVersion = "1",
    ),
    pf_src = "packedPFCandidates",
)
