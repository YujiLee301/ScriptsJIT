import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatExamples.SSLPuppiProducer_cff import SSLPuppiProducer

from Configuration.ProcessModifiers.enableSonicTriton_cff import enableSonicTriton
process = cms.Process('PuppiGNN',enableSonicTriton)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input=cms.untracked.int32(10)
)

# Input source
process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(
                                '/store/mc/RunIISummer20UL18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v11_L1v1-v1/260000/03794341-C401-CC45-B5FC-D11264E449CE.root'
                            ),
                            )

process.load("HeterogeneousCore.SonicTriton.TritonService_cff")
process.TritonService.verbose = False
#process.TritonService.fallback.useDocker = True
process.TritonService.fallback.verbose = False
# uncomment this part if there is one server running at 0.0.0.0 with grpc port 8001
process.TritonService.servers.append(
    cms.PSet(
        name = cms.untracked.string("default"),
        address = cms.untracked.string("0.0.0.0"),
        port = cms.untracked.uint32(8001),
    )
)


process.SSLPuppiProducer = SSLPuppiProducer.clone(
)
process.SSLPuppiProducer.Client.verbose = cms.untracked.bool(True)

process.p = cms.Path()
process.p += process.SSLPuppiProducer

process.output = cms.OutputModule("PoolOutputModule",
                                  outputCommands=cms.untracked.vstring(
                                      'keep *'),
                                  fileName=cms.untracked.string(
                                      "SSLPuppiTestSonic.root"
                                      )
                                  )
process.outpath  = cms.EndPath(process.output)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool( True ),
    numberOfThreads = cms.untracked.uint32( 1 ),
    numberOfStreams = cms.untracked.uint32( 0 ),
    sizeOfStackForThreadsInKB = cms.untracked.uint32( 10*1024 )
)
