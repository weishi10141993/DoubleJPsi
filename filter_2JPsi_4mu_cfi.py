import FWCore.ParameterSet.Config as cms

def filter_2JPsi_Pt3_4mu(process):

    # Double J/psi filter
    process.oniagenfilter = cms.EDFilter("MCParticlePairFilter",
                                         Status = cms.untracked.vint32(2, 2),
                                         ParticleCharge = cms.untracked.int32(0),
                                         ParticleID1 = cms.untracked.vint32(443),
                                         ParticleID2 = cms.untracked.vint32(443)
                                         )
    # Four muons filter
    process.fourmugenfilter = cms.EDFilter("MCMultiParticleFilter",
                                           src = cms.untracked.InputTag("generatorSmeared"),
                                           NumRequired = cms.int32(4),
                                           AcceptMore = cms.bool(True),
                                           ParticleID = cms.vint32(13),
                                           Status = cms.vint32(1),
                                           PtMin = cms.vdouble(3.0),
                                           EtaMax = cms.vdouble(2.5),
                                           )

    ## filter sequence
    process.gen_mu_select = cms.Sequence(process.oniagenfilter * process.fourmugenfilter)
    process.generation_step += process.gen_mu_select
    return process
