import FWCore.ParameterSet.Config as cms

def customizeForAlpakaSplitter(process):
    process.trial = cms.EDProducer("trial",  # Plugin name
        configString=cms.string("This is my configuration string"),  # Matches configString_
        nHits=cms.uint32(100),                                       # Matches nHits_
        offset=cms.int32(10),                                        # Matches offset_
        ptMin=cms.double(200),                                       # Matches ptMin_
        deltaR=cms.double(0.05),                                     # Matches deltaR_
        chargeFracMin=cms.double(2.0),                               # Matches chargeFracMin_
        tanLorentzAngle=cms.double(0.02),                            # Matches tanLorentzAngle_
        tanLorentzAngleBarrelLayer1=cms.double(0.015),               # Matches tanLorentzAngleBarrelLayer1_
        expSizeXAtLorentzAngleIncidence=cms.double(0.1),             # Matches expSizeXAtLorentzAngleIncidence_
        expSizeXDeltaPerTanAlpha=cms.double(0.02),                   # Matches expSizeXDeltaPerTanAlpha_
        expSizeYAtNormalIncidence=cms.double(0.1),                   # Matches expSizeYAtNormalIncidence_
        centralMIPCharge=cms.double(26000),                          # Matches centralMIPCharge_
        chargePerUnit=cms.double(2000),                              # Matches chargePerUnit_
        forceXError=cms.double(100),                                 # Matches forceXError_
        forceYError=cms.double(150),                                 # Matches forceYError_
        fractionalWidth=cms.double(0.4),                             # Matches fractionalWidth_
        verbose=cms.bool(False)                                      # Matches verbose_
    )

    process.trial_step = cms.Path(process.trial)

    if hasattr(process, 'schedule'):
        process.schedule.extend([process.trial_step])
    else:
        process.schedule = cms.Schedule(process.trial_step)

    return process

