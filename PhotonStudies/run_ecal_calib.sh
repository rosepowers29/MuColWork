#!/bin/bash
apptainer run --no-home -B /ospool/uc-shared/project/futurecolliders/data/fmeloni/:/data -B /home/$USER -B /scratch/$USER /ospool/uc-shared/project/futurecolliders/rosep8/k4toroid.sif
source /setup.sh
python /scratch/rosep8/PhotonStudies/theta_calib.py -m /ospool/uc-shared/project/futurecolliders/data/fmeloni/DataMuC-MuColl10_v0A/v2/reco/photonGun_E_0_50/ -n /ospool/uc-shared/project/futurecolliders/data/fmeloni/DataMuC-MuColl10_v0A/v2/reco/photonGun_E_50_250/ -o /ospool/uc-shared/project/futurecolliders/data/fmeloni/DataMuC_MuColl10_v0A/v2/reco/photonGun_E_250_1000/

