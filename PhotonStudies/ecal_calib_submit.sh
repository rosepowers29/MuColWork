Universe = Vanilla
+SingularityImage = "osdf:///ospool/uc-shared/project/futurecolliders/rosep8/k4toroid.sif"
Executable = run_ecal_calib.sh
Requirements = ( HAS_SINGULARITY ) && ( HAS_CVMFS_unpacked_cern_ch )
should_transfer_files = YES
Output = output.out.$(Cluster)-$(Process)
Log = log.$(Cluster)
Error = error.out.$(Cluster)-$(Process)
transfer_input_files = osdf:///ospool/uc_shared/project/futurecolliders/rosep8/k4toroid.sif, osdf:///ospool/uc-shared/project/futurecolliders/data/fmeloni/DataMuC_MuColl10_v0A/v2/reco/photonGun_E_0_50/, osdf:///ospool/uc-shared/project/futurecolliders/data/fmeloni/DataMuC_MuColl10_v0A/v2/reco/photonGun_E_50_250/, osdf:///ospool/uc-shared/project/futurecolliders/data/fmeloni/DataMuC_MuColl10_v0A/v2/reco/photonGun_E_250_1000/, /scratch/rosep8/PhotonStudies/theta_calib.py
when_to_transfer_output = ON_EXIT
request_cpus = 1
request_disk = 4 GB
request_memory = 4 GB
+ProjectName = "collab.futurecolliders"
Queue 1
