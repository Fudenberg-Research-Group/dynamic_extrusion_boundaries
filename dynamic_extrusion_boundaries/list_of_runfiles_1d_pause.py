import glob
import ast 
import datetime

LEF_pause= 0.9

pause_multiplier = 1/(1-LEF_pause)

lifetimes_b = [16.6, 50, 150]
lifetimes = [pause_multiplier*life for life in lifetimes_b]

velocities = [0.33, 1, 3]

face_stalls = [1.0]
back_stalls = [0]

CTCF_lifetimes_b =  [0.5, 0.75, 1, 1.25, 1.5, 3, 10, 15, 30, 50, 150]
CTCF_lifetimes = [pause_multiplier*clife for clife in CTCF_lifetimes_b]
CTCF_offtimes_b = [0.05, 0.075, 0.1, 0.125, 0.15, 0.3, 1, 1.5, 3, 5, 15]
CTCF_offtimes = [pause_multiplier*cof for cof in CTCF_offtimes_b]

stall_dists = [100]
LEF_birth = 0.1

LEF_separations = [100]
sites_per_monomer = 10
replication_number = 10
monomer_per_replica = 1000

steps_b = 200
steps = steps_b/pause_multiplier

already_processed = []
for fname  in glob.glob('/home1/start/polychrom/projects/Dynamic_boundary_elements/simulations/sims/*_STALL*'):
    already_processed.append(fname.split('/')[-1])

with open(str(datetime.date.today())+'_runfile.txt','w') as f:
    for lifetime in lifetimes:
        for velocity_multiplier in velocities:
            for face_stall in face_stalls:
                for back_stall in back_stalls:
                    for LEF_separation in LEF_separations:
                        for CTCF_lifetime in CTCF_lifetimes:
                            for CTCF_offtime in CTCF_offtimes:
                                paramset = (
                                    'folder_face_'+str(face_stall)+
                                    '_back_'+str(back_stall)+
                                    '_Clife_'+str(CTCF_lifetime) +
                                    '_Cof_'+str(CTCF_offtime)+
                                    '_life_'+str(lifetime)+
                                    '_slife_'+str(lifetime)+
                                    '_birth_'+str(LEF_birth)+
                                    '_pause_'+str(LEF_pause)+
                                    '_sep_'+str(LEF_separation)+
                                    '_site_'+str(sites_per_monomer)+
                                    '_monomer_'+str(monomer_per_replica)+
                                    '_replica_'+str(replication_number)+
                                    '_steps_'+str(steps)+
                                    '_vel_'+str(velocity_multiplier)
                                    )
                                if paramset not in already_processed:
                                    f.write( paramset +'\n')
                                else:
                                    print('already done')
                                
            

paramdict_keys={
                'CTCF_facestall':'face',
                'CTCF_backstall':'back',
                'CTCF_lifetime':'Clife',
                'CTCF_offtime':'Cof',
                'LEF_lifetime':'life',
                'LEF_stalled_lifetime':'slife',
                'LEF_birth':'birth',
                'LEF_pause':'pause',
                'LEF_separation':'sep',
                'sites_per_monomer':'site',
                'monomers_per_replica':'monomer',
                'number_of_replica':'replica',
                'steps':'steps',
                'velocity_multiplier':'vel'
                }
