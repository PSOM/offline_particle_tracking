## This routine prints the particle-tracking experiment configuration to a
# text file. The filename includes the date and time of the simulation.
def saveconfiguration(model,particle):
    print('SAVE EXPERIMENT CONFIGURATION... ')
    import datetime

    # Get current time with datetime
    nowtime = datetime.datetime.now().strftime(' %d%m%y-%H%M')

    # Open the file
    file = open('%s/%s_%s_configuration.out' %(particle.outputdir,nowtime,particle.outputfilename),'w')

    # Write the data
    file.write('******** OCEAN MODEL PARAMETERS ********\n')
    file.write('Start day =  %f \n' %model['start_day'])
    file.write('Ocean model timestep =  %f seconds \n' %model['timestep'])
    file.write('Directory of ocean model outputs = %s \n' %model['path'])
    file.write('Ocean model resolution used = every  %f file(s) \n' %model['outputskip'])
    file.write('Ocean model periodicities: E-W =  %f  N-S =  %f \n' %(model['periodic_ew'],model['periodic_ns']))
    file.write('\n')
    file.write('\n')
    file.write('******** PARTICLE TRACKING PARAMETERS ********\n')
    file.write('\n')
    file.write('#----------------\n')
    file.write('SEEDING\n')
    file.write('#----------------\n')
    file.write('Particle initialization (1st seeding) = Day  %f \n' %particle.initime)
    file.write('Particle seeding frequency = every  %f day(s) \n' %particle.inifreq)
    file.write('Number of seeding events =  %f \n' %particle.ininumber)
    file.write('Type of seeding (static vs dynamic) = %s \n' %particle.initype)
    file.write('\n')
    file.write('Number of particle classes =  %f \n' %particle.numofclasses)
    file.write('\n')
    file.write('Particle seeding strategy in x:\n')
    file.write('istart =  %f \n' %particle.istart)
    file.write('irange =  %f \n' %particle.irange)
    file.write('irez =  %f \n' %particle.irez)
    file.write('Particle seeding strategy in y:\n')
    file.write('jstart =  %f \n' %particle.jstart)
    file.write('jrange =  %f \n' %particle.jrange)
    file.write('jrez =  %f \n' %particle.jrez)
    file.write('Particle seeding strategy in z:\n')
    file.write('kstart =  %f \n' %particle.kstart)
    file.write('krange =  %f \n' %particle.krange)
    file.write('krez =  %f \n' %particle.krez)
    file.write('\n')
    file.write('#----------------\n')
    file.write('ADVECTION\n')
    file.write('#----------------\n')
    file.write('Timestep for particle tracking =  %f day(s)\n' %particle.timestep)
    file.write('Length of particle tracking simulation =  %f day(s)\n' %particle.length)
    file.write('Particle tracking direction: %s \n' %particle.direction)
    file.write('\n')
    file.write('#----------------\n')
    file.write('OUTPUT\n')
    file.write('#----------------\n')
    file.write('Particle output frequency =  %f day(s)\n' %particle.outfreq)
    file.write('Particle output format = %s \n' %particle.outputformat)
    file.write('Particle output directory =%s \n' %particle.outputdir)
    file.write('Particle output filename = %s \n' %particle.outputfilename)

    # Close the file
    file.close()

    print('DONE\n')
    print('#################\n')
    print('\n')
    return