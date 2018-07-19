#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def saveparti(tt,particle,parti):
    """
    This routine saves the particle data either in a CSV file, or in an SQLite database.
    If saved in a CSV file, one file is created for every recorded timestep.
    If saved as sqlite, data is added to the database at every recorded timestep.
    
    """
    import os
    from sqlalchemy import create_engine
    
    parti2 = parti.copy()
    parti2.columns = ['DOY', 'ID','x','y','z','u','v','w','wsink','wtotal','salinity','temperature','density','PV','vorticity']
    
    ## CSV file format
    if particle.outputformat=='csv':
        print('SAVE PARTICLE DATA in CSV-file... ')
        
        # Define the variable to print in output file
        # exclude velocity field at previous timestep
        #vartoprint = [1,2,6:16]
        
        # Define filename
        FILENAME = os.path.join(particle.outputdir,particle.outputfilename+'_doy'+str(tt)+'.csv')
        
        # test if file already exist
        if os.path.isfile(FILENAME):
            print('File already exist - This code is set to never overwrite a file')
            exit
        
        parti2.to_csv(FILENAME, index=False)
            
    
        ## SQLite database
    elif particle.outputformat=='sqlite':
        print('SAVE PARTICLE DATA in SQLite... ')
        
        # Set database filename
        DBFILENAME = os.path.join(particle.outputdir,particle.outputfilename+'_'+particle.direction+'.db')
        # Set database table name
        tablename = 'particles'
        
    
        # Creates the SQLite database if does not exist
        engine =  create_engine('sqlite:///'+DBFILENAME, echo=False)
        # Creates the table in the database
        
        parti2.to_sql(tablename, con=engine,if_exists='append',index=False)
    
        # Close connection
        engine.dispose()
        print('Disonnected to SQLite database')
        
    else:
        print('Cannot recognize particle output format.')
        exit
    
    print('DONE')
    print('#################')
    print(' ')
    return
