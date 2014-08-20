#!/usr/local/bin/python

# This is adapted from Christian Blouin's
# qscript code and edited in order to have
# a String passed to it from another Python script.

# imports
import os

def returnScript(command, jobname='default_job', np=1, npmin=1, nomail=True, sync=False): 
    # template, default values
    scriptDir = os.path.dirname(os.path.realpath(__file__))
    template = os.path.join(scriptDir,'template.sh')

    try:
        email = open("%s/.qsubemail"%(os.environ["HOME"])).read()
    except:
        email = ''
        nomail = True

    # Fill in the template
    script = open(template).read()
    script = script.replace("#NAME#",jobname)
    script = script.replace("#EMAIL#",email)
    script = script.replace("#NPMIN#", str(npmin))
    script = script.replace("#NP#", str(np))

    if np == 1:
        script = script.replace("#CMD#", command)
    else:
        script = script.replace("#mpirun -np $NSLOTS", "mpirun -np %d %s"%(np,command))
        script = script.replace("#CMD#", "")
            
    if nomail:
        script = script.replace('#$ -M', '##$ -M')
        script = script.replace('#$ -m e', '#$ -m n')

    # Write to file
    fout = open("%s.sh"%(jobname),'w')
    fout.write(script)
    fout.close()

    if sync: return 'qsub -sync y %s.sh'%(jobname)
    return 'qsub %s.sh'%(jobname)
