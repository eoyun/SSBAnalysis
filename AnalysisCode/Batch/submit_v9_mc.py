import os
import sys
import subprocess


def Make_CondorScr(sampleName,outputName,NumFileList,NumJob) :
    path = os.getcwd()
    os.system("mkdir -p %s/%s"%(outputName,sampleName))
    path=path.replace("/Batch", "/")
    os.getcwd()
    os.system ('pwd')
    os.system ('mkdir -p Log')
    os.system ('ls')
    print path
    numjobs = NumJob
    numfiles = NumFileList
    SampleFile = sampleName
    SubFile = SampleFile + "_"
    #qsubcmd = 'qsub -q short '
    Lists = MakeSeparateList(numfiles,numjobs)
    for i in range(0,numjobs):
        idx_ = i +1
        runing = outputName + "/"+SampleFile + "_%s.sh" % (idx_)
        f = open( runing, 'w')
        f.write('#!/bin/tcsh \n')
        f.write('setenv SCRAM_ARCH slc6_amd64_gcc530 \n')
        f.write('source /cvmfs/cms.cern.ch/cmsset_default.csh \n')
        f.write('setenv X509_USER_PROXY $1\n')
	f.write('voms-proxy-info -all\n')
	f.write('voms-proxy-info -all -file $1\n')
	f.write('setenv LD_PRELOAD "/usr/lib64/libpdcap.so" \n')
        f.write('cd '+ path + " \n")
        #f.write('mkdir -p ./output/%s/%s \n'%(sampleName,outputName))
        f.write('mkdir -p ./output/%s/%s \n'%(outputName,sampleName))
        f.write('cmsenv \n')
        SampList = MakeSampleIdxList(SampleFile,Lists[i])
        print "SampList %s"%(SampList)
        f.write('set inputlists = (%s) \n'%(SampList) )
        f.write('foreach i ( $inputlists )\n')
        f.write('   mkdir -p output \n')
        f.write('   ./ssb_analysis %s/${i}.list %s/%s/${i}.root 0 %s \n'%(sampleName,outputName,sampleName,sampleName))
        f.write('end \n')
        #runcmd = "./ssb_analysis "+ SampleFile + "_%s.list " % (i) + SampleFile + "_%s.root " % (i) + " 0 "
        #f.write(runcmd)
        f.close()
        runchMod = "chmod 755 " + runing
        os.system (runchMod)

        ### making submit File ###
        f1 = open( "./%s/%s/condor_%s_%s.submit"%(outputName,sampleName,SampleFile,idx_) , 'w')
        f1.write('# Unix submit description file\n')
        f1.write('Universe = vanilla\n')
        f1.write('Executable = ./%s/%s_%s.sh\n'%(outputName,SampleFile,idx_))
        f1.write('request_memory = 500 \n')
        f1.write('should_transfer_files   = Yes\n')
        f1.write('Output      = ./%s/%s/%s_%s.output\n'%(outputName,sampleName,SampleFile,idx_))
        f1.write('error       = ./%s/%s/errors_%s_%s.txt\n'%(outputName,sampleName,SampleFile,idx_))
        f1.write('log         = ./%s/%s/test_%s_%s.log\n'%(outputName,sampleName,SampleFile,idx_))
        f1.write('Queue 1\n')
        f1.close()
        subchMod = "condor_submit " + "./%s/%s/condor_%s_%s.submit"%(outputName,sampleName,SampleFile,idx_)
        os.system (subchMod)
    pass

def MakeSampleIdxList(Sample,Lists):
    SampleList =""
    for index_ in Lists:
        print "index_ %s in MakeSampleIdxList "%(index_)
        SampleList += '"%s_%s"'%(Sample,index_)
        SampleList += " "
        pass
        #SamleList
    print SampleList
    return SampleList
    pass

def MakeSeparateList(NumFiles, NumJob):
    quo = NumFiles/NumJob
    seplists=[]
    for inx_ in range(NumJob) :
        emptyarray =[]
        seplists.append(emptyarray)
    print "size of seplists : %s "%len(seplists)
    for inumfile in  range (1,NumFiles+1):
        #print "%s"%(inumfile%NumJob)
        seplists[inumfile%NumJob].append(inumfile)
    for idx_ in seplists:
        print "content of %s "%(idx_)
    return seplists
    pass


#Make_CondorScr("SingleNeutrino_v1",40,40)
if __name__ == '__main__':
    Study = "Test_v20"
    #Make_CondorScr("TTJets_Signal",Study,20,10)
    Make_CondorScr("TTJets_others",Study,621,100)
    Make_CondorScr("DYJetsToLL_M_50",Study,732,100)
    Make_CondorScr("ST_tW_antitop",Study,70,70)
    Make_CondorScr("ST_tW_top",Study,71,71)
    Make_CondorScr("DYJetsToLL_M_10To50",Study,178,100)
    Make_CondorScr("WJetsToLNu",Study,120,100)
    Make_CondorScr("TTbar_WJetToLNu",Study,22,22)
    Make_CondorScr("TTbar_WQQ",Study,4,4)
    Make_CondorScr("TTbar_ZToLLNuNu",Study,56,56)
    Make_CondorScr("TTbar_ZQQ",Study,4,4)
    Make_CondorScr("ZZ",Study,21,21)
    Make_CondorScr("WW",Study,81,81)
    Make_CondorScr("WZ",Study,41,41)
    Make_CondorScr("TTJets_Signal",Study,621,621)

    Make_CondorScr("Data_DoubleEG_Run2016B",Study,264,100)
    Make_CondorScr("Data_DoubleEG_Run2016C",Study,87,87)
    Make_CondorScr("Data_DoubleEG_Run2016D",Study,146,100)
    Make_CondorScr("Data_DoubleEG_Run2016E",Study,124,100)
    Make_CondorScr("Data_DoubleEG_Run2016F",Study,91,91)
    Make_CondorScr("Data_DoubleEG_Run2016G",Study,214,100)
    Make_CondorScr("Data_DoubleEG_Run2016HV2",Study,232,100)
    Make_CondorScr("Data_DoubleEG_Run2016HV3",Study,7,7)
    Make_CondorScr("Data_SingleElectron_Run2016B",Study,1055,100)
    Make_CondorScr("Data_SingleElectron_Run2016C",Study,348,87)
    Make_CondorScr("Data_SingleElectron_Run2016D",Study,584,100)
    Make_CondorScr("Data_SingleElectron_Run2016E",Study,496,100)
    Make_CondorScr("Data_SingleElectron_Run2016F",Study,359,91)
    Make_CondorScr("Data_SingleElectron_Run2016G",Study,854,100)
    Make_CondorScr("Data_SingleElectron_Run2016HV2",Study,925,100)
    Make_CondorScr("Data_SingleElectron_Run2016HV3",Study,25,7)
