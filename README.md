Testing in CMSSW_13_3_1:

```
cd CMSSW_13_3_1/src
cmsenv
git clone https://github.com/YujiLee301/ScriptsJIT.git
mv ScriptsJIT/* .
cmsRun PhysicsTools/PatExamples/test/SSLPuppi_cfg.py
```
