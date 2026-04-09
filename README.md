Based on the KinFitter package by Thomas Goepfert:
KinFitter - A Kinematic Fit with Constraints

Wiki: https://github.com/goepfert/KinFitter/wiki/KinFitter---A-Kinematic-Fit-with-Constraints
Documentation: http://www.iktp.tu-dresden.de/~goepfert/KinFitter.pdf

# Setup:
Clone the repo:
```bash
git clone git@github.com:LuSchaller/pyKinFit.git
```
Then setup your python environment and do
```bash
cd pyKinFit
make
```
Edit the paths in path.sh or manually export the environment variables  
  
in Python:
```
import pyKinFit
pyKinFit.setBestCombi(inputpt, inputeta, inputphi, inputmass)
```
You can execute
```bash
python3 cfest.py
```
to confirm everything works
