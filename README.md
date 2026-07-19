# LEVEL_MEASURMENT
A simple program to calculate straightness measurements using a precision level

When rebuilding machine tools or other tasks that require careful linear measurment this can be done inexpensively with an accurate level (around 0,02 mm/m) 

This simple program lets you enter the measurments are they are made and calculates the deviations etc. 

This is 100% Claude generated, if anyone wants to make inprovments, developments etc that is most welcome

I removed the EXE file, it turns out that if it is downloaded it will be blocked by the windows anti-virus. 

I would suggest installing python 3, install pip, then numpy matplotlib and running the python version. You could also convert the python to exe locally and that will run just fine.  
One plus of having python installed is you can  make your own scripts.

# 1. Install Python (close and REOPEN the terminal afterwards so it's found)
winget install -e --id Python.Python.3.13

# 2. Confirm it installed
python --version

# 3. Install the two required libraries (tkinter is already built in)
python -m pip install --upgrade pip
python -m pip install numpy matplotlib

# 4. Download the script into your Downloads folder
cd $HOME\Downloads
Invoke-WebRequest -Uri "https://raw.githubusercontent.com/lukerickert/LEVEL_MEASURMENT/main/level_straightness.py" -OutFile "level_straightness.py"

# 5. Run it
python level_straightness.py
