# Update and Upgrade Ubuntu/Debian System
sudo apt-get update 
sudo apt-get upgrade 

# Install Python3, pip, and venv
sudo apt-get install python3.8
sudo apt-get install python3-pip
sudo apt-get install python3-venv 


python3 -m venv virtualEnvETL
source virtualEnvETL/bin/activate 
pip install -r requirements.txt 

deactivate 





