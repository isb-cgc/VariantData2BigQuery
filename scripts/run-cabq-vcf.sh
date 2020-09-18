#!/usr/bin/env bash

# Copyright 2019, Institute for Systems Biology
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.




# source ~/setEnvVars.sh

export MY_VENV=~/virtualEnvETL
export PYTHONPATH=.:${MY_VENV}/lib:~/extlib


# Directory created to store the intemediary files 
mkdir -p ~/NextGenETL/intermediateFiles 

source bin/activate
popd > /dev/null
cd ..
python3 /Users/oshahzada98/Desktop/NextGenETL/BQ_Table_Building/build_cabq_vcf.py /Users/oshahzada98/Desktop/NextGenETL/ConfigFiles/VcfCBQBuild.yml

deactivate
