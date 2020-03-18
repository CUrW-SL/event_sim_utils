#!/usr/bin/env bash

echo `date`

echo "Changing into ~/event_sim_utils"
cd /home/uwcc-admin/event_sim_utils
echo "Inside `pwd`"


# If no venv (python3 virtual environment) exists, then create one.
if [ ! -d "venv" ]
then
    echo "Creating venv python3 virtual environment."
    virtualenv -p python3 venv
fi

# Activate venv.
echo "Activating venv python3 virtual environment."
source venv/bin/activate

# Install dependencies using pip.
if [ ! -f "db.log" ]
then
    echo "Installing PyMySQL"
    pip install PyMySQL
    echo "Installing db adapter"
    pip install git+https://github.com/shadhini/curw_db_adapter.git
fi

# Deactivating virtual environment
echo "Deactivating virtual environment"
deactivate
