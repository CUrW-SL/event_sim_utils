#!/home/curw/event_sim_utils/venv/bin/python3

import csv
import traceback
import os
import sys
import getopt
from datetime import datetime, timedelta

from db_adapter.base import get_Pool, destroy_Pool
from db_adapter.constants import set_db_config_file_path
from db_adapter.constants import connection as con_params

ROOT_DIR = '/home/curw/event_sim_utils'


def create_csv(file_name, data):
    """
    Create new csv file using given data
    :param file_name: <file_path/file_name>.csv
    :param data: list of lists
    e.g. [['Person', 'Age'], ['Peter', '22'], ['Jasmine', '21'], ['Sam', '24']]
    :return:
    """
    with open(file_name, 'w') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerows(data)


def extract_active_curw_obs_rainfall_stations(start_time, end_time):
    """
        Extract currently active (active within last week) rainfall obs stations
        :return:
        """
    # Connect to the database
    pool = get_Pool(host=con_params.CURW_OBS_HOST, port=con_params.CURW_OBS_PORT, user=con_params.CURW_OBS_USERNAME,
                    password=con_params.CURW_OBS_PASSWORD, db=con_params.CURW_OBS_DATABASE)

    obs_stations = [['hash_id', 'station_id', 'station_name', 'latitude', 'longitude']]

    connection = pool.connection()

    try:

        with connection.cursor() as cursor1:
            cursor1.callproc('getActiveRfStationsAtGivenTime', (start_time, end_time))
            results = cursor1.fetchall()

            for result in results:
                obs_stations.append([result.get('hash_id'), result.get('station_id'), result.get('station_name'),
                                     result.get('latitude'), result.get('longitude')])

        # Write to csv file
        create_csv(os.path.join(ROOT_DIR,'grids/obs_stations/rainfall/curw_active_rainfall_obs_stations.csv'), obs_stations)

    except Exception as ex:
        traceback.print_exc()
    finally:
        connection.close()
        destroy_Pool(pool)


def usage():
    usageText = """
    ----------------------------------------------------------
    Find active rainfall observation stations at a given time
    ----------------------------------------------------------

    Usage: ./grids/obs_stations/rainfall/update_active_curw_rainfall_stations.py [-s "YYYY-MM-DD HH:MM:SS"] [-e "YYYY-MM-DD HH:MM:SS"]

    -h  --help          Show usage
    -s  --start_time    Rain timeseries start time (e.g: "2019-06-05 00:00:00"). Default is 00:00:00, yesterday.
    -e  --end_time      Rain timeseries end time (e.g: "2019-06-05 23:30:00"). Default is 00:00:00, tomorrow.
    """
    print(usageText)


if __name__=="__main__":

    set_db_config_file_path(os.path.join(ROOT_DIR, 'db_adapter_config.json'))

    try:

        start_time = None
        end_time = None

        try:
            opts, args = getopt.getopt(sys.argv[1:], "h:s:e:",
                                       ["help", "start_time=", "end_time="])
        except getopt.GetoptError:
            usage()
            sys.exit(2)
        for opt, arg in opts:
            if opt in ("-h", "--help"):
                usage()
                sys.exit()
            elif opt in ("-s", "--start_time"):
                start_time = arg.strip()
            elif opt in ("-e", "--end_time"):
                end_time = arg.strip()

        if start_time is None:
            start_time = (datetime.now() - timedelta(days=1)).strftime('%Y-%m-%d 00:00:00')

        if end_time is None:
            end_time = (datetime.now() + timedelta(days=1)).strftime('%Y-%m-%d 00:00:00')

        extract_active_curw_obs_rainfall_stations(start_time, end_time)

    except Exception as e:
        traceback.print_exc()
    finally:
        print("Process finished.")
