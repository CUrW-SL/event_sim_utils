import csv
import traceback
import pymysql

from db_adapter.base import get_Pool, destroy_Pool
from db_adapter.constants import CURW_OBS_USERNAME, CURW_OBS_DATABASE, CURW_OBS_HOST, CURW_OBS_PASSWORD, CURW_OBS_PORT


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


def extract_active_curw_obs_rainfall_stations():
    """
        Extract currently active (active within last week) rainfall obs stations
        :return:
        """
    # Connect to the database
    pool = get_Pool(host=CURW_OBS_HOST, port=CURW_OBS_PORT, user=CURW_OBS_USERNAME, password=CURW_OBS_PASSWORD,
                    db=CURW_OBS_DATABASE)

    obs_stations = [['hash_id', 'station_id', 'station_name', 'latitude', 'longitude']]

    connection = pool.connection()

    try:

        with connection.cursor() as cursor1:
            cursor1.callproc(procname='getActiveRainfallObsStations')
            results = cursor1.fetchall()

            for result in results:
                obs_stations.append([result.get('hash_id'), result.get('station_id'), result.get('station_name'),
                                     result.get('latitude'), result.get('longitude')])

        # Write to csv file
        create_csv('grids/obs_stations/rainfall/curw_active_rainfall_obs_stations.csv', obs_stations)

    except Exception as ex:
        traceback.print_exc()
    finally:
        connection.close()
        destroy_Pool(pool)


extract_active_curw_obs_rainfall_stations()