#!/home/curw/event_sim_utils/venv/bin/python3

import traceback
import os, getopt, sys
from datetime import datetime, timedelta

from db_adapter.csv_utils import read_csv
from db_adapter.base import get_Pool, destroy_Pool
from db_adapter.constants import set_db_config_file_path
from db_adapter.constants import connection as con_params
from db_adapter.curw_sim.timeseries import Timeseries
from db_adapter.curw_sim.constants import HecHMS
from db_adapter.curw_sim.grids import GridInterpolationEnum
from db_adapter.curw_sim.timeseries import MethodEnum
from db_adapter.curw_sim.common import process_5_min_ts, process_15_min_ts, \
    extract_obs_rain_5_min_ts, extract_obs_rain_15_min_ts
from db_adapter.logger import logger

ROOT_DIR = '/home/curw/event_sim_utils'
DATE_TIME_FORMAT = '%Y-%m-%d %H:%M:%S'


def check_time_format(time):
    try:
        time = datetime.strptime(time, DATE_TIME_FORMAT)

        if time.strftime('%S') != '00':
            print("Seconds should be always 00")
            exit(1)
        if time.strftime('%M') not in ('05', '10', '15', '20', '25', '30', '35', '40', '45', '50', '55', '00'):
            print("Minutes should be multiple of 5.")
            exit(1)

        return True
    except Exception:
        traceback.print_exc()
        print("Time {} is not in proper format".format(time))
        exit(1)


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

        return obs_stations

    except Exception as ex:
        traceback.print_exc()
    finally:
        connection.close()
        destroy_Pool(pool)


# for bulk insertion for a given one grid interpolation method
def update_rainfall_obs(target_model, method, timestep, start_time, end_time):

    """
    Update rainfall observations for flo2d models
    :param model: target model
    :param method: value interpolation method
    :param grid_interpolation: grid interpolation method
    :param timestep:
    :return:
    """
    obs_start = datetime.strptime(start_time, '%Y-%m-%d %H:%M:%S')
    try:

        # Connect to the database
        curw_obs_pool = get_Pool(host=con_params.CURW_OBS_HOST, user=con_params.CURW_OBS_USERNAME,
                                 password=con_params.CURW_OBS_PASSWORD, port=con_params.CURW_OBS_PORT,
                                 db=con_params.CURW_OBS_DATABASE)

        curw_obs_connection = curw_obs_pool.connection()

        curw_sim_pool = get_Pool(host=con_params.CURW_SIM_HOST, user=con_params.CURW_SIM_USERNAME,
                                 password=con_params.CURW_SIM_PASSWORD, port=con_params.CURW_SIM_PORT,
                                 db=con_params.CURW_SIM_DATABASE)

        TS = Timeseries(pool=curw_sim_pool)

        # [hash_id, station_id, station_name, latitude, longitude]
        active_obs_stations = extract_active_curw_obs_rainfall_stations(start_time=start_time, end_time=end_time)[1:]
        obs_stations_dict = { }  # keys: obs station id , value: [hash id, name, latitude, longitude]

        for obs_index in range(len(active_obs_stations)):
            obs_stations_dict[active_obs_stations[obs_index][1]] = [active_obs_stations[obs_index][0],
                                                                    active_obs_stations[obs_index][2],
                                                                    active_obs_stations[obs_index][3],
                                                                    active_obs_stations[obs_index][4]]

        for obs_id in obs_stations_dict.keys():
            meta_data = {
                    'latitude': float('%.6f' % float(obs_stations_dict.get(obs_id)[2])),
                    'longitude': float('%.6f' % float(obs_stations_dict.get(obs_id)[3])),
                    'model': target_model, 'method': method,
                    'grid_id': 'rainfall_{}_{}'.format(obs_id, obs_stations_dict.get(obs_id)[1])
                    }

            tms_id = TS.get_timeseries_id_if_exists(meta_data=meta_data)

            if tms_id is None:
                tms_id = TS.generate_timeseries_id(meta_data=meta_data)
                meta_data['id'] = tms_id
                TS.insert_run(meta_data=meta_data)

            TS.update_grid_id(id_=tms_id, grid_id=meta_data['grid_id'])

            obs_hash_id = obs_stations_dict.get(obs_id)[0]

            obs_timeseries = []

            if timestep == 5:
                ts = extract_obs_rain_5_min_ts(connection=curw_obs_connection, start_time=obs_start, end_time=end_time,
                                               id=obs_hash_id)
                if ts is not None and len(ts) > 1:
                    obs_timeseries.extend(process_5_min_ts(newly_extracted_timeseries=ts, expected_start=obs_start)[1:])
                    # obs_start = ts[-1][0]
            elif timestep == 15:
                ts = extract_obs_rain_15_min_ts(connection=curw_obs_connection, start_time=obs_start, end_time=end_time,
                                                id=obs_hash_id)
                if ts is not None and len(ts) > 1:
                    obs_timeseries.extend(process_15_min_ts(newly_extracted_timeseries=ts, expected_start=obs_start)[1:])
                    # obs_start = ts[-1][0]

            # for i in range(len(obs_timeseries)):
            #     if obs_timeseries[i][1] == -99999:
            #         obs_timeseries[i][1] = 0

            if obs_timeseries is not None and len(obs_timeseries) > 0:
                TS.insert_data(timeseries=obs_timeseries, tms_id=tms_id, upsert=True)

    except Exception as e:
        traceback.print_exc()
        logger.error("Exception occurred while updating obs rainfalls in curw_sim.")
    finally:
        curw_obs_connection.close()
        destroy_Pool(pool=curw_sim_pool)
        destroy_Pool(pool=curw_obs_pool)


def usage():
    usageText = """
    ---------------------------------
    Populate rainfall HecHMS :: OBS
    ---------------------------------

    Usage: ./rain/hechms_OBS.py [-s "YYYY-MM-DD HH:MM:SS"] [-e "YYYY-MM-DD HH:MM:SS"]

    -h  --help          Show usage
    -s  --start_time    Rain timeseries start time (e.g: "2019-06-05 00:00:00"). Default is 23:30:00, 3 days before today.
    -e  --end_time      Rain timeseries end time (e.g: "2019-06-05 23:30:00"). Default is 23:30:00, tomorrow.
    """
    print(usageText)


if __name__=="__main__":

    set_db_config_file_path(os.path.join(ROOT_DIR, 'db_adapter_config.json'))

    try:
        start_time = None
        end_time = None
        model = "hechms"
        method = "OBS"

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
            start_time = (datetime.now() - timedelta(days=3)).strftime('%Y-%m-%d 23:30:00')
        else:
            check_time_format(time=start_time)

        if end_time is None:
            end_time = (datetime.now() + timedelta(days=1)).strftime('%Y-%m-%d 23:30:00')
        else:
            check_time_format(time=end_time)


        print("{} : ####### Insert grid based obs rainfall for hechms.".format(datetime.now()))
        update_rainfall_obs(target_model=HecHMS, method=method, timestep=5,
                            start_time=start_time, end_time=end_time)

    except Exception as e:
        traceback.print_exc()
    finally:
        print("Process finished.")