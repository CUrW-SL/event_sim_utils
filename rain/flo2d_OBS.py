#!/home/uwcc-admin/event_sim_utils/venv/bin/python3

import sys
import getopt
import os
import traceback
from datetime import datetime, timedelta

from db_adapter.csv_utils import read_csv

from db_adapter.constants import set_db_config_file_path
from db_adapter.constants import connection as con_params
from db_adapter.base import get_Pool, destroy_Pool
from db_adapter.curw_sim.grids import get_flo2d_cells_to_wrf_grid_mappings, get_flo2d_cells_to_obs_grid_mappings
from db_adapter.curw_sim.timeseries import Timeseries as Sim_Timeseries
from db_adapter.curw_sim.common import convert_15_min_ts_to_5_mins_ts, append_value_for_timestamp, summed_timeseries, \
    process_5_min_ts, process_15_min_ts, fill_missing_values, \
    extract_obs_rain_5_min_ts, extract_obs_rain_15_min_ts
from db_adapter.curw_sim.grids import GridInterpolationEnum
from db_adapter.curw_sim.timeseries import MethodEnum
from db_adapter.curw_sim.constants import FLO2D_250, FLO2D_150, FLO2D_150_V2

from db_adapter.logger import logger


DATE_TIME_FORMAT = '%Y-%m-%d %H:%M:%S'
ROOT_DIR = '/home/uwcc-admin/event_sim_utils'


def check_time_format(time, model):
    try:
        time = datetime.strptime(time, DATE_TIME_FORMAT)

        if time.strftime('%S') != '00':
            print("Seconds should be always 00")
            exit(1)
        if model=="flo2d_250" and time.strftime('%M') not in ('05', '10', '15', '20', '25', '30', '35', '40', '45', '50', '55', '00'):
            print("Minutes should be multiple of 5 fro flo2d_250")
            exit(1)
        if model=="flo2d_150" and time.strftime('%M') not in ('15', '30', '45', '00'):
            print("Minutes should be multiple of 15 for flo2d_150")
            exit(1)

        return True
    except Exception:
        traceback.print_exc()
        print("Time {} is not in proper format".format(time))
        exit(1)


# for bulk insertion for a given one grid interpolation method
def update_rainfall_obs(flo2d_model, method, grid_interpolation, timestep, start_time, end_time):

    """
    Update rainfall observations for flo2d models
    :param flo2d_model: flo2d model
    :param method: value interpolation method
    :param grid_interpolation: grid interpolation method
    :param timestep: output timeseries timestep
    :return:
    """

    obs_start = datetime.strptime(start_time, '%Y-%m-%d %H:%M:%S')

    try:

        # Connect to the database
        curw_obs_pool = get_Pool(host=con_params.CURW_OBS_HOST, user=con_params.CURW_OBS_USERNAME,
                                 password=con_params.CURW_OBS_PASSWORD,
                                 port=con_params.CURW_OBS_PORT, db=con_params.CURW_OBS_DATABASE)

        curw_obs_connection = curw_obs_pool.connection()

        curw_sim_pool = get_Pool(host=con_params.CURW_SIM_HOST, user=con_params.CURW_SIM_USERNAME,
                                 password=con_params.CURW_SIM_PASSWORD,
                                 port=con_params.CURW_SIM_PORT, db=con_params.CURW_SIM_DATABASE)

        TS = Sim_Timeseries(pool=curw_sim_pool)

        # [hash_id, station_id, station_name, latitude, longitude]
        active_obs_stations = read_csv('grids/obs_stations/rainfall/curw_active_rainfall_obs_stations.csv')
        flo2d_grids = read_csv('grids/flo2d/{}m.csv'.format(flo2d_model))  # [Grid_ ID, X(longitude), Y(latitude)]

        stations_dict_for_obs = { }  # keys: obs station id , value: hash id

        for obs_index in range(len(active_obs_stations)):
            stations_dict_for_obs[active_obs_stations[obs_index][1]] = active_obs_stations[obs_index][0]

        flo2d_obs_mapping = get_flo2d_cells_to_obs_grid_mappings(pool=curw_sim_pool, grid_interpolation=grid_interpolation, flo2d_model=flo2d_model)

        for flo2d_index in range(len(flo2d_grids)):
            lat = flo2d_grids[flo2d_index][2]
            lon = flo2d_grids[flo2d_index][1]
            cell_id = flo2d_grids[flo2d_index][0]
            meta_data = {
                    'latitude': float('%.6f' % float(lat)), 'longitude': float('%.6f' % float(lon)),
                    'model': flo2d_model, 'method': method,
                    'grid_id': '{}_{}_{}'.format(flo2d_model, grid_interpolation, (str(cell_id)).zfill(10))
                    }

            tms_id = TS.get_timeseries_id(grid_id=meta_data.get('grid_id'), method=meta_data.get('method'))

            if tms_id is None:
                tms_id = TS.generate_timeseries_id(meta_data=meta_data)
                meta_data['id'] = tms_id
                TS.insert_run(meta_data=meta_data)

            obs1_station_id = str(flo2d_obs_mapping.get(meta_data['grid_id'])[0])
            obs2_station_id = str(flo2d_obs_mapping.get(meta_data['grid_id'])[1])
            obs3_station_id = str(flo2d_obs_mapping.get(meta_data['grid_id'])[2])

            obs_timeseries = []

            if timestep == 5:
                if obs1_station_id != str(-1):
                    obs1_hash_id = stations_dict_for_obs.get(obs1_station_id)

                    ts = extract_obs_rain_5_min_ts(connection=curw_obs_connection, start_time=obs_start, id=obs1_hash_id,
                                                   end_time=end_time)
                    if ts is not None and len(ts) > 1:
                        obs_timeseries.extend(process_5_min_ts(newly_extracted_timeseries=ts, expected_start=obs_start)[1:])
                        # obs_start = ts[-1][0]

                    if obs2_station_id != str(-1):
                        obs2_hash_id = stations_dict_for_obs.get(obs2_station_id)

                        ts2 = extract_obs_rain_5_min_ts(connection=curw_obs_connection, start_time=obs_start, id=obs2_hash_id,
                                                        end_time=end_time)
                        if ts2 is not None and len(ts2) > 1:
                            obs_timeseries = fill_missing_values(newly_extracted_timeseries=ts2, OBS_TS=obs_timeseries)
                            if obs_timeseries is not None and len(obs_timeseries) > 0:
                                expected_start = obs_timeseries[-1][0]
                            else:
                                expected_start= obs_start
                            obs_timeseries.extend(process_5_min_ts(newly_extracted_timeseries=ts2, expected_start=expected_start)[1:])
                            # obs_start = ts2[-1][0]

                        if obs3_station_id != str(-1):
                            obs3_hash_id = stations_dict_for_obs.get(obs3_station_id)

                            ts3 = extract_obs_rain_5_min_ts(connection=curw_obs_connection, start_time=obs_start, id=obs3_hash_id,
                                                            end_time=end_time)
                            if ts3 is not None and len(ts3) > 1 and len(obs_timeseries) > 0:
                                obs_timeseries = fill_missing_values(newly_extracted_timeseries=ts3, OBS_TS=obs_timeseries)
                                if obs_timeseries is not None:
                                    expected_start = obs_timeseries[-1][0]
                                else:
                                    expected_start= obs_start
                                obs_timeseries.extend(process_5_min_ts(newly_extracted_timeseries=ts3, expected_start=expected_start)[1:])
            elif timestep == 15:
                if obs1_station_id != str(-1):
                    obs1_hash_id = stations_dict_for_obs.get(obs1_station_id)

                    ts = extract_obs_rain_15_min_ts(connection=curw_obs_connection, start_time=obs_start, id=obs1_hash_id,
                                                    end_time=end_time)
                    if ts is not None and len(ts) > 1:
                        obs_timeseries.extend(process_15_min_ts(newly_extracted_timeseries=ts, expected_start=obs_start)[1:])
                        # obs_start = ts[-1][0]

                    if obs2_station_id != str(-1):
                        obs2_hash_id = stations_dict_for_obs.get(obs2_station_id)

                        ts2 = extract_obs_rain_15_min_ts(connection=curw_obs_connection, start_time=obs_start, id=obs2_hash_id,
                                                         end_time=end_time)
                        if ts2 is not None and len(ts2) > 1:
                            obs_timeseries = fill_missing_values(newly_extracted_timeseries=ts2, OBS_TS=obs_timeseries)
                            if obs_timeseries is not None and len(obs_timeseries) > 0:
                                expected_start = obs_timeseries[-1][0]
                            else:
                                expected_start = obs_start
                            obs_timeseries.extend(process_15_min_ts(newly_extracted_timeseries=ts2, expected_start=expected_start)[1:])
                            # obs_start = ts2[-1][0]

                        if obs3_station_id != str(-1):
                            obs3_hash_id = stations_dict_for_obs.get(obs3_station_id)

                            ts3 = extract_obs_rain_15_min_ts(connection=curw_obs_connection, start_time=obs_start, id=obs3_hash_id,
                                                             end_time=end_time)
                            if ts3 is not None and len(ts3) > 1 and len(obs_timeseries) > 0:
                                obs_timeseries = fill_missing_values(newly_extracted_timeseries=ts3, OBS_TS=obs_timeseries)
                                if obs_timeseries is not None:
                                    expected_start = obs_timeseries[-1][0]
                                else:
                                    expected_start = obs_start
                                obs_timeseries.extend(process_15_min_ts(newly_extracted_timeseries=ts3, expected_start=expected_start)[1:])

            for i in range(len(obs_timeseries)):
                if obs_timeseries[i][1] == -99999:
                    obs_timeseries[i][1] = 0

            if obs_timeseries is not None and len(obs_timeseries) > 0:
                TS.insert_data(timeseries=obs_timeseries, tms_id=tms_id, upsert=True)

    except Exception as e:
        traceback.print_exc()
        logger.error("Exception occurred while updating obs rainfalls in curw_sim.")
    finally:
        curw_obs_connection.close()
        destroy_Pool(pool=curw_sim_pool)
        destroy_Pool(pool=curw_obs_pool)
        logger.info("Process finished")


def usage():
    usageText = """
    -----------------------------------------
    Populate rainfall Flo2D 250 & 150 :: OBS
    -----------------------------------------

    Usage: ./rain/flo2d_OBS.py [-m flo2d_XXX][-s "YYYY-MM-DD HH:MM:SS"] [-e "YYYY-MM-DD HH:MM:SS"] [-g XXXX]

    -h  --help          Show usage
    -m  --model         FLO2D model (e.g. flo2d_250, flo2d_150). Default is flo2d_250.
    -s  --start_time    Rain timeseries start time (e.g: "2019-06-05 00:00:00"). Default is 23:30:00, 3 days before today.
    -e  --end_time      Rain timeseries end time (e.g: "2019-06-05 23:30:00"). Default is 23:30:00, tomorrow.
    -g  --grid_tag      Grid mapping method (e.g: "MDPA", "TP"). Default is "MDPA".
    """
    print(usageText)


if __name__=="__main__":

    set_db_config_file_path(os.path.join(ROOT_DIR, 'db_adapter_config.json'))

    try:
        start_time = None
        end_time = None
        flo2d_model = None
        # method = "OBS"
        grid_tag = "MDPA"  # note - grid tag "TP" has not handled yet

        try:
            opts, args = getopt.getopt(sys.argv[1:], "h:m:s:e:g:",
                                       ["help", "flo2d_model=", "start_time=", "end_time=", "grid_tag="])
        except getopt.GetoptError:
            usage()
            sys.exit(2)
        for opt, arg in opts:
            if opt in ("-h", "--help"):
                usage()
                sys.exit()
            elif opt in ("-m", "--flo2d_model"):
                flo2d_model = arg.strip()
            elif opt in ("-s", "--start_time"):
                start_time = arg.strip()
            elif opt in ("-e", "--end_time"):
                end_time = arg.strip()
            elif opt in ("-g", "--grid_tag"):
                grid_tag = arg.strip()

        if grid_tag == "MDPA":
            grid_interpolation = GridInterpolationEnum.getAbbreviation(GridInterpolationEnum.MDPA)
            method = MethodEnum.getAbbreviation(MethodEnum.OBS)
        else:
            exit(0)

        if flo2d_model is None:
            print("Flo2d model is not specified.")
            exit(1)
        elif flo2d_model not in ("flo2d_250", "flo2d_150", "flo2d_150_v2"):
            print("Flo2d model should be either \"flo2d_250\" or \"flo2d_150\" or \"flo2d_150_v2\"")
            exit(1)

        if start_time is None:
            start_time = (datetime.now() - timedelta(days=3)).strftime('%Y-%m-%d 23:30:00')
        else:
            check_time_format(time=start_time, model=flo2d_model)

        if end_time is None:
            end_time = (datetime.now() + timedelta(days=1)).strftime('%Y-%m-%d 23:30:00')
        else:
            check_time_format(time=end_time, model=flo2d_model)

        if flo2d_model == FLO2D_250:
            timestep = 5
        elif flo2d_model == FLO2D_150:
            timestep = 15
        elif flo2d_model == FLO2D_150_V2:
            timestep = 15

        # find actove curw weather stations during the specified time window
        os.system("./grids/obs_stations/rainfall/update_active_curw_rainfall_stations.py -s {} -e {}"
                  .format(start_time, end_time))

        # prepare and populate flo2d grid maps
        # 1. flo2d grids to weather stations
        # 2. flo2d grids to d03 stations
        os.system("./grid_maps/flo2d/update_flo2d_grid_maps.py -m {} -g {}".format(flo2d_model, grid_interpolation))

        print("{} : ####### Insert obs rainfall for FLO2D 250 grids".format(datetime.now()))
        update_rainfall_obs(flo2d_model=flo2d_model, method=method, grid_interpolation=grid_interpolation,
                            timestep=timestep, start_time=start_time, end_time=end_time)
    except Exception as e:
        traceback.print_exc()
    finally:
        print("Process finished.")

