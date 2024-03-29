#!/home/curw/event_sim_utils/venv/bin/python3
"""use pandas, replace into"""

import operator
import collections
import traceback
import os, getopt, sys
import pandas as pd
import numpy as np
import json
from datetime import datetime, timedelta

from math import acos, cos, sin, radians

from db_adapter.csv_utils import read_csv

from db_adapter.constants import set_db_config_file_path
from db_adapter.constants import connection as con_params
from db_adapter.base import get_Pool, destroy_Pool
from db_adapter.curw_sim.timeseries import Timeseries as Sim_Timeseries
from db_adapter.curw_sim.grids import GridInterpolationEnum
from db_adapter.curw_sim.timeseries import MethodEnum
from db_adapter.curw_sim.common import extract_obs_rain_5_min_ts
from db_adapter.curw_sim.constants import FLO2D_250, FLO2D_150, FLO2D_150_V2

from db_adapter.curw_sim.common import DelTimeseries, get_curw_sim_hash_ids

from db_adapter.logger import logger


DATE_TIME_FORMAT = '%Y-%m-%d %H:%M:%S'
ROOT_DIR = '/home/curw/event_sim_utils'


def check_time_format(time, model):
    try:
        time = datetime.strptime(time, DATE_TIME_FORMAT)

        if time.strftime('%S') != '00':
            print("Seconds should be always 00")
            exit(1)
        if model=="flo2d_250" and time.strftime('%M') not in ('05', '10', '15', '20', '25', '30', '35', '40', '45', '50', '55', '00'):
            print("Minutes should be multiple of 5 fro flo2d_250")
            exit(1)
        if model in ("flo2d_150", "flo2d_150_v2") and time.strftime('%M') not in ('15', '30', '45', '00'):
            print("Minutes should be multiple of 15 for flo2d_150 and flo2d_150_v2")
            exit(1)

        return True
    except Exception:
        traceback.print_exc()
        print("Time {} is not in proper format".format(time))
        exit(1)


def list_of_lists_to_df_first_row_as_columns(data):
    return pd.DataFrame.from_records(data[1:], columns=data[0])


def clear_database(curw_sim_pool, method, model):
    print("clearing  database ... ")
    run_table = "run"
    data_table = "data"

    hash_ids = get_curw_sim_hash_ids(pool=curw_sim_pool, run_table=run_table, model=model, method=method,
                                     obs_end_start=None, obs_end_end=None, grid_id=None)
    DTS = DelTimeseries(pool=curw_sim_pool, run_table=run_table, data_table=data_table)
    DTS.bulk_delete_all_by_hash_id(ids=hash_ids)

    print("database cleared !")


# extract curw active rainfall stations within a given perios
def extract_active_curw_obs_rainfall_stations(curw_obs_pool, start_time, end_time):
    """
        Extract currently active (active within last week) rainfall obs stations
        :return:
        """

    obs_stations = [['hash_id', 'station_id', 'station_name', 'latitude', 'longitude']]

    connection = curw_obs_pool.connection()

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
        destroy_Pool(curw_obs_pool)


# map nearest observational stations to flo2d grids
def find_nearest_obs_stations_for_flo2d_stations(flo2d_stations_csv, obs_stations, flo2d_model):
    # obs_stations : [hash_id,station_id,station_name,latitude,longitude]

    flo2d_station = read_csv(flo2d_stations_csv)  # [Grid_ID,X,Y]

    flo2d_obs_mapping_dict = {}

    for flo2d_index in range(len(flo2d_station)):

        grid_id = flo2d_station[flo2d_index][0]

        flo2d_obs_mapping = []

        flo2d_lat = float(flo2d_station[flo2d_index][2])
        flo2d_lng = float(flo2d_station[flo2d_index][1])

        distances = {}

        for obs_index in range(len(obs_stations)):
            lat = float(obs_stations[obs_index][3])
            lng = float(obs_stations[obs_index][4])

            intermediate_value = cos(radians(flo2d_lat)) * cos(radians(lat)) * cos(
                radians(lng) - radians(flo2d_lng)) + sin(radians(flo2d_lat)) * sin(radians(lat))
            if intermediate_value < 1:
                distance = 6371 * acos(intermediate_value)
            else:
                distance = 6371 * acos(1)

            distances[obs_stations[obs_index][1]] = distance

        sorted_distances = collections.OrderedDict(sorted(distances.items(), key=operator.itemgetter(1))[:10])

        count = 0
        for key in sorted_distances.keys():
            if count < 3 and sorted_distances.get(key) <= 25:
                flo2d_obs_mapping.append(key)
                count += 1
            else:
                break

        # print(flo2d_obs_mapping)
        flo2d_obs_mapping_dict[grid_id] = flo2d_obs_mapping

    # flo2d_grid_mappings[dict.get("grid_id")] = [dict.get("obs1"), dict.get("obs2"), dict.get("obs3")]
    flo2d_grid_mappings_dict = {}

    return flo2d_obs_mapping_dict

# for bulk insertion for a given one grid interpolation method
def update_rainfall_obs(curw_obs_pool, curw_sim_pool, flo2d_model, method, grid_interpolation, timestep, start_time, end_time):

    """
    Update rainfall observations for flo2d models
    :param flo2d_model: flo2d model
    :param method: value interpolation method
    :param grid_interpolation: grid interpolation method
    :param timestep: output timeseries timestep
    :return:
    """

    # obs_start = datetime.strptime(start_time, '%Y-%m-%d %H:%M:%S')

    try:

        curw_obs_connection = curw_obs_pool.connection()

        # [hash_id, station_id, station_name, latitude, longitude]
        # active_obs_stations = read_csv(os.path.join(ROOT_DIR,'grids/obs_stations/rainfall/curw_active_rainfall_obs_stations.csv'))
        active_obs_stations = extract_active_curw_obs_rainfall_stations(curw_obs_pool=curw_obs_pool, start_time=start_time, end_time=end_time)[1:]
        flo2d_grids = read_csv(os.path.join(ROOT_DIR,'grids/flo2d/{}m.csv'.format(flo2d_model)))  # [Grid_ ID, X(longitude), Y(latitude)]

        stations_dict_for_obs = { }  # keys: obs station id , value: hash id

        for obs_index in range(len(active_obs_stations)):
            stations_dict_for_obs[active_obs_stations[obs_index][1]] = active_obs_stations[obs_index][0]

        # flo2d_obs_mapping = get_flo2d_cells_to_obs_grid_mappings(pool=curw_sim_pool, grid_interpolation=grid_interpolation, flo2d_model=flo2d_model)
        flo2d_obs_mapping = find_nearest_obs_stations_for_flo2d_stations(
            flo2d_stations_csv=os.path.join(ROOT_DIR,'grids/flo2d/{}m.csv'.format(flo2d_model)),
            obs_stations=active_obs_stations, flo2d_model=flo2d_model)

        # retrieve observed timeseries
        obs_df = pd.DataFrame()
        obs_df['time'] = pd.date_range(start=start_time, end=end_time, freq='5min')

        for obs_id in stations_dict_for_obs.keys():
            ts = extract_obs_rain_5_min_ts(connection=curw_obs_connection, start_time=start_time,
                                           id=stations_dict_for_obs.get(obs_id), end_time=end_time)
            ts.insert(0, ['time', obs_id])
            ts_df = list_of_lists_to_df_first_row_as_columns(ts)
            ts_df[obs_id] = ts_df[obs_id].astype('float64')

            obs_df = pd.merge(obs_df, ts_df, how="left", on='time')

        obs_df.set_index('time', inplace=True)
        obs_df['0'] = 0
        if timestep == 15:
            obs_df = obs_df.resample('15min', label='right', closed='right').sum()

        TS = Sim_Timeseries(pool=curw_sim_pool)

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

            print(datetime.now().strftime(DATE_TIME_FORMAT))
            print("grid_id:", cell_id)
            obs_station_ids = flo2d_obs_mapping.get(cell_id)

            if len(obs_station_ids) == 1:
                obs_ts_df = obs_df[obs_station_ids].to_frame(name='final')
            elif len(obs_station_ids) == 2:
                obs_ts_df = obs_df[obs_station_ids]
                obs_ts_df[obs_station_ids[0]] = obs_ts_df[obs_station_ids[0]].fillna(obs_ts_df[obs_station_ids[1]])
                obs_ts_df['final'] = obs_ts_df[obs_station_ids[0]]
            elif len(obs_station_ids) == 3:
                obs_ts_df = obs_df[obs_station_ids]
                obs_ts_df[obs_station_ids[1]] = obs_ts_df[obs_station_ids[1]].fillna(obs_ts_df[obs_station_ids[2]])
                obs_ts_df[obs_station_ids[0]] = obs_ts_df[obs_station_ids[0]].fillna(obs_ts_df[obs_station_ids[1]])
                obs_ts_df['final'] = obs_ts_df[obs_station_ids[0]]
            else:
                obs_ts_df = obs_df['0'].to_frame(name='final')

            final_ts_df = obs_ts_df['final'].reset_index()
            final_ts_df['time'] = final_ts_df['time'].dt.strftime(DATE_TIME_FORMAT)
            final_ts = final_ts_df.values.tolist()

            if final_ts is not None and len(final_ts) > 0:
                TS.replace_data(timeseries=final_ts, tms_id=tms_id)
                TS.update_latest_obs(id_=tms_id, obs_end=(final_ts[-1][1]))

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
    -------------------------------------------------
    Populate rainfall Flo2D 250, 150 & 150_v2 :: OBS
    -------------------------------------------------

    Usage: ./rain/flo2d_OBS.py [-m flo2d_XXX][-s "YYYY-MM-DD HH:MM:SS"] [-e "YYYY-MM-DD HH:MM:SS"] [-g XXXX]

    -h  --help          Show usage
    -m  --flo2d_model   FLO2D model (e.g. flo2d_250, flo2d_150). Default is flo2d_250.
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
        # os.system("{}/grids/obs_stations/rainfall/update_active_curw_rainfall_stations.py -s {} -e {}"
        #           .format(ROOT_DIR, start_time, end_time))

        # prepare and populate flo2d grid maps
        # 1. flo2d grids to weather stations
        # 2. flo2d grids to d03 stations
        # os.system("{}/grid_maps/flo2d/update_flo2d_grid_maps.py -m {} -g {}".format(ROOT_DIR, flo2d_model, grid_interpolation))

        curw_obs_pool = get_Pool(host=con_params.CURW_OBS_HOST, user=con_params.CURW_OBS_USERNAME,
                                 password=con_params.CURW_OBS_PASSWORD,
                                 port=con_params.CURW_OBS_PORT, db=con_params.CURW_OBS_DATABASE)

        curw_sim_pool = get_Pool(host=con_params.CURW_SIM_HOST, user=con_params.CURW_SIM_USERNAME,
                                 password=con_params.CURW_SIM_PASSWORD,
                                 port=con_params.CURW_SIM_PORT, db=con_params.CURW_SIM_DATABASE)

        print("{} : ####### Insert obs rainfall for {} grids".format(datetime.now(), flo2d_model))
        update_rainfall_obs(curw_obs_pool=curw_obs_pool, curw_sim_pool=curw_sim_pool,
                            flo2d_model=flo2d_model, method=method, grid_interpolation=grid_interpolation,
                            timestep=timestep, start_time=start_time, end_time=end_time)
    except Exception as e:
        traceback.print_exc()
    finally:
        print("Process finished.")


