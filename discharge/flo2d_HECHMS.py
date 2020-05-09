#!/home/curw/event_sim_utils/venv/bin/python3

import sys
import getopt
import os
import traceback
import pandas as pd
from datetime import datetime, timedelta

from db_adapter.csv_utils import read_csv

from db_adapter.constants import set_db_config_file_path
from db_adapter.constants import connection as con_params
from db_adapter.base import get_Pool, destroy_Pool
from db_adapter.curw_sim.constants import FLO2D_250, FLO2D_150, FLO2D_150_V2
from db_adapter.curw_sim.timeseries.discharge import Timeseries as DTimeseries
from db_adapter.curw_fcst.timeseries import Timeseries as Fcst_Timeseries
from db_adapter.curw_fcst.source import get_source_id


DATE_TIME_FORMAT = '%Y-%m-%d %H:%M:%S'
ROOT_DIR = '/home/curw/event_sim_utils'


def check_time_format(time, model):
    # hourly timeseries
    try:
        time = datetime.strptime(time, DATE_TIME_FORMAT)

        if time.strftime('%S') != '00':
            print("Seconds should be always 00")
            exit(1)
        # if model=="flo2d_250" and time.strftime('%M') not in ('05', '10', '15', '20', '25', '30', '35', '40', '45', '50', '55', '00'):
        #     print("Minutes should be multiple of 5 fro flo2d_250")
        #     exit(1)
        if model in ("flo2d_250", "flo2d_150", "flo2d_150_v2") and time.strftime('%M') != '00':
            print("Minutes should be always 00")
            exit(1)

        return True
    except Exception:
        traceback.print_exc()
        print("Time {} is not in proper format".format(time))
        exit(1)


def list_of_lists_to_df_first_row_as_columns(data):
    return pd.DataFrame.from_records(data[1:], columns=data[0])


def process_fcst_ts_from_hechms_outputs(curw_fcst_pool, extract_stations, i, start, end, sim_tag=None):

    FCST_TS = Fcst_Timeseries(curw_fcst_pool)

    try:
        # [station_name,latitude,longitude,target,model,version,sim_tag,station]
        source_model = extract_stations[i][4]
        version = extract_stations[i][5]
        station_id = extract_stations[i][7]

        if sim_tag is None:
            sim_tag = extract_stations[i][6]

        variable_id = 3 # Discharge
        unit_id = 3 # m3/s | Instantaneous

        source_id = get_source_id(pool=curw_fcst_pool, model=source_model, version=version)
        print(sim_tag, station_id, source_id, variable_id, unit_id)
        fcst_series = FCST_TS.get_latest_timeseries(sim_tag, station_id, source_id, variable_id, unit_id, start=None)
        print(fcst_series)
        if (fcst_series is None) or (len(fcst_series)<1):
            return None

        fcst_series.insert(0, ['time', 'value'])
        fcst_df = list_of_lists_to_df_first_row_as_columns(fcst_series)

        if start is None:
            start = (fcst_df['time'].min()).strftime(DATE_TIME_FORMAT)
        if end is None:
            end = (fcst_df['time'].max()).strftime(DATE_TIME_FORMAT)


        df = (pd.date_range(start=start, end=end, freq='60min')).to_frame(name='time')

        processed_df = pd.merge(df, fcst_df, on='time', how='left')

        processed_df.interpolate(method='linear', limit_direction='both', limit=100)
        processed_df.fillna(inplace=True, value=0)

        processed_df['time'] = processed_df['time'].dt.strftime(DATE_TIME_FORMAT)

        return processed_df.values.tolist()

    except Exception as e:
        traceback.print_exc()


def update_discharge_from_hechms(curw_sim_pool, curw_fcst_pool, flo2d_model, method, start_time, end_time, sim_tag):
    try:
        TS = DTimeseries(pool=curw_sim_pool)

        # [station_name,latitude,longitude,target,model,version,sim_tag,station]
        extract_stations = read_csv(os.path.join(ROOT_DIR, 'grids/discharge_stations/flo2d_stations.csv'))

        for i in range(len(extract_stations)):
            station_name = extract_stations[i][0]
            latitude = extract_stations[i][1]
            longitude = extract_stations[i][2]
            target_model = extract_stations[i][3]

            if target_model == flo2d_model:
                if station_name in ():
                    meta_data = {
                        'latitude': float('%.6f' % float(latitude)),
                        'longitude': float('%.6f' % float(longitude)),
                        'model': target_model, 'method': method,
                        'grid_id': 'discharge_{}'.format(station_name)
                    }

                    tms_id = TS.get_timeseries_id_if_exists(meta_data=meta_data)

                    if tms_id is None:
                        tms_id = TS.generate_timeseries_id(meta_data=meta_data)
                        meta_data['id'] = tms_id
                        TS.insert_run(meta_data=meta_data)

                    processed_discharge_ts = process_fcst_ts_from_hechms_outputs(curw_fcst_pool=curw_fcst_pool,
                                                                                 extract_stations= extract_stations,
                                                                                 i=i, start=start_time, end=end_time,
                                                                                 sim_tag=sim_tag)

                    if processed_discharge_ts is not None and len(processed_discharge_ts) > 0:
                        TS.insert_data(timeseries=processed_discharge_ts, tms_id=tms_id, upsert=True)

                else:
                    continue # skip the current iteration

    except Exception as e:
        traceback.print_exc()


def usage():
    usageText = """
    --------------------------------------------------
    Populate discharge Flo2D 250, 150 & 150_v2 :: OBS
    --------------------------------------------------

    Usage: ./discharge/flo2d_HECHMS.py [-m flo2d_XXX][-s "YYYY-MM-DD HH:MM:SS"] [-e "YYYY-MM-DD HH:MM:SS"]

    -h  --help          Show usage
    -m  --flo2d_model   FLO2D model (e.g. flo2d_250, flo2d_150). Default is flo2d_250.
    -s  --start_time    Discharge timeseries start time (e.g: "2019-06-05 00:00:00"). Default is 23:00:00, 3 days before today.
    -e  --end_time      Discharge timeseries end time (e.g: "2019-06-05 23:00:00"). Default is 23:00:00, tomorrow.
    -t  --hec_sim_tag   Simulation Tag of the desired HecHMS run from which you wanna prepare the Discharge timeseries.
                        (e.g: "event_run")
    """
    print(usageText)


if __name__=="__main__":

    set_db_config_file_path(os.path.join(ROOT_DIR, 'db_adapter_config.json'))

    try:
        start_time = None
        end_time = None
        flo2d_model = None
        method = "HECHMS"
        hechms_sim_tag = None

        try:
            opts, args = getopt.getopt(sys.argv[1:], "h:m:s:e:t:",
                                       ["help", "flo2d_model=", "start_time=", "end_time=", "hec_sim_tag="])
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
            elif opt in ("-t", "--hec_sim_tag"):
                hechms_sim_tag = arg.strip()

        if flo2d_model is None:
            print("Flo2d model is not specified.")
            exit(1)
        elif flo2d_model not in ("flo2d_250", "flo2d_150", "flo2d_150_v2"):
            print("Flo2d model should be either \"flo2d_250\" or \"flo2d_150\" or \"flo2d_150_v2\"")
            exit(1)

        if start_time is None:
            start_time = (datetime.now() - timedelta(days=3)).strftime('%Y-%m-%d 23:00:00')
        else:
            check_time_format(time=start_time, model=flo2d_model)

        if end_time is None:
            end_time = (datetime.now() + timedelta(days=1)).strftime('%Y-%m-%d 23:00:00')
        else:
            check_time_format(time=end_time, model=flo2d_model)

        timestep = 60

        curw_sim_pool = get_Pool(host=con_params.CURW_SIM_HOST, user=con_params.CURW_SIM_USERNAME,
                                 password=con_params.CURW_SIM_PASSWORD, port=con_params.CURW_SIM_PORT,
                                 db=con_params.CURW_SIM_DATABASE)

        curw_fcst_pool = get_Pool(host=con_params.CURW_FCST_HOST, user=con_params.CURW_FCST_USERNAME,
                                 password=con_params.CURW_FCST_PASSWORD, port=con_params.CURW_FCST_PORT,
                                 db=con_params.CURW_FCST_DATABASE)


        print("{} : ####### Insert hechms discharge series for {}.".format(datetime.now(), flo2d_model))
        update_discharge_from_hechms(curw_sim_pool=curw_sim_pool, curw_fcst_pool=curw_fcst_pool, flo2d_model=flo2d_model,
                                     method=method, start_time=start_time, end_time=end_time, sim_tag=hechms_sim_tag)

    except Exception as e:
        traceback.print_exc()
    finally:
        destroy_Pool(pool=curw_sim_pool)
        print("Process finished.")


