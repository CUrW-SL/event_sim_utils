#!/home/curw/event_sim_utils/venv/bin/python3

import sys
import getopt
import os
import traceback
from datetime import datetime, timedelta

from db_adapter.csv_utils import read_csv

from db_adapter.constants import set_db_config_file_path
from db_adapter.constants import connection as con_params
from db_adapter.base import get_Pool, destroy_Pool
from db_adapter.curw_sim.timeseries.waterlevel import Timeseries
from db_adapter.curw_sim.common import fill_ts_missing_entries
from db_adapter.curw_sim.constants import FLO2D_250, FLO2D_150, FLO2D_150_V2

DATE_TIME_FORMAT = '%Y-%m-%d %H:%M:%S'
ROOT_DIR = '/home/curw/event_sim_utils'


RANWALA_WL_ID = 100064


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


def calculate_hanwella_wl_from_ranwala(ranwala_ts):
    hanwella_ts = []
    "Han[i]=1.623194Ran[i]-5.16108(Ran[i]-Ran[i-1])-10.847356"

    for i in range(len(ranwala_ts) - 1):
        # x = Ranwala
        # DX Ranwala (x[x] - x[t-1]}
        # Hanwella = X 1.642174610188251` -   DX 3.8585516925010444` -
        #    8.810870547723741`;
        x = float(ranwala_ts[i + 1][1])
        dx = float(ranwala_ts[i + 1][1] - ranwala_ts[i][1])
        hanwella_wl = x * 1.623194 - dx * 5.16108 - 10.847356
        hanwella_ts.append([ranwala_ts[i + 1][0], '%.3f' % hanwella_wl])

    for i in range(len(hanwella_ts)):
        if float(hanwella_ts[i][1]) < 0.2:
            hanwella_ts[i][1] = 0.2

    return hanwella_ts


def calculate_glencourse_wl_from_ranwala(ranwala_ts):
    "-5.908532+2.2784865x-0.0309476x2"
    glencourse_ts = []

    for i in range(len(ranwala_ts)):
        ranwala_wl = float(ranwala_ts[i][1])
        glencourse_wl = -5.908532 + (2.2784865 * ranwala_wl) - (0.0309476 * (ranwala_wl ** 2))
        glencourse_ts.append([ranwala_ts[i][0], '%.3f' % glencourse_wl])

    return glencourse_ts


def update_waterlevel_obs(obs_connection, curw_sim_pool, flo2d_model, method, timestep, start_time, end_time):
    try:

        # [station_name,latitude,longitude,target]
        extract_stations = read_csv('grids/waterlevel_stations/extract_stations.csv')
        extract_stations_dict = {}  # keys: target_model , value: [latitude, longitude, station_name]
        # older version ::: keys: station_name , value: [latitude, longitude, target_model]

        for obs_index in range(len(extract_stations)):
            extract_stations_dict[extract_stations[obs_index][3]] = [extract_stations[obs_index][1],
                                                                     extract_stations[obs_index][2],
                                                                     extract_stations[obs_index][0]]

        station_name = extract_stations_dict.get(flo2d_model)[2]
        meta_data = {
            'latitude': float('%.6f' % float(extract_stations_dict.get(flo2d_model)[0])),
            'longitude': float('%.6f' % float(extract_stations_dict.get(flo2d_model)[1])),
            'model': flo2d_model, 'method': method,
            'grid_id': 'waterlevel_{}'.format(station_name)
        }

        TS = Timeseries(pool=curw_sim_pool)

        tms_id = TS.get_timeseries_id_if_exists(meta_data=meta_data)

        ranwala_ts = []

        if tms_id is None:
            tms_id = TS.generate_timeseries_id(meta_data=meta_data)
            meta_data['id'] = tms_id
            TS.insert_run(meta_data=meta_data)

        with obs_connection.cursor() as cursor1:
            cursor1.callproc('getWL', (RANWALA_WL_ID, start_time, end_time))
            results = cursor1.fetchall()
            for result in results:
                ranwala_ts.append([result.get('time'), result.get('value')])

        interpolated_ranwala_ts = fill_ts_missing_entries(start=start_time, end=end_time, timeseries=ranwala_ts,
                                                          interpolation_method='linear', timestep=60)

        estimated_wl_ts = []

        if station_name == 'hanwella':
            estimated_wl_ts = calculate_hanwella_wl_from_ranwala(interpolated_ranwala_ts)
        elif station_name == 'glencourse':
            estimated_wl_ts = calculate_glencourse_wl_from_ranwala(interpolated_ranwala_ts)

        if estimated_wl_ts is not None and len(estimated_wl_ts) > 0:
            TS.insert_data(timeseries=estimated_wl_ts, tms_id=tms_id, upsert=True)

    except Exception as e:
        traceback.print_exc()


def usage():
    usageText = """
    --------------------------------------------------
    Populate waterlevel Flo2D 250, 150 & 150_v2 :: OBS
    --------------------------------------------------

    Usage: ./waterlevel/flo2d_OBS.py [-m flo2d_XXX][-s "YYYY-MM-DD HH:MM:SS"] [-e "YYYY-MM-DD HH:MM:SS"]

    -h  --help          Show usage
    -m  --flo2d_model   FLO2D model (e.g. flo2d_250, flo2d_150). Default is flo2d_250.
    -s  --start_time    Waterlevel timeseries start time (e.g: "2019-06-05 00:00:00"). Default is 23:00:00, 3 days before today.
    -e  --end_time      Waterlevel timeseries end time (e.g: "2019-06-05 23:00:00"). Default is 23:00:00, tomorrow.
    """
    print(usageText)


if __name__=="__main__":

    set_db_config_file_path(os.path.join(ROOT_DIR, 'db_adapter_config.json'))

    try:
        start_time = None
        end_time = None
        flo2d_model = None
        method = "OBS"

        try:
            opts, args = getopt.getopt(sys.argv[1:], "h:m:s:e:",
                                       ["help", "flo2d_model=", "start_time=", "end_time="])
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

        if flo2d_model == FLO2D_250:
            timestep = 60
        elif flo2d_model == FLO2D_150:
            timestep = 60
        elif flo2d_model == FLO2D_150_V2:
            timestep = 60


        curw_sim_pool = get_Pool(host=con_params.CURW_SIM_HOST, user=con_params.CURW_SIM_USERNAME,
                                 password=con_params.CURW_SIM_PASSWORD, port=con_params.CURW_SIM_PORT,
                                 db=con_params.CURW_SIM_DATABASE)

        curw_obs_pool = get_Pool(host=con_params.CURW_OBS_HOST, user=con_params.CURW_OBS_USERNAME,
                                 password=con_params.CURW_OBS_PASSWORD, port=con_params.CURW_OBS_PORT,
                                 db=con_params.CURW_OBS_DATABASE)

        obs_connection = curw_obs_pool.connection()

        print("{} : ####### Insert obs waterlevel series for {}.".format(datetime.now(), flo2d_model))
        update_waterlevel_obs(obs_connection=obs_connection, curw_sim_pool=curw_sim_pool, flo2d_model=flo2d_model,
                              method=method, timestep=timestep, start_time=start_time, end_time=end_time)

    except Exception as e:
        traceback.print_exc()
    finally:
        obs_connection.close()
        destroy_Pool(pool=curw_obs_pool)
        destroy_Pool(pool=curw_sim_pool)
        print("Process finished.")


