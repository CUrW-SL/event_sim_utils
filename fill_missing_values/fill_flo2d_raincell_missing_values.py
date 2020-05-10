#!/home/curw/event_sim_utils/venv/bin/python3
import pymysql
import traceback
import getopt
import sys, os
from datetime import datetime, timedelta

from db_adapter.constants import set_db_config_file_path
from db_adapter.constants import connection as con_params

DATE_TIME_FORMAT = '%Y-%m-%d %H:%M:%S'
ROOT_DIR = '/home/curw/event_sim_utils'

# connection params
HOST = ""
USER = ""
PASSWORD = ""
DB ="curw_sim"
PORT = 3306

FLO2D_250 = "flo2d_250"
FLO2D_150 = "flo2d_150"
FLO2D_30 = "flo2d_30"
HecHMS = "hechms"

OPTIONS = ['OBS', 'FCST', 'CHECK']


def check_for_missing_values(start, end, model, method):

    # Connect to the database
    connection = pymysql.connect(host=HOST, user=USER, password=PASSWORD, db=DB,
            cursorclass=pymysql.cursors.DictCursor)

    print("Connected to database")

    if model == FLO2D_250:
        timestep = 5
    elif model == FLO2D_150:
        timestep = 15
    elif model == HecHMS:
        timestep = 5
    elif model == FLO2D_30:
        timestep = 15

    start = datetime.strptime(start, DATE_TIME_FORMAT)
    end = datetime.strptime(end, DATE_TIME_FORMAT)

    expected_count = int(((end-start).total_seconds()/60)/timestep)

    ids = []

    try:
        with connection.cursor() as cursor1:
            sql_statement = "select `id` from `run` where `model`=%s and `method`=%s;" # " \
                            # "and (`obs_end` is null or `obs_end` < %s);"
            cursor1.execute(sql_statement, (model, method))
            # cursor1.execute(sql_statement, (model, expected_obs_end))
            results = cursor1.fetchall()
            for result in results:
                ids.append([result.get('id')])

        for id in ids:
            with connection.cursor() as cursor2:
                sql_statement = "select count(`time`) as `count` from `data` where id=%s and `time`>%s and `time`<=%s;"
                cursor2.execute(sql_statement, (id, start, end))
                count = cursor2.fetchone()['count']

            if count < expected_count:
                print(id)
                print("{} timesteps are missing".format(str(expected_count-count)))


    except Exception as ex:
        traceback.print_exc()
    finally:
        connection.close()
        print("{} process completed".format(datetime.now()))


def fill_missing_obs_with_0s(start, end, model, method):

    # Connect to the database
    connection = pymysql.connect(host=HOST, user=USER, password=PASSWORD, db=DB,
            cursorclass=pymysql.cursors.DictCursor)

    print("Connected to database")

    if model == FLO2D_250:
        timestep = 5
    elif model == FLO2D_150:
        timestep = 15
    elif model == HecHMS:
        timestep = 5
    elif model == FLO2D_30:
        timestep = 15

    start = datetime.strptime(start, DATE_TIME_FORMAT)
    end = datetime.strptime(end, DATE_TIME_FORMAT)

    expected_count = int(((end-start).total_seconds()/60)/timestep)

    ids = []

    try:

        with connection.cursor() as cursor1:
            sql_statement = "select `id` from `run` where `model`=%s and `method`=%s;"  # " \
            # "and (`obs_end` is null or `obs_end` < %s);"
            cursor1.execute(sql_statement, (model, method))
            results = cursor1.fetchall()
            for result in results:
                ids.append([result.get('id')])

        for id in ids:
            with connection.cursor() as cursor2:
                sql_statement = "select count(`time`) as `count` from `data` where id=%s and `time`>%s and `time`<=%s;"
                cursor2.execute(sql_statement, (id, start, end))
                count = cursor2.fetchone()['count']

            if count < expected_count:
                print(id)
                timestamp = start
                while timestamp <= end:
                    print(timestamp)

                    try:
                        with connection.cursor() as cursor3:
                            sql_statement = "INSERT INTO `data` (`id`,`time`,`value`) VALUES (%s,%s,%s);"
                            cursor3.execute(sql_statement, (id, timestamp, 0))
                        connection.commit()
                    except Exception as ex:
                        connection.rollback()
                        traceback.print_exc()

                    timestamp = timestamp + timedelta(minutes=timestep)

    except Exception as ex:
        traceback.print_exc()
    finally:
        connection.close()
        print("{} process completed".format(datetime.now()))


def fill_missing_fcsts(end, model, method):

    #Connect to the database
    connection = pymysql.connect(host=HOST, user=USER, password=PASSWORD, db=DB,
            cursorclass=pymysql.cursors.DictCursor)

    print("Connected to database")

    if model == FLO2D_250:
        timestep = 5
    elif model == FLO2D_150:
        timestep = 15
    elif model == HecHMS:
        timestep = 5
    elif model == FLO2D_30:
        timestep = 15

    end = datetime.strptime(end, DATE_TIME_FORMAT)

    ids = []

    try:

        with connection.cursor() as cursor1:
            sql_statement = "select `id` from `run` where `model`=%s and `method`=%s;"
            cursor1.execute(sql_statement, (model, method))
            results = cursor1.fetchall()
            for result in results:
                ids.append([result.get('id')])

        for id in ids:
            with connection.cursor() as cursor2:
                sql_statement = "select max(`time`) as `time` from `data` where id=%s;"
                cursor2.execute(sql_statement, id)
                start = cursor2.fetchone()['time']

            print(id, start, end)

            ts0 = []

            timestamp = start
            while timestamp < end:
                timestamp = timestamp + timedelta(minutes=timestep)
                ts0.append([id, timestamp, 0])

            try:
                with connection.cursor() as cursor3:
                    sql_statement = "INSERT INTO `data` (`id`,`time`,`value`) VALUES (%s,%s,%s);"
                    cursor3.executemany(sql_statement, ts0)
                connection.commit()
            except Exception as ex:
                connection.rollback()
                traceback.print_exc()


    except Exception as ex:
        traceback.print_exc()
    finally:
        connection.close()
        print("{} process completed".format(datetime.now()))


def usage():
    usageText = """
    Usage: ./fill_missing_values/fill_flo2d_raincell_missing_values.py [-m XXXX][-s "YYYY-MM-DD HH:MM:SS"] [-e "YYYY-MM-DD HH:MM:SS"] [-o [OBS|FCST]] [-E]

    -h  --help          Show usage
    -m  --model         Model (e.g. flo2d_250, flo2d_150, hechms). Default is flo2d_250.
    -M  --method        Timeseries calculation method. (e.g. OBS, MME) Default is MME.
    -s  --start_time    First timestamp to  be filled with dummy values if value is missing (e.g: "2019-06-05 00:00:00"). 
                        Default is 00:00:00, yesterday.
    -e  --end_time      Last timestamp to be filled with dummy values if value is missing (e.g: "2019-06-05 23:30:00"). 
                        Default is 23:30:00, tomorrow.
    -o  --option        OBS | FCST. If OBS, fill gaps in between start_time and end_time. If FCST, fill fcsts until end_time.
                        Default is FCST.
    -E  --event_sim     Considered databse, event database or not (e.g. -E, --event_sim)
    """
    print(usageText)


def check_time_format(time, model):
    try:
        time = datetime.strptime(time, DATE_TIME_FORMAT)

        if time.strftime('%S') != '00':
            print("Seconds should be always 00")
            exit(1)
        if (model==FLO2D_250 or model==HecHMS) and time.strftime('%M') not in ('05', '10', '15', '20', '25', '30', '35', '40', '45', '50', '55', '00'):
            print("Minutes should be multiple of 5 fro flo2d_250 and hechms")
            exit(1)
        if (model==FLO2D_150 or model==FLO2D_30) and time.strftime('%M') not in ('15', '30', '45', '00'):
            print("Minutes should be multiple of 15 for flo2d_150 and flo2d_30")
            exit(1)

        return True
    except Exception:
        traceback.print_exc()
        print("Time {} is not in proper format".format(time))
        exit(1)


if __name__=="__main__":


    set_db_config_file_path('/home/curw/event_sim_utils/db_adapter_config.json')

    HOST = con_params.CURW_SIM_HOST
    USER = con_params.CURW_SIM_USERNAME
    PASSWORD = con_params.CURW_SIM_PASSWORD
    DB = con_params.CURW_SIM_DATABASE
    PORT = con_params.CURW_SIM_PORT

    try:
        start_time = None
        end_time = None
        model = None
        method=None
        option = None
        event_sim = False

        try:
            opts, args = getopt.getopt(sys.argv[1:], "h:m:M:s:e:o:E",
                    ["help", "model=", "method=", "start_time=", "end_time=", "option=", "event_sim"])
        except getopt.GetoptError:
            usage()
            sys.exit(2)
        for opt, arg in opts:
            if opt in ("-h", "--help"):
                usage()
                sys.exit()
            elif opt in ("-m", "--model"):
                model = arg.strip()
            elif opt in ("-M", "--method"):
                method = arg.strip()
            elif opt in ("-s", "--start_time"):
                start_time = arg.strip()
            elif opt in ("-e", "--end_time"):
                end_time = arg.strip()
            elif opt in ("-o", "--option"):
                option = arg.strip()
            elif opt in ("-E", "--event_sim"):
                event_sim = True

        if event_sim:
            set_db_config_file_path(os.path.join(ROOT_DIR, 'db_adapter_config_event_sim.json'))

        print(model, start_time, end_time, option)

        if model is None:
            model = FLO2D_250
        elif model not in (FLO2D_250, FLO2D_150, HecHMS, FLO2D_30):
            print("Model should be one of theses; \"flo2d_250\", \"flo2d_150\", \"hechms\", \"flo2d_30\". ")
            exit(1)

        if end_time is None:
            end_time = (datetime.now() + timedelta(days=3)).strftime('%Y-%m-%d 00:00:00')
        else:
            check_time_format(time=end_time, model=model)

        if method is None:
            method="MME"

        if option is None:
            option = OPTIONS[1]


        ################################################
        # Fill missing fcsts until the given end_time. #
        ################################################
        if option==OPTIONS[1]:
            fill_missing_fcsts(end=end_time, model=model)

        ################################################
        # Fill gaps in between start_time and end_time #
        ################################################
        elif option== OPTIONS[0]:
            if start_time is None:
                start_time = (datetime.now() - timedelta(days=1)).strftime('%Y-%m-%d 00:00:00')
            else:
                check_time_format(time=start_time, model=model)

            fill_missing_obs_with_0s(start=start_time, end=end_time, model=model)

        ################################################
        # Check whether  #
        ################################################
        elif option == OPTIONS[2]:
            if start_time is None:
                print("You haven't specified the start time")
                exit(1)

            check_for_missing_values(start=start_time, end=end_time, model=model)

    except Exception:
        traceback.print_exc()
