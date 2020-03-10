import csv
import operator
import collections
import traceback

from math import acos, cos, sin, radians

from db_adapter.base import get_Pool, destroy_Pool
from db_adapter.constants import CURW_SIM_HOST, CURW_SIM_PORT, CURW_SIM_USERNAME, CURW_SIM_PASSWORD, CURW_SIM_DATABASE
from db_adapter.curw_sim.grids import add_flo2d_raincell_grid_mappings, \
    get_flo2d_cells_to_obs_grid_mappings, get_flo2d_cells_to_wrf_grid_mappings, \
    GridInterpolationEnum
from db_adapter.curw_sim.constants import FLO2D_250, FLO2D_150, FLO2D_30, FLO2D_150_V2

flo2d_models_list = [FLO2D_250, FLO2D_150, FLO2D_150_V2]


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


def read_csv(file_name):
    """
    Read csv file
    :param file_name: <file_path/file_name>.csv
    :return: list of lists which contains each row of the csv file
    """

    with open(file_name, 'r') as f:
        data = [list(line) for line in csv.reader(f)][1:]

    return data


def find_nearest_obs_stations_for_flo2d_stations(flo2d_stations_csv, obs_stations_csv, flo2d_model):

    obs_stations = read_csv(obs_stations_csv)  # [hash_id,station_id,station_name,latitude,longitude]

    flo2d_station = read_csv(flo2d_stations_csv)  # [Grid_ID,X,Y]

    flo2d_obs_mapping_list = [['{}_station_id'.format(flo2d_model), 'ob_1_id', 'ob_1_dist', 'ob_2_id', 'ob_2_dist', 'ob_3_id',
                               'ob_3_dist']]

    for flo2d_index in range(len(flo2d_station)):

        flo2d_obs_mapping = [flo2d_station[flo2d_index][0]]

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
                flo2d_obs_mapping.extend([key, sorted_distances.get(key)])
                count += 1
            elif count < 3:
                flo2d_obs_mapping.extend([-1, -1])
                count += 1

        # print(flo2d_obs_mapping)
        flo2d_obs_mapping_list.append(flo2d_obs_mapping)

    create_csv('grid_maps/flo2d/MDPA_{}_obs_mapping.csv'.format(flo2d_model), flo2d_obs_mapping_list)

# find_nearest_obs_stations_for_flo2d_stations('flo2d_30m.csv', 'curw_active_rainfall_obs_stations.csv')


try:

    pool = get_Pool(host=CURW_SIM_HOST, port=CURW_SIM_PORT, user=CURW_SIM_USERNAME, password=CURW_SIM_PASSWORD,
                    db=CURW_SIM_DATABASE)

    grid_interpolation_method = GridInterpolationEnum.getAbbreviation(GridInterpolationEnum.MDPA)

    for flo2d_model in flo2d_models_list:
        print("Update {} grid mappings".format(flo2d_model))

        find_nearest_obs_stations_for_flo2d_stations(flo2d_stations_csv='grids/flo2d/{}m.csv'.format(flo2d_model),
                                                     obs_stations_csv='grids/obs_stations/rainfall/curw_active_rainfall_obs_stations.csv',
                                                     flo2d_model=flo2d_model)

        add_flo2d_raincell_grid_mappings(pool=pool, flo2d_model=flo2d_model, grid_interpolation=grid_interpolation_method,
                                         obs_map_file_path='grid_maps/flo2d/{}_{}_obs_mapping.csv'
                                         .format(grid_interpolation_method, flo2d_model),
                                         d03_map_file_path='grid_maps/flo2d/{}_{}_d03_stations_mapping.csv'
                                         .format(grid_interpolation_method, flo2d_model))
        print("{} {} grids added".format(len(get_flo2d_cells_to_wrf_grid_mappings(
            pool=pool, flo2d_model=flo2d_model, grid_interpolation=grid_interpolation_method).keys()), flo2d_model))
        print("{} {} grids added".format(len(get_flo2d_cells_to_obs_grid_mappings(
            pool=pool, flo2d_model=flo2d_model, grid_interpolation=grid_interpolation_method).keys()), flo2d_model))


except Exception as e:
    traceback.print_exc()
finally:
    destroy_Pool(pool=pool)
    print("Process Finished.")