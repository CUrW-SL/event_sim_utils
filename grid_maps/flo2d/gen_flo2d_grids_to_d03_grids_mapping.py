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


def find_nearest_d03_station_for_flo2d_grids(flo2d_stations_csv, d03_stations_csv, flo2d_model):

    flo2d_grids = read_csv(flo2d_stations_csv)

    d03_stations = read_csv(d03_stations_csv)

    nearest_d03_stations_list = [['flo2d_grid_id', 'nearest_d03_station_id', 'dist']]

    for origin_index in range(len(flo2d_grids)):

        nearest_d03_station = [flo2d_grids[origin_index][0]]

        origin_lat = float(flo2d_grids[origin_index][2])
        origin_lng = float(flo2d_grids[origin_index][1])

        distances = {}

        for d03_index in range(len(d03_stations)):
            lat = float(d03_stations[d03_index][1])
            lng = float(d03_stations[d03_index][2])

            intermediate_value = cos(radians(origin_lat)) * cos(radians(lat)) * cos(
                    radians(lng) - radians(origin_lng)) + sin(radians(origin_lat)) * sin(radians(lat))
            if intermediate_value < 1:
                distance = 6371 * acos(intermediate_value)
            else:
                distance = 6371 * acos(1)

            distances[d03_stations[d03_index][0]] = distance

        sorted_distances = collections.OrderedDict(sorted(distances.items(), key=operator.itemgetter(1))[:10])

        count = 0
        for key in sorted_distances.keys():
            if count < 1:
                nearest_d03_station.extend([key, sorted_distances.get(key)])
                count += 1

        print(nearest_d03_station)
        nearest_d03_stations_list.append(nearest_d03_station)

    create_csv('MDPA_{}_d03_stations_mapping.csv'.format(flo2d_model), nearest_d03_stations_list)


try:

    pool = get_Pool(host=CURW_SIM_HOST, port=CURW_SIM_PORT, user=CURW_SIM_USERNAME, password=CURW_SIM_PASSWORD,
                    db=CURW_SIM_DATABASE)

    grid_interpolation_method = GridInterpolationEnum.getAbbreviation(GridInterpolationEnum.MDPA)

    for flo2d_model in flo2d_models_list:
        print("Update {} grid mappings".format(flo2d_model))

        find_nearest_d03_station_for_flo2d_grids(flo2d_stations_csv='grids/flo2d/{}m.csv'.format(flo2d_model),
                                                 d03_stations_csv='grids/wrf_d03/d03_stations.csv',
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