import csv
import operator
import collections
import traceback

from math import acos, cos, sin, radians

from db_adapter.base import get_Pool, destroy_Pool
from db_adapter.curw_sim.grids import add_obs_to_d03_grid_mappings_for_rainfall, \
    get_obs_to_d03_grid_mappings_for_rainfall, \
    GridInterpolationEnum
from db_adapter.constants import CURW_SIM_HOST, CURW_SIM_PORT, CURW_SIM_USERNAME, CURW_SIM_PASSWORD, CURW_SIM_DATABASE


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


def find_nearest_d03_station_for_obs_grids(obs_stations_csv, d03_stations_csv):

    obs_grids = read_csv(obs_stations_csv)  # [hash_id,station_id,station_name,latitude,longitude]

    d03_stations = read_csv(d03_stations_csv)  # [id,latitude,longitude]

    nearest_d03_stations_list = [['obs_grid_id', 'd03_1_id', 'd03_1_dist', 'd03_2_id', 'd03_2_dist', 'd03_3_id', 'd03_3_dist']]

    for origin_index in range(len(obs_grids)):

        nearest_d03_station = [obs_grids[origin_index][1]]

        origin_lat = float(obs_grids[origin_index][3])
        origin_lng = float(obs_grids[origin_index][4])

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
            if count < 3 and sorted_distances.get(key) < 15:
                nearest_d03_station.extend([key, sorted_distances.get(key)])
                count += 1

        print(nearest_d03_station)
        nearest_d03_stations_list.append(nearest_d03_station)

    create_csv('grid_maps/obs_stations/rainfall/MDPA_obs_d03_stations_mapping.csv', nearest_d03_stations_list)


try:

    find_nearest_d03_station_for_obs_grids('grids/obs_stations/rainfall/curw_active_rainfall_obs_stations.csv',
                                           'grids/wrf_d03/d03_stations.csv')

    print(" Add obs to wrf_d03 grid mappings")

    pool = get_Pool(host=CURW_SIM_HOST, port=CURW_SIM_PORT, user=CURW_SIM_USERNAME, password=CURW_SIM_PASSWORD, db=CURW_SIM_DATABASE)

    grid_interpolation_method = GridInterpolationEnum.getAbbreviation(GridInterpolationEnum.MDPA)

    print(" Add obs to wrf_d03 grid mappings")
    add_obs_to_d03_grid_mappings_for_rainfall(pool=pool, grid_interpolation=grid_interpolation_method,
                                              obs_to_d03_map_path='grid_maps/obs_stations/rainfall/{}_obs_d03_stations_mapping.csv'
                                              .format(grid_interpolation_method),
                                              active_obs_path='grids/obs_stations/rainfall/curw_active_rainfall_obs_stations.csv')
    print("{} rainfall observed station grids added".format(len(get_obs_to_d03_grid_mappings_for_rainfall(pool=pool, grid_interpolation=grid_interpolation_method).keys())))

except Exception as e:
    traceback.print_exc()
finally:
    destroy_Pool(pool=pool)
    print("Process Finished.")


