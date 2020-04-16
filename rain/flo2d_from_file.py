#!/home/curw/event_sim_utils/venv/bin/python3

import sys
import getopt
import os
import json
import traceback
from datetime import datetime, timedelta

import operator
import collections
from math import acos, cos, sin, radians

import pandas as pd
import numpy as np
import geopandas as gpd
from scipy.spatial import Voronoi
from shapely.geometry import Polygon, Point

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


def find_nearest_stations_for_flo2d_grids(flo2d_grids_csv, stations_dict):
    """

    :param flo2d_stations_csv:
    :param stations_dict: key = point_name, value = [latitude,longitude]
    :param flo2d_model: 
    :return:
    """

    flo2d_station = read_csv(flo2d_grids_csv)  # [Grid_ID,X,Y]

    # key: flo2d_grid_id   value: [lat, lon]
    flo2d_stations_mapping_dict = {}

    for flo2d_index in range(len(flo2d_station)):

        flo2d_lat = float(flo2d_station[flo2d_index][2])
        flo2d_lng = float(flo2d_station[flo2d_index][1])

        distances = {}

        for key in stations_dict.keys():
            lat = float(stations_dict.get(key)[1])
            lng = float(stations_dict.get(key)[0])

            intermediate_value = cos(radians(flo2d_lat)) * cos(radians(lat)) * cos(
                    radians(lng) - radians(flo2d_lng)) + sin(radians(flo2d_lat)) * sin(radians(lat))
            if intermediate_value < 1:
                distance = 6371 * acos(intermediate_value)
            else:
                distance = 6371 * acos(1)

            distances[key] = distance

        sorted_distances = collections.OrderedDict(sorted(distances.items(), key=operator.itemgetter(1))[:10])
        flo2d_stations_mapping = []

        count = 0
        for key in sorted_distances.keys():
            if count < 1 and sorted_distances.get(key) <= 25:
                flo2d_stations_mapping = stations_dict.get(key)
                count += 1
            elif count < 1:
                flo2d_stations_mapping = [-1, -1]
                count += 1
            else:
                continue

        # print(flo2d_obs_mapping)
        flo2d_stations_mapping_dict[flo2d_station[flo2d_index][0]] = flo2d_stations_mapping

    print(json.dumps(flo2d_stations_mapping_dict))
    return flo2d_stations_mapping_dict


def _voronoi_finite_polygons_2d(vor, radius=None):
    """
    Reconstruct infinite voronoi regions in a 2D diagram to finite
    regions.
    Parameters
    ----------
    vor : Voronoi
        Input diagram
    radius : float, optional
        Distance to 'points at infinity'.
    Returns
    -------
    regions : list of tuples
        Indices of vertices in each revised Voronoi regions.
    vertices : list of tuples
        Coordinates for revised Voronoi vertices. Same as coordinates
        of input vertices, with 'points at infinity' appended to the
        end.
    from: https://stackoverflow.com/questions/20515554/colorize-voronoi-diagram
    """
    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()
    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max()
    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]
        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue
        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue
            # Compute the missing endpoint of an infinite ridge
            t = vor.points[p2] - vor.points[p1]  # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius
            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:, 1] - c[1], vs[:, 0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]
        # finish
        new_regions.append(new_region.tolist())
    return new_regions, np.asarray(new_vertices)


def get_voronoi_polygons(points_dict, shape_file, shape_attribute=None, output_shape_file=None, add_total_area=True):
    """
    :param points_dict: dict of points {'id' --> [lon, lat]}
    :param shape_file: shape file path of the area
    :param shape_attribute: attribute list of the interested region [key, value]
    :param output_shape_file: if not none, a shape file will be created with the output
    :param add_total_area: if true, total area shape will also be added to output
    :return:
    geo_dataframe with voronoi polygons with columns ['id', 'lon', 'lat','area', 'geometry'] with last row being the area of the
    shape file
    """
    if shape_attribute is None:
        shape_attribute = ['OBJECTID', 1]

    shape_df = gpd.GeoDataFrame.from_file(shape_file)
    shape_polygon_idx = shape_df.index[shape_df[shape_attribute[0]] == shape_attribute[1]][0]
    shape_polygon = shape_df['geometry'][shape_polygon_idx]

    ids = [p if type(p) == str else np.asscalar(p) for p in points_dict.keys()]
    points = np.array(list(points_dict.values()))[:, :2]
    vor = Voronoi(points)

    regions, vertices = _voronoi_finite_polygons_2d(vor)

    data = []
    for i, region in enumerate(regions):
        polygon = Polygon([tuple(x) for x in vertices[region]])
        if polygon.intersects(shape_polygon):
            intersection = polygon.intersection(shape_polygon)
            data.append({'id': ids[i], 'lon': vor.points[i][0], 'lat': vor.points[i][1], 'area': intersection.area,
                         'geometry': intersection
                         })
    df = gpd.GeoDataFrame(data, columns=['id', 'lon', 'lat', 'area', 'geometry'], crs=shape_df.crs)
    if output_shape_file is not None:
        df.to_file(output_shape_file)

    return df


def divide_flo2d_grids_to_polygons(flo2d_model, polygons):

    flo2d_grids = read_csv(
        os.path.join(ROOT_DIR, 'grids/flo2d/{}m.csv'.format(flo2d_model)))  # [Grid_ ID, X(longitude), Y(latitude)]

    for grid in flo2d_grids:
        point = Point(float(grid[1]), float(grid[2]))

        for index, row in polygons.iterrows():
            polygon = polygons.iloc[index]['geometry']
            if point.within(polygon):
                grid.append(polygons.iloc[index]['id'])
                continue

    return flo2d_grids


# for bulk insertion for a given one grid interpolation method
def update_rainfall_from_file(flo2d_grid_polygon_map, stations_dict, rainfall_df, mean_rf, flo2d_model, method,
                              grid_interpolation, timestep, start_time=None, end_time=None):

    """
    Update rainfall observations for flo2d models
    :param flo2d_model: flo2d model
    :param method: value interpolation method
    :param grid_interpolation: grid interpolation method
    :param timestep: output timeseries timestep
    :return:
    """

    # start = datetime.strptime(start_time, '%Y-%m-%d %H:%M:%S')

    try:

        # Connect to the database
        curw_sim_pool = get_Pool(host=con_params.CURW_SIM_HOST, user=con_params.CURW_SIM_USERNAME,
                                 password=con_params.CURW_SIM_PASSWORD,
                                 port=con_params.CURW_SIM_PORT, db=con_params.CURW_SIM_DATABASE)

        TS = Sim_Timeseries(pool=curw_sim_pool)

        # # [hash_id, station_id, station_name, latitude, longitude]
        # flo2d_grid_polygon_map :: [Grid_ ID, X(longitude), Y(latitude), matching_point]

        # stations_dict_for_obs = { }  # keys: obs station id , value: hash id

        for grid in flo2d_grid_polygon_map:
            lat = grid[2]
            lon = grid[1]
            cell_id = grid[0]
            meta_data = {
                    'latitude': float('%.6f' % float(lat)), 'longitude': float('%.6f' % float(lon)),
                    'model': flo2d_model, 'method': method,
                    'grid_id': '{}_{}_{}'.format(flo2d_model, grid_interpolation, (str(cell_id)).zfill(10))
                    }

            if len(grid) > 3:
                polygon = grid[3]

                poly_lat = stations_dict.get(polygon)[1]
                poly_lon = stations_dict.get(polygon)[0]

                timeseries = rainfall_df.loc[
                    (rainfall_df['latitude'] == poly_lat) & (rainfall_df['longitude'] == poly_lon)].values.tolist()

            else:
                timeseries = mean_rf

            tms_id = TS.get_timeseries_id(grid_id=meta_data.get('grid_id'), method=meta_data.get('method'))

            if tms_id is None:
                tms_id = TS.generate_timeseries_id(meta_data=meta_data)
                meta_data['id'] = tms_id
                TS.insert_run(meta_data=meta_data)

            print("grid_id:", meta_data['grid_id'])

            # for i in range(len(obs_timeseries)):
            #     if obs_timeseries[i][1] == -99999:
            #         obs_timeseries[i][1] = 0

            if timeseries is not None and len(timeseries) > 0:
                TS.insert_data(timeseries=timeseries, tms_id=tms_id, upsert=True)

    except Exception as e:
        traceback.print_exc()
        logger.error("Exception occurred while updating obs rainfalls in curw_sim.")
    finally:
        destroy_Pool(pool=curw_sim_pool)
        logger.info("Process finished")


def usage():
    usageText = """
    -------------------------------------------------
    Populate rainfall Flo2D 250, 150 & 150_v2 :: FF
    -------------------------------------------------

    Usage: ./rain/flo2d_from_file.py [-m flo2d_XXX] [-M XXX] [-g XXXX] [-f <file_path>] 
    [-s "YYYY-MM-DD HH:MM:SS"] [-e "YYYY-MM-DD HH:MM:SS"]

    -h  --help          Show usage
    -m  --flo2d_model   FLO2D model (e.g. flo2d_250, flo2d_150). Default is flo2d_250.
    -M  --method        Rain interpolation method ("BC"-Bias Corrected, "MME"-Multi Model Ensemble, "FF"-From File). Default is FF.
    -g  --grid_tag      Grid mapping method (e.g: "MDPA", "TP"). Default is "MDPA".
    -f  --file          Path to the CSV file containing the rainfall data (e.g: "/home/uwcc-admin/event_sim_utils/corrected_rf.csv")
    """
    # -s  --start_time    Rain timeseries start time (e.g: "2019-06-05 00:00:00"). Default is 23:30:00, 3 days before today.
    # -e  --end_time      Rain timeseries end time (e.g: "2019-06-05 23:30:00"). Default is 23:30:00, tomorrow.

    print(usageText)


if __name__=="__main__":

    set_db_config_file_path(os.path.join(ROOT_DIR, 'db_adapter_config.json'))

    try:
        # start_time = None
        # end_time = None
        flo2d_model = None
        method = "FF"
        grid_tag = "MDPA"  # note - grid tag "TP" has not handled yet
        file_path = None

        try:
            opts, args = getopt.getopt(sys.argv[1:], "h:m:s:e:M:g:f:",
                                       ["help", "flo2d_model=", "start_time=", "end_time=", "method=", "grid_tag=", "file="])
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
            elif opt in ("-M", "--method"):
                method = arg.strip()
            elif opt in ("-g", "--grid_tag"):
                grid_tag = arg.strip()
            elif opt in ("-f", "--file"):
                file_path = arg.strip()

        if file_path is not None:
            if not os.path.exists(file_path):
                print('Unable to find file : ', file_path)
                traceback.print_exc()
                exit(1)
        else:
            print("File path is not specified.")
            exit(1)

        if grid_tag == "MDPA":
            grid_interpolation = GridInterpolationEnum.getAbbreviation(GridInterpolationEnum.MDPA)
        elif grid_tag == "TP":
            grid_interpolation = GridInterpolationEnum.getAbbreviation(GridInterpolationEnum.TP)
            method = method + "_" + grid_interpolation
        else:
            print("Grid interpolation method should be either \"MDPA\" or \"TP\"")
            exit(1)

        if flo2d_model is None:
            print("Flo2d model is not specified.")
            exit(1)
        elif flo2d_model not in ("flo2d_250", "flo2d_150", "flo2d_150_v2"):
            print("Flo2d model should be either \"flo2d_250\" or \"flo2d_150\" or \"flo2d_150_v2\"")
            exit(1)

        # if start_time is None:
        #     start_time = (datetime.now() - timedelta(days=3)).strftime('%Y-%m-%d 23:30:00')
        # else:
        #     check_time_format(time=start_time, model=flo2d_model)
        #
        # if end_time is None:
        #     end_time = (datetime.now() + timedelta(days=1)).strftime('%Y-%m-%d 23:30:00')
        # else:
        #     check_time_format(time=end_time, model=flo2d_model)

        if flo2d_model == FLO2D_250:
            timestep = 5
        elif flo2d_model == FLO2D_150:
            timestep = 15
        elif flo2d_model == FLO2D_150_V2:
            timestep = 15

        corrected_rf_df = pd.read_csv(file_path, delimiter=',')
        mean_rf = corrected_rf_df.groupby('time').mean()['WRF_A'].reset_index().values.tolist()

        distinct_stations = corrected_rf_df.groupby(['longitude', 'latitude']).size()

        points_dict = {}
        print('distinct_stations')
        count = 1000
        for index, row in distinct_stations.iteritems():
            points_dict['point_{}'.format(count)] = [index[0], index[1]]
            count += 1

        if grid_tag == "MDPA":
            print("{} : ####### Insert rainfall from file to {} grids".format(datetime.now(), flo2d_model))
            # update_rainfall_from_file(flo2d_model=flo2d_model, method=method, grid_interpolation=grid_interpolation,
            #                           timestep=timestep, start_time=start_time, end_time=end_time)
        elif grid_tag == "TP":
            shape_file_path = os.path.join(ROOT_DIR, 'shape_files/Kalani_basin_hec_wgs/Kalani_basin_hec_wgs.shp')

            output_shape_file_path = os.path.join(ROOT_DIR, 'shape_files/output', "{}_out_shp.shp".format(
                (datetime.now()).strftime("%Y-%m-%d_%H-%M-%S")))

            polygons = get_voronoi_polygons(points_dict=points_dict, shape_file=shape_file_path, shape_attribute=['OBJECTID_1', 1],
                                 output_shape_file=output_shape_file_path, add_total_area=True)

            flo2d_grid_polygon_map = divide_flo2d_grids_to_polygons(flo2d_model=flo2d_model, polygons=polygons)

            print("{} : ####### Insert rainfall from file to {} grids".format(datetime.now(), flo2d_model))
            update_rainfall_from_file(flo2d_grid_polygon_map=flo2d_grid_polygon_map, stations_dict=points_dict,
                                      rainfall_df=corrected_rf_df, mean_rf=mean_rf, flo2d_model=flo2d_model, method=method,
                                      grid_interpolation=grid_interpolation, timestep=timestep)

    except Exception as e:
        traceback.print_exc()
    finally:
        print("Process finished.")
