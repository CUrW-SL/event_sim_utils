import os
import json
from datetime import datetime, timedelta
import pandas as pd
import numpy as np
import geopandas as gpd
from scipy.spatial import Voronoi
from shapely.geometry import Polygon, Point

from db_adapter.csv_utils import read_csv

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


# def assign_timeseries_to_flo2d_grids():


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

    pd.set_option('display.max_rows', df.shape[0] + 1)
    pd.set_option('display.max_columns', df.shape[1] + 1)
    print(df)
    return df

ROOT_DIR = "/home/shadhini/dev/repos/curw-sl/event_sim_utils"
file_path = "/home/shadhini/dev/repos/curw-sl/event_sim_utils/corrected_rf.csv"

corrected_rf_df = pd.read_csv(file_path, delimiter=',')

# mean_rf_df = corrected_rf_df.groupby('time').mean()['WRF_A'].reset_index().values.tolist()
#
# # print((mean_rf_df['WRF_A']).reset_index().values.tolist())
# print(mean_rf_df)
# df = pd.DataFrame({'A': [1, 1, 2, 1, 2],
#                    'B': [np.nan, 2, 3, 4, 5],
#                    'C': [1, 2, 1, 1, 2]}, columns=['A', 'B', 'C'])
# print(df)
#
# grouped_df = df.groupby('A').mean()
# print(grouped_df)
#
# grouped_df['mean'] =
distinct_stations = corrected_rf_df.groupby(['longitude', 'latitude']).size()

points_dict = {}
print('distinct_stations')
count = 1000
for index, row in distinct_stations.iteritems():
    points_dict['point_{}'.format(count)] = [index[0], index[1]] # ['longitude', 'latitude']
    count += 1

shape_file_path = os.path.join(ROOT_DIR, 'shape_files/250m_model/250m_model.shp')

output_shape_file_path = os.path.join(ROOT_DIR, 'shape_files/output_temp', "{}_out_shp.shp".format(
    (datetime.now()).strftime("%Y-%m-%d_%H-%M-%S")))

polygons = get_voronoi_polygons(points_dict=points_dict, shape_file=shape_file_path, shape_attribute=['Id', 0],
                                 output_shape_file=output_shape_file_path, add_total_area=True)

flo2d_model = "flo2d_250"
flo2d_grids = read_csv(
    os.path.join(ROOT_DIR, 'grids/flo2d/{}m.csv'.format(flo2d_model)))  # [Grid_ ID, X(longitude), Y(latitude)]

for grid in flo2d_grids:
    point = Point(float(grid[1]), float(grid[2]))

    for index, row in polygons.iterrows():
        polygon = polygons.iloc[index]['geometry']
        if point.within(polygon):
            grid.append(polygons.iloc[index]['id'])
            continue

print(flo2d_grids)

count=0
for grid in flo2d_grids:

    if count > 3:
        exit(0)
    lat = grid[2]
    lon = grid[1]
    cell_id = grid[0]

    meta_data = {
        'latitude': float('%.6f' % float(lat)), 'longitude': float('%.6f' % float(lon)),
        'model': flo2d_model, 'method': "MME",
        'grid_id': '{}_{}_{}'.format(flo2d_model, "TP", (str(cell_id)).zfill(10))
    }


    if len(grid) > 3:
        count += 1
        print(json.dumps(meta_data))

        polygon = grid[3]

        poly_lat = points_dict.get(polygon)[1]
        poly_lon = points_dict.get(polygon)[0]

        timeseries = corrected_rf_df.loc[(corrected_rf_df['latitude'] == poly_lat) & (corrected_rf_df['longitude'] == poly_lon)]

        print(timeseries.values.tolist())