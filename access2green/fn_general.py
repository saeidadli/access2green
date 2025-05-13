"""
Access to green spaces function.
"""

import pathlib

import networkx as nx
import osmnx as ox
import pandas as pd
import geopandas as gpd

from shapely.geometry import Point, LineString, Polygon, MultiPolygon
from shapely.ops import transform


#give network as a variable to cachment function
def access2green(
    tazs, #give TAZs polygons as geopandas
    uid, #Column with unique id for each zone
    buffer_distance_meters = 1000, #buffers the traffic zones by 1000 meter before downloading osm
    trip_time = 15, #in minutes
    travel_speed = 4.5, #walking speed in km/hour
):
    """
    0- find the CRS of the TAZs and find the UTM of the TAZs
    1- downlaod network based on the TAZs in UTM system
    2- convert TAZs to centroid
    3- create catchments for each centroid and convert it to geodatabase
    4- download the parks in utm, combine them, give them uid and area
    5- download the walking network in utm
    6- intersect the walking and parks and keep the uid and area. This will be our entrances
    7- find all entrances per catchment
    8 - remove duplicate entrances in each catchmet based on UiD
    9- Sum area of all unique entrances for each catchment
    10- map out the results
    """

    assert tazs.crs.name=='WGS 84', "Error message: input GeoDataFrame does't have CRS"

    #find utm of TAZs
    utm_crs = tazs.estimate_utm_crs()
    
    #define the download area
    bbox = tazs.total_bounds
    bbox_polygon = Polygon([(bbox[0], bbox[1]), (bbox[0], bbox[3]), (bbox[2], bbox[3]), (bbox[2], bbox[1])])

    #Buffer the TAZ
    buffer_distance_meters = 1000
    buffer_distance_degrees = buffer_distance_meters  / (111 * 1000)
    osm_download_area = bbox_polygon.buffer(buffer_distance_degrees)

    #Download the OSM under the buffer area
    cf = """
         ["area"!~"yes"]
         ["highway"]
         ["highway"!~"motor|proposed|construction|abandoned|platform|raceway"]
         ["foot"!~"no"]
         ["service"!~"private"]
         ["access"!~"private"]
         """
    G = ox.graph_from_polygon(
        polygon = osm_download_area, 
        network_type = 'walk',
        custom_filter = cf,
        retain_all=False,
    )

    #project the graph to UTM
    G = ox.project_graph(G, to_crs=utm_crs)

    # add an edge attribute for time in minutes required to traverse each edge
    meters_per_minute = travel_speed * 1000 / 60 #km per hour to m per minute
    for u, v, k, data in G.edges(data=True, keys=True):
        data['time'] = data['length'] / meters_per_minute

    # Create a new GeoDataFrame with centroids
    centroids_gdf = tazs.copy()
    centroids_gdf = centroids_gdf.to_crs(utm_crs) #use a projected crs for accuracy
    centroids_gdf['geometry'] = centroids_gdf['geometry'].centroid
    centroids_gdf = centroids_gdf.assign(nearest_node = lambda r: (ox.nearest_nodes(G, r['geometry'].x, r['geometry'].y)))

    #for each taz centroid find isopoly
    def make_iso_polys(G, center_node, trip_times, edge_buff=70, node_buff=70, infill=True):
        isochrone_polys = []
        for trip_time in sorted(trip_times, reverse=True):
            subgraph = nx.ego_graph(G, center_node, radius=trip_time, distance="time")

            node_points = [Point((data["x"], data["y"])) for node, data in subgraph.nodes(data=True)]
            nodes_gdf = gpd.GeoDataFrame({"id": list(subgraph.nodes)}, geometry=node_points)
            nodes_gdf = nodes_gdf.set_index("id")

            edge_lines = []
            for n_fr, n_to in subgraph.edges():
                f = nodes_gdf.loc[n_fr].geometry
                t = nodes_gdf.loc[n_to].geometry
                edge_lookup = G.get_edge_data(n_fr, n_to)[0].get("geometry", LineString([f, t]))
                edge_lines.append(edge_lookup)

            n = nodes_gdf.buffer(node_buff).geometry
            e = gpd.GeoSeries(edge_lines).buffer(edge_buff).geometry
            all_gs = list(n) + list(e)
            new_iso = gpd.GeoSeries(all_gs).unary_union

            # try to fill in surrounded areas so shapes will appear solid and
            # blocks without white space inside them
            if infill:
                new_iso = Polygon(new_iso.exterior)
            isochrone_polys.append(new_iso)
        return isochrone_polys
    
    centroids_gdf['geometry'] = centroids_gdf['nearest_node'].map(lambda n: make_iso_polys(G, n, [trip_time])[0])

    # download green space boundaries
    tags={
        'leisure': ['park', 'recreation_ground'],
        'natural': ['wood', 'scrub'],
        'landuse': ['grass']
    }

    green_spaces_gdf = ox.features_from_polygon(osm_download_area, tags=tags)
    dissolved  = green_spaces_gdf.unary_union
    green_spaces_gdf = gpd.GeoDataFrame(geometry=[dissolved], crs=4326)
    green_spaces_gdf = green_spaces_gdf.explode(index_parts=False)
    green_spaces_gdf = green_spaces_gdf.to_crs(utm_crs)
    green_spaces_gdf['area_sqm'] = green_spaces_gdf.geometry.area
    green_spaces_gdf = green_spaces_gdf.reset_index(drop=True)
    green_spaces_gdf['unique_id'] = range(1, len(green_spaces_gdf) + 1)
    green_spaces_gdf['geometry'] = green_spaces_gdf.boundary
    green_spaces_gdf = green_spaces_gdf[green_spaces_gdf.geometry.type.isin(['MultiLineString', 'LineString'])].copy()

    #create the walking netwrok from ox graph
    walking_network_gdf = ox.graph_to_gdfs(G, nodes=False, edges=True)

    # Intersect walking network and green spaces
    overlap_gdf = gpd.overlay(green_spaces_gdf, walking_network_gdf[['geometry']], how='intersection', keep_geom_type=False)
    overlap_gdf = overlap_gdf.explode(index_parts=False)
    overlap_gdf['geometry'] = overlap_gdf.representative_point()

    joined_df = gpd.sjoin(overlap_gdf, centroids_gdf , how="inner")
    joined_df = joined_df.drop_duplicates(subset=[uid,'unique_id'])

    aggregated_df = joined_df.groupby(uid)["area_sqm"].sum()  # Assuming "index_right" is the polygon index
    aggregated_df = aggregated_df.reset_index(name="total_green_area")

    return tazs.merge(aggregated_df, how='left')
