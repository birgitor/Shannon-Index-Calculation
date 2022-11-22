import numpy as np
import geopandas as gpd


landuse = gpd.read_file('aspern_landuse.geojson') #read in the file with the landuse polygons


########################################################################
#Initialization
########################################################################


########################
#reproject the geodataframe
#######################


def reproject (geodf, crs = "EPSG:3857"):
    # geodf = input geodataframe
    # crs = a projected crs (unit must be meter), default is EPSG:3857

    if geodf.crs == "EPSG:4326":
        geodf_reprj = geodf.to_crs(crs) # reproject, so the unit is meter and not degree


    return geodf_reprj


########################
# create centroids and buffer
#######################


def indexCalculation (geodf_reprj, radius = 500):

    # create centroids
    centroids = geodf_reprj.centroid

    polyCenter = gpd.GeoDataFrame.copy(geodf_reprj)  # duplicate polygon geodataframe
    polyCenter['geometry'] = centroids # assign the centroids as new geometry


    # create buffer areas with defined radius
    centerbuf = gpd.GeoDataFrame.copy(polyCenter)  # copy center points GeoDataFrame (--> variable for buffer polygons)
    centerbuf['geometry'] = centerbuf.buffer(radius)  # calculate 500m buffer around center points, assign new geometry to the copy
    centerbuf = gpd.GeoDataFrame(centerbuf[['OBJECTID', 'geometry']]) # leave only two columns with ID and geometry, drop the rest
    centerbuf.rename(columns={"OBJECTID": "bufID"}, inplace=True)  # rename unique ID, so that it is clear for the intersection

    # do the intersection
    _overlay = gpd.overlay(centerbuf, geodf_reprj,how='intersection')  # the result is a new gdf, containing the geometry of the overlapping features


    # calculate the area
    _overlay['area'] = _overlay['geometry'].area # area is important to weigh the landuse types


    #preperation for index calc
    geodf_reprj['div_index'] = ""  # create a column for the index



    ### INDEX CALCULATION ###

    for point in polyCenter['OBJECTID']: # for each polygon center point
        summands = []
        neighborhoodArea = sum(_overlay.loc[_overlay['bufID'] == point, ['area']].values[0:])  # calcualte the whole area of a neighboorhood (all areas within radius)
        neighborhoodObjects = _overlay[(_overlay['bufID'] == point)] #put all the polygons within the specific neighborhood of the current loop in a new GDF
        for landUseType in geodf_reprj['use_lvl2'].unique(): #iterate through laduse types
            landUseArea= []
            for row in neighborhoodObjects.itertuples(): #iterate throogh all objects (polygons) within the specific neighborhood of curent loop
                landUse = row.use_lvl2
                if landUse == landUseType: #if landuse of object equals the landuse type of the current loop
                    landUseArea.append(row.area)
            landUseSum = sum(landUseArea)
            if landUseSum != 0:
                shareLanduse = landUseSum / neighborhoodArea
                p1 = shareLanduse * np.log(shareLanduse)  # calculate the summand of the land Use type for the index
                summands.append(p1)
        shannonValue = (-1) * (np.sum(summands))  # calculate the index value for the specific gridID (resp. point)
        geodf_reprj.loc[geodf_reprj['OBJECTID'] == point, 'div_index'] = shannonValue  # add value to the polygon

    print (geodf_reprj.head())
    #geodf_reprj.to_file('shan_layer_new.geojson', driver='GeoJSON')




    return geodf_reprj

print("finish")
######################################################
### call the functions #####


geodf_reprj = reproject(landuse)

result = indexCalculation (geodf_reprj)

print ("here is the result")
print (result)