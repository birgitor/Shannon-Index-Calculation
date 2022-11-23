import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon
import math
import time

########################
#calculate the grid
#######################

landuse = gpd.read_file('aspern_landuse.geojson') #read in the file with the landuse polygons

totalTime = time.time()




def reproject (geodf, crs = "EPSG:3857"):
    # geodf = input geodataframe
    # crs = a projected crs (unit must be meter), default is EPSG:3857

    if geodf.crs == "EPSG:4326":
        geodf_reprj = geodf.to_crs(crs) # reproject, so the unit is meter and not degree


    return geodf_reprj





def createGrid(geodf_reprj, l=50,w=50):
    # geodf_reproj - a geodataframe, that has a projected crs (unit must be meter, like EPSG:3857)
	# l = lenght of the grid cell (default is 30m)
	# width of the grid cell (default is 50m)

    xmin, ymin, xmax, ymax = geodf_reprj.total_bounds #set the boundaries for the grid file

    length = l #set the grids lenght and width
    wide = w

    cols = list(np.arange(xmin, xmax + wide, wide))  #create a list of coordinates
    rows = list(np.arange(ymin, ymax + length, length))

    polygons = []
    for x in cols[:-1]:
        for y in rows[:-1]:
            polygons.append(Polygon([(x,y), (x+wide, y), (x+wide, y+length), (x, y+length)]))


    grid = gpd.GeoDataFrame(crs=geodf_reprj.crs, geometry=polygons) #form a geoDataFrame out of the polygons coordinates

    id_values=[] #add a unique id to the grid
    for i in range(len(grid.index)):
        id_values.append(i)

    grid.insert(loc=0, column="gridID", value=id_values)

    return grid


def polygon2grid(grid, geodf_reprj, uniqueID = None, removeNaN=True, outputGeo="points"):
    # grid = regular grid
    # geodf_reprj = a (projected) geodataframe containing polygons for intersection and a unique ID
    # uniqueID = name of the column with a unique ID; if nothing is specified, the ID is equal to the row number
    # removeNaN = remove the Null/NaN values, default = True
    # outputGeo = defines output geometry, gri cells ("gridCells") or points ("points", default)

    # set unique ID and remove all unnecessary columns
    if uniqueID == None:
        geodf_reprj.reset_index(inplace=True) # create a column with index as unique ID
        geodf_reprj = geodf_reprj.rename(columns={'index': 'polyID'}) # rename index column
        geodf_reprj = geodf_reprj[['polyID', 'geometry']]
    else:
        geodf_reprj = geodf_reprj[[uniqueID, 'geometry']]

    # do the intersection
    _overlay = gpd.overlay(grid, geodf_reprj, how='intersection') #the result is a new gdf, containing the geometry of the overlapping features

    # pick the largest share
    _overlay['area'] = _overlay['geometry'].area  # calculate the area of geometry
    _overlay.sort_values(by=['area'], inplace=True)  # sort the feature by size
    _overlay.drop_duplicates(subset='gridID', keep='last', inplace=True)  #drop features with same ID, only keep last (= largest)
    #_overlay.drop(['area', 'geometry'], inplace=True, axis=1) #drop area and geometry column, which are not needed anymore

    if uniqueID == None:
        _overlay.sort_values(by=['polyID'], inplace=True)
        matchTable = pd.DataFrame(_overlay[['polyID', 'gridID']])
        #matchDict = matchTable.groupby('polyID')['gridID'].apply(list).to_dict()  # convert the pd to a dict
    else:
        _overlay.sort_values(by=[uniqueID], inplace=True)
        matchTable = pd.DataFrame(_overlay[[uniqueID, 'gridID']])
        #matchDict = matchTable.groupby(uniqueID)['gridID'].apply(list).to_dict()  # convert the pd to a dict



    return matchTable



def pointsWithinDist (grid, matchTable, radius = 500):
    #grid
    #matching table of intersection, pandas dataframe
    #radius = radius for buffer areas


    gridMatch = grid.loc[grid['gridID'].isin(matchTable['gridID'].values.tolist())] # only maintain grid cells that have a corresponding ID in the matching table


    center = gridMatch.centroid  # create centroids
    centroids = gpd.GeoDataFrame.copy(gridMatch)  # duplicate grid
    centroids['geometry'] = center  # assign center points as geometry to duplicated grid
    # --> result is a geodataframe with centroids as geometry and the attributes of the according grid cell

    '''
    centerbuf = gpd.GeoDataFrame.copy(centroids)  # variable for buffer polygons
    centerbuf['geometry'] = centerbuf.buffer(radius)  # calculate 500m buffer around center points
    centerbuf.rename(columns={"gridID": "bufID"}, inplace=True) # rename unique ID
    '''

    centroids['x'] = centroids.geometry.x
    centroids['y'] = centroids.geometry.y

    # --> could probably use geoDataFrame functionality  
    # --> sofar it takes 6 seconds 

    centroids = pd.DataFrame(centroids[['gridID', 'x', 'y' ]])
    pointWithinRadDict={}
    for rowA in centroids.itertuples():
        validPoints =[]
        for rowB in centroids.itertuples():
            diffx = rowA.x - rowB.x
            diffy = rowA.y - rowB.y
            dist = math.sqrt(diffx ** 2 + diffy ** 2)
            if dist <= radius:
                validPoints.append(rowB.gridID)
        pointWithinRadDict.update({rowA.gridID :validPoints})


    return pointWithinRadDict




def ShannonIndex (grid, geodf_reprj, matchTable, pointWithinRadDict, landUseCol,uniqueID = None,  outputGeom = "polygon"):
    # grid = regular grid
    # geodf_reprj = the reprojected geodataframe (polygon layer) with the landuse information
    # matchTable = geodataframe with the matching IDs between grid and geodf_reprj
    # uniqueID as definded in the matchTable
    # pointWithingRadDict = Dictionary with all point IDs as keys and all neighboring points within a certain radius
    # lanUseCol = column of geodf_reprj with landUse information that should be used for the Shannon calculation
    # outputGeom = which geometry is usd for Shannon Index calc, block/polygon level ("polygon", default) or grid cells ("grid")

    # something with grid -- > 0,0003 s
    grid = grid.loc[grid['gridID'].isin(matchTable['gridID'].values.tolist())]  # only maintain grid cells that have a corresponding ID in the matching table


    grid['div_index'] = "" # create empty column for shanin index values

    # matchDict transform --> 0.00007s
    matchDict = pd.Series(matchTable.OBJECTID.values, index=matchTable.gridID).to_dict() # convert Dataframe to Dict
    
    #geodf_reprj.set_index('OBJECTID', inplace=True) # set the column with the unique ID as Index
    print(geodf_reprj.shape)
    lu_vals = geodf_reprj[landUseCol].values
    print(len(list(lu_vals)))
    print ('for loops start here')

    landUseCol_ID = 0
    for i, curCol in enumerate(geodf_reprj.columns):
        print(curCol, landUseCol, i)
        if curCol == landUseCol:
            landUseCol_ID = i

    uniqueLanduseList = geodf_reprj[landUseCol].unique()
    uniqueLanduseDict = {str(lu):0 for lu in uniqueLanduseList}
    for key in pointWithinRadDict: # iterating through all keys ( means iterating through all gridIDs, for which a Shannon Index value will be calculated
        summands = []
        landUses = [] # could be a dictionary with "residential":450 <- squaremeter number 
        curLandUseDict = dict(uniqueLanduseDict)
        # takes 0.01s 
        
        for pointID in pointWithinRadDict[key]: # iterating through all the neighboring points (stored in the dict)
            polygonID = matchDict[pointID] # for the gridID (resp. point ID), get the corresponding polygon ID

            # much faster
            pointLandUse = lu_vals[25]
            #pointLandUse =  geodf_reprj.loc[polygonID, landUseCol] # geodf_reprj.loc[polygonID][landUseCol] #get the landuse type of the polygonID
            #landUses.append(pointLandUse) # append the landuse of the point to a list

            curLandUseDict[pointLandUse] += 1
    
        
        landuseCnt = sum(curLandUseDict.values()) #len(landUses)
        for landUseType in uniqueLanduseList: # for each landuse type in the original geodf_reprj
            landUseCount = curLandUseDict[landUseType] #landUses.count(landUseType) # count the number of points in the neighborhood with this landuse type
            if landUseCount != 0:
                shareLanduse = landUseCount / landuseCnt
                p1 = shareLanduse * np.log(shareLanduse) #calculate the summand of the land Use type for the index
                summands.append(p1)

        shannonValue = (-1) * (np.sum(summands)) #calculate the index value for the specific gridID (resp. point)
        grid.loc[grid['gridID'] == key, 'div_index'] = shannonValue #add value to the grid
        
    print(grid.head())



#######

# call the functions

#######

geodf_reprj = reproject(landuse)

curTime = time.time()
grid = createGrid(geodf_reprj)

# print time
print("grid creation took: ", curTime - time.time())


curTime = time.time()
matchTable = polygon2grid(grid,geodf_reprj, 'OBJECTID')

# print time
print("matchTable creation took: ", curTime - time.time())

curTume = time.time()
pointWithinRadDict = pointsWithinDist(grid,matchTable)

# print time
print("pointWithinRadDict creation took: ", curTime - time.time())


curTime = time.time()
ShannonIndex(grid, geodf_reprj, matchTable,pointWithinRadDict,'use_lvl2', 'OBJECTID')

# print time
print("shannon index loop took: ", curTime - time.time())


print("totalTime, ", time.time() - totalTime)