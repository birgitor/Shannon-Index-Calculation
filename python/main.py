import numpy as np
import geopandas as gpd
from shapely.geometry import Polygon


########################
#calculate the grid
#######################

landuse = gpd.read_file('aspern_landuse.geojson') #read in the file with the landuse polygons

def index_calculation (landuse):

    landuse_repr = landuse.to_crs("EPSG:3857") #reproject it to Webmercator, so the unit is meter and not degree

    xmin, ymin, xmax, ymax = landuse_repr.total_bounds #set the boundaries for the grid file

    length = 50 #set the grids lenght and width (50m)
    wide = 50

    cols = list(np.arange(xmin, xmax + wide, wide))  #create a list of coordinates
    rows = list(np.arange(ymin, ymax + length, length))


    polygons = []
    for x in cols[:-1]:
        for y in rows[:-1]:
            polygons.append(Polygon([(x,y), (x+wide, y), (x+wide, y+length), (x, y+length)]))


    grid = gpd.GeoDataFrame(crs="EPSG:3857", geometry=polygons) #form a geoDataFrame out of the polygons coordinates
    #print (grid.head())


    ########################
    #intersect the grid with the landuse, so that each grid cell has a landuse attribute
    #######################

    #add a unique id
    id_values=[]
    for i in range(len(grid.index)):
        id_values.append(i)

    grid.insert(loc=0, column="id", value=id_values)

    #intersect grid and layer

    _overlay = gpd.overlay(grid, landuse_repr, how='intersection')

    _overlay['area'] = _overlay['geometry'].area/10**6
    _overlay.sort_values(by=['area'], inplace=True)
    _overlay.drop_duplicates(subset='id', keep='last', inplace=True)

    _overlay.drop(['area', 'geometry'], inplace=True, axis=1)

    grid_use = grid.join(_overlay.set_index('id'), on='id')

    # Drop rows with None/NaN

    grid_use = grid_use[grid_use.NUTZUNG_LEVEL2.notnull()]



    ########################
    #form centeroids out of grid
    #######################

    center = grid_use.centroid #create centroids

    use_points= gpd.GeoDataFrame.copy(grid_use)  #duplicate grid_use

    use_points['geometry'] = center #assign center points as geometry to duplicated grid use
    # --> result is a geodataframe with centroids as geometry and the attributes of the according grid cell




    ########################
    # create buffer around centroids
    #######################

    centerbuf = gpd.GeoDataFrame.copy(use_points) #variable for buffer polygons
    centerp = gpd.GeoDataFrame.copy(centerbuf) #variable for cell center points

    centerbuf['geometry'] =centerbuf.buffer(500) #calculate 500m buffer around center points


    ########################
    #calculate index value for each center point now, and then add the values to the cell (grid use)
    ########################


    grid_use['div_index']=""
    #print(grid_use.head())

    print ("last step starts here")

    for i in centerp['id']: #iterate through all buffers (buffered points)
        pointsInPoly = gpd.tools.sjoin(centerp, centerbuf.loc[centerbuf['id']==i], predicate="within", how='inner') #for each buffered polygon, check which points are in there, add this point to a new variable
        erholung = (sum(pointsInPoly['NUTZUNG_LEVEL2_left']== 'Erholungs- u. Freizeiteinrichtungen'))/(len(pointsInPoly))
        betrieb = (sum(pointsInPoly['NUTZUNG_LEVEL2_left']== 'Geschäfts,- Kern- und Mischnutzung (Schwerpunkt betriebl. Tätigkeit)'))/(len(pointsInPoly))
        gewasser = (sum(pointsInPoly['NUTZUNG_LEVEL2_left']== 'Gewässer'))/(len(pointsInPoly))
        industrie = (sum(pointsInPoly['NUTZUNG_LEVEL2_left']== 'Industrie- und Gewerbenutzung'))/(len(pointsInPoly))
        landwirtschaft = (sum(pointsInPoly['NUTZUNG_LEVEL2_left']== 'Landwirtschaft'))/(len(pointsInPoly))
        natur = (sum(pointsInPoly['NUTZUNG_LEVEL2_left']== 'Naturraum'))/(len(pointsInPoly))
        sozInf = (sum(pointsInPoly['NUTZUNG_LEVEL2_left']== 'soziale Infrastruktur'))/(len(pointsInPoly))
        street = (sum(pointsInPoly['NUTZUNG_LEVEL2_left']== 'Straßenraum'))/(len(pointsInPoly))
        technik = (sum(pointsInPoly['NUTZUNG_LEVEL2_left']== 'Technische Infrastruktur/Kunstbauten/Sondernutzung'))/(len(pointsInPoly))
        verkehr = (sum(pointsInPoly['NUTZUNG_LEVEL2_left']== 'weitere verkehrliche Nutzungen'))/(len(pointsInPoly))
        wohnen = (sum(pointsInPoly['NUTZUNG_LEVEL2_left']== 'Wohn- u. Mischnutzung (Schwerpunkt Wohnen)'))/(len(pointsInPoly))
        if erholung != 0: #if there are no points assigned to a certain land use type, this causes a division by zero (logarithm of zero is not defined)
            p1= erholung *np.log(erholung)
        else: p1 = 0
        if betrieb != 0:
            p2= betrieb *np.log(betrieb)
        else: p2 = 0
        if gewasser != 0:
            p3= gewasser *np.log(gewasser)
        else: p3 =0
        if industrie !=0:
            p4= industrie *np.log(industrie)
        else: p4 = 0
        if landwirtschaft != 0:
            p5=landwirtschaft *np.log(landwirtschaft)
        else: p5 =0
        if natur !=0:
            p6= natur *np.log(natur)
        else: p6=0
        if sozInf !=0:
            p7=sozInf *np.log(sozInf)
        else: p7 =0
        if street !=0:
            p8= street*np.log(street)
        else: p8 =0
        if technik != 0:
            p9= technik*np.log(technik)
        else: p9 =0
        if verkehr != 0:
            p10= verkehr *np.log(verkehr)
        else: p10 =0
        if wohnen !=0:
            p11= wohnen*np.log(wohnen)
        else: p11 =0
        shannonValue= (-1)*(np.sum([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11]))
        grid_use.loc[grid_use['id']==i, ['div_index']]= shannonValue


    grid_use = grid_use.to_crs("EPSG:4326") #convert it back to WGS 84 so that it can be used in Mapbox

    grid_use['div_index'] = grid_use['div_index'].replace('.',',',regex=True)

    #grid_use.to_file('shan_index_grid.geojson', driver='GeoJSON')

    ################
    #aggregate on block level
    ################


    landuse['shannon']=""

    shan_average = grid_use.groupby(['blockID'], as_index=False).mean()


    for i in shan_average['blockID']:
        mean_value = shan_average.loc[shan_average['blockID'] == i, 'div_index'].values[0]
        landuse.loc[landuse['blockID'] == i, 'shannon'] = mean_value



    landuse['shannon'] = landuse['shannon'].replace('.',',')

    return landuse

endresult = index_calculation(landuse)

    #landuse.to_file('shan_layer.geojson', driver='GeoJSON')
