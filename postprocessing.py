import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as patches
import numpy as np
import cartopy.crs as ccrs
import timeit
from matplotlib.patches import Rectangle
import matplotlib.colors
import csv
import shapefile
import cartopy

#check runtime
start = timeit.default_timer()

#import simulation data
data_xarray0 = xr.concat(objs=(xr.open_dataset('E:/Run_27-05_Tides/Domburg_2020-01-01_Stokes=0_7300.nc'), xr.open_dataset('E:/Run_27-05_Tides/Domburg_2020-01-01_Stokes=0_73002.nc')), dim='obs')
data_xarray1 = xr.concat(objs=(xr.open_dataset('E:/Run_27-05_Tides/Domburg_2020-01-01_Stokes=0_7301.nc'), xr.open_dataset('E:/Run_27-05_Tides/Domburg_2020-01-01_Stokes=0_73012.nc')), dim='obs')
data_xarray2 = xr.concat(objs=(xr.open_dataset('E:/Run_27-05_Tides/Domburg_2020-01-01_Stokes=0_7302.nc'), xr.open_dataset('E:/Run_27-05_Tides/Domburg_2020-01-01_Stokes=0_73022.nc')), dim='obs')

#data_xarray0 and data_xarray1 can be concatenated directly (same size)
#data_xarray2 has not the same amount of particle released (release for only one year)
#so data_xarray2 will be concatenated after some processing (line 86)
data_xarray = xr.concat(objs=(data_xarray0, data_xarray1), dim='traj')

runtime = 730

basepath = r'C:/Users/bramv/OneDrive - Universiteit Utrecht/UU/Jaar 3/BONZ/Datafiles/'
outpath = r'C:/Users/bramv/OneDrive - Universiteit Utrecht/UU/Jaar 3/BONZ/Outputfiles/'
#load required input data
file_landmask = '/datafile_landMask_297x_375y'
landdata = np.genfromtxt(file_landmask, delimiter=None)

file_coast = basepath + 'datafile_coastMask_297x_375y'
coastdata  = np.genfromtxt(file_coast, delimiter=None)

file_fishing = basepath + 'datafile_fishingInputMatrices_297x_375y'
prior_fishing = np.genfromtxt(file_fishing, delimiter=None)
fishingmask = (landdata == 0)

file_popmatrices = basepath + 'netcdf_populationInputMatrices_thres50_297x_375y.nc'
popmatrix = xr.open_mfdataset(file_popmatrices)
popmatrix_2020 = (popmatrix['pop_input'].values)[4,:,:]
prior_coastal = popmatrix_2020 * coastdata

current_data = xr.open_mfdataset('C:/Users/bramv/documents/CMEMS/*.nc')
lons = current_data.coords['longitude'].values
lats = current_data.coords['latitude'].values
fieldMesh_x,fieldMesh_y = np.meshgrid(lons,lats)

#load division of target regions
#ordered arounda axis 2; order: UK_E, UK_SW, UK_SE, NL, BE, FR_N, FR_Brit, SC, IR, Other
coastalregions = np.load('coastalregions.npy')

#ordered arounda axis 2; order: Ch_W, Ch_E, NL_Sea, NorthSea, Other
fisheryregions = np.load('fisheryregions.npy')

#create bins based on grid cells, for histograms
xbins = np.linspace(-20, 13, 298)
ybins = np.linspace(40, 65, 376)   

#create colormap for plotting river data
num_colors = 9
cmap_r = plt.get_cmap('Greys', num_colors)
cmap_r2 = matplotlib.colors.ListedColormap(['white', 'black'])

#ignore observations after 730 days, make sure you follow every particle for the same time
x_raw = data_xarray['lon'].values
y_raw = data_xarray['lat'].values
x2_raw = data_xarray2['lon'].values
y2_raw = data_xarray2['lat'].values

x1 = np.empty((len(x_raw), int(runtime)))
y1 = np.empty((len(x_raw), int(runtime)))
x2 = np.empty((len(x2_raw), int(runtime)))
y2 = np.empty((len(x2_raw), int(runtime)))

for i in range(len(x_raw)):
    try: 
        x1[i] = (x_raw[i, np.isfinite(x_raw[i])==True])[:runtime]
        y1[i] = (y_raw[i, np.isfinite(y_raw[i])==True])[:runtime]
    #this ignores particles that are out of bounds
    except ValueError:
        pass

for i in range(len(x2_raw)):
    try: 
        x2[i] = (x2_raw[i, np.isfinite(x2_raw[i])==True])[:runtime]
        y2[i] = (y2_raw[i, np.isfinite(y2_raw[i])==True])[:runtime] 
    #this ignores particles that are out of bounds
    except ValueError:
        pass

#concatenate all trajectories
x = np.concatenate([x1[:,:runtime], x2[:,:runtime]])
y = np.concatenate([y1[:,:runtime], y2[:,:runtime]])

stop = timeit.default_timer()
print('Time cell loading variables: ', stop - start)  
#%% Load riverdata
start = timeit.default_timer()

def riverData():
    riverShapeFile    = basepath + 'Riverdata_2021/Meijer2021_midpoint_emissions/Meijer2021_midpoint_emissions.shp'
    pollutionFile        = basepath + 'Riverdata_2021/Meijer2021_midpoint_emissions.csv'
    dataArray_ID = 1 #column with yearly waste discharged by river
    
    sf = shapefile.Reader(riverShapeFile)
    
    #extract files within NorthSea
    plottingDomain = [-8.3, 5, 47, 57]
    
    rivers = {}
    rivers['longitude'] = np.array([])
    rivers['latitude'] = np.array([])
    rivers['ID'] = np.array([],dtype=int)
    rivers['dataArray'] = np.array([])
    
    for i1 in range(len(sf.shapes())):
        long = sf.shape(i1).points[0][0]
        lat = sf.shape(i1).points[0][1]
        
        if plottingDomain[0] < long <plottingDomain[1] and plottingDomain[2] < lat < plottingDomain[3]:
            rivers['longitude'] = np.append(rivers['longitude'],long)
            rivers['latitude'] = np.append(rivers['latitude'],lat)
            rivers['ID'] = np.append(rivers['ID'],i1)
            
            
    with open(pollutionFile, 'r',encoding='ascii') as csvfile:
        filereader = csv.reader(csvfile, delimiter=',')
        i1 = 0
        for row in filereader:
            if i1 > 0:
                data_ID = i1-1 
                if i1 == 1:
                    dataArray = [float(row[i2].replace(',','.')) for i2 in range(len(row))]
                    rivers['dataArray'] = dataArray
                else:
                    if data_ID in rivers['ID']:
                        dataArray = [float(row[i2].replace(',','.')) for i2 in range(len(row))]
                        rivers['dataArray'] = np.vstack([rivers['dataArray'],dataArray])
            i1 += 1
    
    coastIndices = np.where(coastdata == 1)
    assert(np.shape(coastIndices)[0] == 2), "coastMask.data should be an array where the first dimension of the three is empty"
    
    # array containing indices of rivers not belonging to North Sea, which are to be deleted
    deleteEntries = np.array([],dtype=int)
    
    # matrix corresponding to fieldmesh, with per coastal cell the amount of river pollution
    riverInputMatrix = np.zeros(fieldMesh_x.shape)
    
    # for every river
    for i1 in range(len(rivers['longitude'])):   
        lon_river = rivers['longitude'][i1]
        lat_river = rivers['latitude'][i1]
        dist = 1e10
        # check which point is closest
        for i2 in range(np.shape(coastIndices)[1]):
            lon_coast = lons[coastIndices[1][i2]]
            lat_coast = lats[coastIndices[0][i2]]
        
            lat_dist = (lat_river - lat_coast) * 1.11e2
            lon_dist = (lon_river - lon_coast) * 1.11e2 * np.cos(lat_river * np.pi / 180)
            dist_tmp = np.sqrt(np.power(lon_dist, 2) + np.power(lat_dist, 2))
            
            # save closest distance
            if dist_tmp < dist:
                dist = dist_tmp
                lat_ID = coastIndices[0][i2]
                lon_ID = coastIndices[1][i2]
            
        # if distance to closest point > threshold (3*approx cell length), delete entry
        if dist > 3*0.125*1.11e2:
            deleteEntries = np.append(deleteEntries,i1)
        # else: get pollution river, and add to releasematrix
        else:
            # add plastic input as obtained from the dataset
            riverInputMatrix[lat_ID,lon_ID] += rivers['dataArray'][i1,dataArray_ID]
            
    # rivers ending in NWSELF
    rivers_NWSHELF = {}
    rivers_NWSHELF['longitude'] = np.delete(rivers['longitude'],deleteEntries)
    rivers_NWSHELF['latitude'] = np.delete(rivers['latitude'],deleteEntries)
    rivers_NWSHELF['ID'] = np.delete(rivers['ID'],deleteEntries)
    rivers_NWSHELF['dataArray'] = np.delete(rivers['dataArray'],deleteEntries,axis=0)
    return rivers_NWSHELF, riverInputMatrix

rivers_NWSHELF, prior_river = riverData()
#ignore rivers with a plastic input smaller than 1 tonne/year, to prevent plots from being too crowded
prior_river[prior_river<1] = 0

stop = timeit.default_timer()
print('Time in cell loading river data: ', stop - start)  
#%% Thicken landborder, for plotting
def thickenCoast(coastalprobs, thickness):
    def getLandBorder(landMask,lon,lat,val_add): 
        n_lat = landMask.shape[0]
        n_lon = landMask.shape[1]
            
        for i1 in range(n_lat):
            for i2 in range(n_lon):
                
                check_bot = True
                check_top = True
                check_left = True
                check_right = True
                
                # check whether land is located at boundary
                if i1 == 0:
                    check_top = False
                if i1 == n_lat-1:
                    check_bot = False
                if i2 == 0:
                    check_left = False
                if i2 == n_lon-1:
                    check_right = False
                    
                # check whether cell is land, if so look for coast
                if landMask[i1,i2] == 1:
                    
                    if check_top:
                        if (landMask[i1-1,i2] == 0) or (landMask[i1-1,i2] >= 2):
                            landMask[i1,i2] = -1
                    if check_bot:
                        if (landMask[i1+1,i2] == 0) or (landMask[i1+1,i2] >= 2):
                            landMask[i1,i2] = -1
                    if check_left:
                        if (landMask[i1,i2-1] == 0) or (landMask[i1,i2-1] >= 2):
                            landMask[i1,i2] = -1
                    if check_right:
                        if (landMask[i1,i2+1] == 0) or (landMask[i1,i2+1] >= 2):
                            landMask[i1,i2] = -1
        landMask[landMask == -1] = val_add         
        return landMask
    
    landMask = landdata.copy()
    coastMask = coastdata.copy()
    
    landBorder = landMask.copy()
    val_add = 2
    for i1 in range(thickness):
        landBorder = getLandBorder(landBorder,lons,lats,val_add)
        val_add += 1
    
    def closest_index(lat,lon,mask_test):
        distMat = 1e5 * np.ones(fieldMesh_x.shape)
        
        test_indices = np.where(mask_test == 1)
        
        distMat_lon = (lon - fieldMesh_x[test_indices[0],test_indices[1]])*1.11e2*0.63 #find distances coastal element w.r.t. ocean cells. 0.63 comes from lat=51deg (Zeeland)
        distMat_lat = (lat - fieldMesh_y[test_indices[0],test_indices[1]])*1.11e2
    
        distMat[test_indices[0],test_indices[1]] = np.sqrt(np.power(distMat_lon, 2) + np.power(distMat_lat, 2))    
        
        return np.where(distMat == distMat.min())[0][0],np.where(distMat == distMat.min())[1][0]
    
    ### interpolate beaching to closest coastal cell
    hist_beaching_coast = np.zeros(fieldMesh_x.shape)
    for i1 in range(len(lats)):
        for i2 in range(len(lons)):
            
            if coastalprobs[i1,i2] > 0:
                
                i_lat,i_lon = closest_index(lats[i1],lons[i2],coastMask)
                hist_beaching_coast[i_lat,i_lon] +=coastalprobs[i1,i2]
    
     ### go through the landborder defined above with increased width
    hist_beaching_extended = np.zeros(fieldMesh_x.shape)
    indices_border = np.where(landBorder > 1)            
    for i1 in range(len(indices_border[0])):
        lon_ = lons[indices_border[1][i1]]
        lat_ = lats[indices_border[0][i1]]
        i_lat,i_lon = closest_index(lat_,lon_,coastMask)
        
        hist_beaching_extended[indices_border[0][i1],indices_border[1][i1]] += hist_beaching_coast[i_lat,i_lon]
    return hist_beaching_extended
#%% find origincells, assuming discrete ages
start = timeit.default_timer()

def find_sourceprobs_age(x, y, lon, lat):
    coastal_problist = [[] for j in range(10)]
    fishery_problist = [[] for j in range(5)]
    for j in range(24):
        print(j)
        fisherycells = np.zeros((375, 297))
        coastalcells = np.zeros((375,297))
        rivercells = np.zeros((375,297))
        likelihood = np.zeros((375,297))
        for i in range(len(x)):
            #calculate likelihood per particle, add to total likelihood
            hist_particle, xedges, yedges = np.histogram2d(x[i, 30*j:30*(j+1) ], y[i,30*j:30*(j+1)], bins=[lon, lat])
            hist_particle = hist_particle.T
            likelihood += hist_particle
        #multiply with prior
        fisherycells = likelihood * prior_fishing
        coastalcells = likelihood * prior_coastal
        rivercells = likelihood * prior_river
        #normalise with P(B)
        fisheryprobs = 40 * fisherycells / np.nansum(fisherycells)
        coastalprobs = 50 * coastalcells / np.nansum(coastalcells)
        riverprobs = 10 * rivercells / np.nansum(rivercells)
        #thicken coast, for plotting
        hist_coastal = thickenCoast(coastalprobs, 3)
        
        #calculate inproduct between target regions and probabilities to find probability per grid cell in source region
        #sum all these grid cell probabilities to get the total per source region, save this probability
        for i in range(10):
            c_i = coastalprobs * coastalregions[:,:,i]
            coastal_problist[i].append(np.sum(c_i))
        for i in range(5):
            f_i = fisheryprobs * fisheryregions[:,:,i]
            fishery_problist[i].append(np.sum(f_i))
        #plotting part        
#        levels_mpw = np.logspace(np.log10(0.001), np.log10(1), 9)
#        levels_fish = np.logspace(np.log10(0.001), np.log10(1), 9)
#        levels_river = np.logspace(np.log10(0.01), np.log10(10), 9)
#        
#        fig,ax = plt.subplots(3)
#        X,Y = np.meshgrid(np.linspace(0,100,100),np.linspace(0,100,100))
#        plt1 = ax[0].contourf(X,Y,np.random.choice(levels_mpw,size=[100,100]),levels_mpw,cmap=plt.cm.Reds, norm=plt.cm.colors.LogNorm(), extend='both')
#        cbar1 = plt.colorbar(plt1)
#        plt2 = ax[1].contourf(X,Y,np.random.choice(levels_fish,size=[100,100]),levels_fish,cmap=plt.cm.Blues, norm=plt.cm.colors.LogNorm(), extend='both')
#        cbar2 = plt.colorbar(plt2)
#        plt.close()
#        plt3 = ax[2].contourf(X,Y,np.random.choice(levels_river,size=[100,100]),levels_river,cmap=plt.cm.Greys, norm=plt.cm.colors.LogNorm(), extend='both')
#        cbar3 = plt.colorbar(plt3)
#        
#        
#        fig = plt.figure(figsize=(15,10))
#        ax = plt.axes(projection=ccrs.PlateCarree())
#        ax.coastlines(resolution='10m')
#        ax.set_extent((-9, 5, 46, 57), ccrs.PlateCarree())
#        ax.add_feature(cartopy.feature.RIVERS)
#        ax.add_feature(cartopy.feature.LAND)
#        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
#                   color='gray', alpha=0.5, linestyle='-')
#        gl.xlabels_top = False
#        gl.ylabels_right = False
#        
#        plt.contourf(fieldMesh_x,fieldMesh_y,hist_coastal,levels=levels_mpw,extend='both',cmap=plt.cm.Reds, linewidths = 2.5, norm=plt.cm.colors.LogNorm())#,norm=plt.cm.colors.LogNorm())
#        
#        plt.contourf(fieldMesh_x,fieldMesh_y,fisheryprobs,levels=levels_fish,extend='both',cmap=plt.cm.Blues,norm=plt.cm.colors.LogNorm())#,norm=plt.cm.colors.LogNorm())
#        
#        riverprobs[riverprobs == 0 ] = 'nan'
#        
#        for i in range(len(lats)):
#            #im is the inside of the river dots, im2 is the black line around it
#            im = ax.scatter(fieldMesh_x[i],fieldMesh_y[i], c = riverprobs[i], cmap=cmap_r, vmin=1e-3, vmax=4, zorder=2, s=80, norm=matplotlib.colors.LogNorm())
#            im2 = ax.scatter(fieldMesh_x[i],fieldMesh_y[i], c = riverprobs[i], cmap=cmap_r2, vmin=1e-3, vmax=1.1e-3, zorder=1, s=100)
#        box = ax.get_position()
#        ax.set_position([1.35*box.x0, 2.75 * box.y0, box.width * 0.8, box.height * 0.8])
#        
#        cax1 = fig.add_axes([0.30, 0.22, 0.4, 0.02])
#        cax2 = fig.add_axes([0.30, 0.14, 0.4, 0.02])
#        cax3 = fig.add_axes([0.30, 0.06, 0.4, 0.02])
#               
#        cbar3 = plt.colorbar(plt1,cax=cax3,orientation='horizontal', ticks=[0.001,0.01,0.1,1])
#        cbar2 = plt.colorbar(plt2,cax=cax2,orientation='horizontal', ticks=[0.001,0.01,0.1,1])
#        cbar1 = plt.colorbar(plt3,cax=cax1,orientation='horizontal',ticks=[0.01, 0.1, 1, 10])
#        #plot the marker, to indicate beaching location
#        ax.scatter(3.4, 51.6, marker='x', c='y', s=80, zorder=3)
#
#        cax3.set_title(r'Coastal probabilities [%]')
#        cax2.set_title(r'Fishery probabilities [%]')  
#        cax1.set_title(r'River probabilities [%]') 
#
#        ax.set_title('Source probabilities for particles beaching near marker on {}'.format(np.datetime64('2020-01-01')-7*j))
#        plt.savefig(outpath + 'plots/{}'.format(j), bbox_inches='tight', dpi=150)
#        plt.close()
    return (fisherycells, coastalcells, riverprobs, coastal_problist, fishery_problist)

fisherycells, coastalcells, riverprobs, coastal_problist_age, fishery_problist_age = find_sourceprobs_age(x, y, xbins, ybins)

stop = timeit.default_timer()
print('Time in cell calculating & plotting probs per age: ', stop - start)  
#%% find origincells, per beaching week
start = timeit.default_timer()

def find_sourceprobs_temp(x, y, lon, lat):
    coastal_problist = [[] for j in range(10)]
    fishery_problist = [[] for j in range(5)]
    #261 the amount of weeks in 5 years
    for j in range(261):
        print(j)
        fisherycells = np.zeros((375, 297))
        coastalcells = np.zeros((375,297))
        rivercells = np.zeros((375,297))
        likelihood = np.zeros((375,297))
        #only consider particles in corresponding beaching week
        for i in range(7*100*j, 7*100*(j+1)): 
            #calculate likelihood per particle, add to total likelihood
            #ignore first 2 observations (days), as particle is then still in Domburg (not the source)
            hist_particle, xedges, yedges = np.histogram2d(x[i,2:], y[i,2:], bins=[lon, lat])
            hist_particle = hist_particle.T
            likelihood += hist_particle
        #multiply with prior
        fisherycells = likelihood * prior_fishing
        coastalcells = likelihood * prior_coastal
        rivercells = likelihood * prior_river
        #normalise with P(B)
        fisheryprobs = 40 * fisherycells / np.nansum(fisherycells)
        coastalprobs = 50 * coastalcells / np.nansum(coastalcells)
        riverprobs = 10 * rivercells / np.nansum(rivercells)
        
        #calculate inproduct between target regions and probabilities to find probability per grid cell in source region
        for i in range(10):
            c_i = coastalprobs * coastalregions[:,:,i]
            coastal_problist[i].append(np.sum(c_i))
        for i in range(5):
            f_i = fisheryprobs * fisheryregions[:,:,i]
            fishery_problist[i].append(np.sum(f_i))
    return (fisherycells, coastalcells, riverprobs, coastal_problist, fishery_problist)

fisherycells_temp, coastalcells_temp, riverprobs_temp, coastal_problist_temp, fishery_problist_temp = find_sourceprobs_temp(x, y, xbins, ybins)

stop = timeit.default_timer()
print('Time in cell calculating & plotting probs per beaching week: ', stop - start)  
#%% Also calculate probabilities over whole runtime without averaging every week, for more accurate pie plot
#The total distribution fishery/coastal/river is 40/50/10, but not necessarily every week
start = timeit.default_timer()

def find_sourceprobs_av(x, y, lon, lat):
    fisherycells = np.zeros((375, 297))
    coastalcells = np.zeros((375,297))
    rivercells = np.zeros((375,297))
    likelihood = np.zeros((375,297))
    for i in range(len(x)): 
        #ignore first 2 observations (days), as particle is then still in Domburg (not the source)
        hist_particle, xedges, yedges = np.histogram2d(x[i, 2:], y[i,2:], bins=[lon, lat])
        hist_particle = hist_particle.T  
        likelihood += hist_particle
    fisherycells = likelihood * prior_fishing
    coastalcells = likelihood * prior_coastal
    rivercells = likelihood * prior_river
    return fisherycells, coastalcells, rivercells

fisherycells_pie, coastalcells_pie, rivercells_pie = find_sourceprobs_av(x, y, xbins, ybins)
#normalise fishery probabilities to 40 (assuming 40% of source plastic comes from fisheries) and coastal to 50  (prior * likelihood / normalisation)
fisheryprobs = 40 * fisherycells_pie / np.sum(fisherycells_pie)
coastalprobs = 50 * coastalcells_pie / np.sum(coastalcells_pie)
riverprobs = 10 * rivercells_pie / np.sum(rivercells_pie)

coastalprobs_total = thickenCoast(coastalprobs, 3)

#calculate probabilities over whole runtime per grid cell
coastal_probs_av = []
fishery_probs_av = []

for i in range(10):
    c_i = coastalprobs * coastalregions[:,:,i]
    coastal_probs_av.append(np.sum(c_i))
for i in range(5):
    f_i = fisheryprobs * fisheryregions[:,:,i]
    fishery_probs_av.append(np.sum(f_i))

stop = timeit.default_timer()
print('Time in cell calculating total probs: ', stop - start)  
#%% Plot total probabilities over whole runtime
start = timeit.default_timer()

#levels_mpw = np.logspace(np.log10(0.001), np.log10(1), 9)
#levels_fish = np.logspace(np.log10(0.001), np.log10(1), 9)
#levels_river = np.logspace(np.log10(0.01), np.log10(10), 9)
#fig,ax = plt.subplots(3)
#X,Y = np.meshgrid(np.linspace(0,100,100),np.linspace(0,100,100))
#plt1 = ax[0].contourf(X,Y,np.random.choice(levels_mpw,size=[100,100]),levels_mpw,cmap=plt.cm.Reds, norm=plt.cm.colors.LogNorm(), extend='both')
#cbar1 = plt.colorbar(plt1)
#plt2 = ax[1].contourf(X,Y,np.random.choice(levels_fish,size=[100,100]),levels_fish,cmap=plt.cm.Blues, norm=plt.cm.colors.LogNorm(), extend='both')
#cbar2 = plt.colorbar(plt2)
#plt.close()
#plt3 = ax[2].contourf(X,Y,np.random.choice(levels_river,size=[100,100]),levels_river,cmap=plt.cm.Greys, norm=plt.cm.colors.LogNorm(), extend='both')
#cbar3 = plt.colorbar(plt3)
# #thicken coast, for plotting
#coastalprobs_total = thickenCoast(coastalprobs, 3)
#
#fig = plt.figure(figsize=(15,10))
#ax = plt.axes(projection=ccrs.PlateCarree())
#ax.coastlines(resolution='10m')
#ax.add_feature(cartopy.feature.RIVERS)
#ax.add_feature(cartopy.feature.LAND)
#gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
#                   color='gray', alpha=0.5, linestyle='-')
#gl.xlabels_top = False
#gl.ylabels_right = False
#
#ax.set_extent((-15, 5, 46, 60), ccrs.PlateCarree())
#plt.contourf(fieldMesh_x,fieldMesh_y,coastalprobs_total,levels=levels_mpw,extend='both',cmap=plt.cm.Reds,norm=plt.cm.colors.LogNorm())
#plt.contourf(fieldMesh_x,fieldMesh_y,fisheryprobs,levels=levels_fish,extend='both',cmap=plt.cm.Blues, norm=plt.cm.colors.LogNorm())
#dont plot river probabilities below 0.1%, to prevent plot from being too crowded
#riverprobs[riverprobs < 1e-3 ] = 'nan'
#
#for i in range(len(lats)):
#    im = ax.scatter(fieldMesh_x[i],fieldMesh_y[i], c = riverprobs[i], cmap=cmap_r, vmin=1e-3, vmax=4, zorder=2, s=80, norm=matplotlib.colors.LogNorm())
#    #plot black border around the river scatter plot, for readability
#    im2 = ax.scatter(fieldMesh_x[i],fieldMesh_y[i], c = riverprobs[i], cmap=cmap_r2, vmin=1e-3, vmax=1.1e-3, zorder=1, s=100)
#
#box = ax.get_position()
#ax.set_position([1.35*box.x0, 2.75 * box.y0, box.width * 0.8, box.height * 0.8])
#
#cax1 = fig.add_axes([0.30, 0.22, 0.4, 0.02])
#cax2 = fig.add_axes([0.30, 0.14, 0.4, 0.02])
#cax3 = fig.add_axes([0.30, 0.06, 0.4, 0.02])
#
#
#cbar3 = plt.colorbar(plt1,cax=cax3,orientation='horizontal', ticks=[0.001,0.01,0.1,1])
#cbar2 = plt.colorbar(plt2,cax=cax2,orientation='horizontal', ticks=[0.001,0.01,0.1,1])
#cbar1 = plt.colorbar(plt3,cax=cax1,orientation='horizontal',ticks=[0.01, 0.1, 1, 10])
#ax.scatter(3.4, 51.6, marker='x', c='y', s=80, zorder=3) 
##ax.set_title('Source probabilities for particles beaching near marker between 2015 and 2020')
#cax3.set_title(r'Coastal probabilities [%]')
#cax2.set_title(r'Fishery probabilities [%]')   
#cax1.set_title(r'River probabilities [%]') 
##plt.savefig(outpath + 'totalprobs.jpg')
#
#stop = timeit.default_timer()
#print('Time in cell plotting total probs: ', stop - start)
#%% Bar chart of probabilities over time
coastal_problist_age = np.asarray(coastal_problist_age)
labels= ["UK E", "UK SE", "UK SW", "NL", "BE", "FR N", "FR Brit.", "SC", "IR", "Other (coastal)"]

fishery_problist_age = np.asarray(fishery_problist_age)
labels2 = ["Channel W", "Channel E", "NL", "North Sea", "Other (fishery)"]

coastal_problist_age[coastal_problist_age < 0] = 0
fishery_problist_age[fishery_problist_age < 0] = 0
data_coastal_norm = np.empty((10,24))
data_fishery_norm = np.empty((5,24))

#normalise to known beaching percentages
for j in range(24):
    data_coastal_norm[:,j] = 50*coastal_problist_age[:,j]/np.sum(coastal_problist_age[:,j])
    data_fishery_norm[:,j] = 40*fishery_problist_age[:,j]/np.sum(fishery_problist_age[:,j])      

fig, ax = plt.subplots()
ax.set_ylim(0,52)
X = np.arange(data_coastal_norm.shape[1])
for i in range(data_coastal_norm.shape[0]):
  ax.bar(X, data_coastal_norm[i],
    bottom = np.sum(data_coastal_norm[:i], axis = 0), label=labels[i])

# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.set_ylabel("Source probability [%]")
ax.set_xlabel("Assumed particle age [months]")

fig, ax = plt.subplots()
ax.set_ylim(0,42)
for i in range(data_fishery_norm.shape[0]):
  ax.bar(X, data_fishery_norm[i],
    bottom = np.sum(data_fishery_norm[:i], axis = 0), label=labels2[i])

# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

ax.set_ylabel("Source probability [%]")
ax.set_xlabel("Assumed particle age [months]")
#%%
data_coastal_year = np.empty((10,52))
data_fishery_year = np.empty((5,52))
data_coastal_norm = np.empty((10,261))
data_fishery_norm = np.empty((5,261))
coastal_problist_temp = np.asarray(coastal_problist_temp)
fishery_problist_temp = np.asarray(fishery_problist_temp)

for j in range(261):
    data_coastal_norm[:,j] = 50*coastal_problist_temp[:,j]/np.sum(coastal_problist_temp[:,j])
    data_fishery_norm[:,j] = 40*fishery_problist_temp[:,j]/np.sum(fishery_problist_temp[:,j])
    
yearly_coastal = np.empty((10,52,5))
yearly_fish = np.empty((5,52,5))

for i in range(260):
    weeknr = i % 52
    j = i // 52
    yearly_coastal[:,weeknr,j] = data_coastal_norm[:,i]
    yearly_fish[:, weeknr, j] = data_fishery_norm[:,i]
    
year_coastal = np.empty((10,52))
year_fish = np.empty((5,52))

for i in range(52):
    for j in range(10):
        year_coastal[j,i] = np.mean(yearly_coastal[j,i,:])

for i in range(52):
    for j in range(5):
        year_fish[j,i] = np.mean(yearly_fish[j,i,:])     
        
fig, ax = plt.subplots()
ax.set_ylim(0,52)
X = np.arange(year_coastal.shape[1])
for i in range(year_coastal.shape[0]):
  ax.bar(X, year_coastal[i],
    bottom = np.sum(year_coastal[:i], axis = 0), label=labels[i]) 
ax.set_xlabel("Beaching date")
ax.set_ylabel("Source probability [%]")
ax.set_xticks([0, 4, 8, 12, 16, 21, 25, 29, 34, 38, 43, 47])
ax.set_xticklabels(['Jan.', 'Feb.', 'Mar.', 'Apr.', 'May', 'Jun.', 'Jul.', 'Aug.', 'Sep.', 'Oct.', 'Nov.', 'Dec.'])
ax.legend()

fig, ax = plt.subplots()
ax.set_ylim(0,42)
X = np.arange(year_fish.shape[1])
for i in range(year_fish.shape[0]):
  ax.bar(X, year_fish[i],
    bottom = np.sum(year_fish[:i], axis = 0), label=labels2[i]) 
ax.set_xlabel("Beaching date")
ax.set_ylabel("Source probability [%]")
ax.set_xticks([0, 4, 8, 12, 16, 21, 25, 29, 34, 38, 43, 47])
ax.set_xticklabels(['Jan.', 'Feb.', 'Mar.', 'Apr.', 'May', 'Jun.', 'Jul.', 'Aug.', 'Sep.', 'Oct.', 'Nov.', 'Dec.'])
ax.legend()