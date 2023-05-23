var startDate = '2019-03-01';
var endDate = '2019-05-31';

// Create an empty feature collection to store the results
var finalFeatureCollection = ee.FeatureCollection([]);

// Create an empty list to store intermediate results
var list = [];

function calculate_indices(image) {
    // Selecting the necessary bands from the image
    var blue = image.select('SR_B2');
    var green = image.select('SR_B3');
    var red = image.select('SR_B4');
    var nir = image.select('SR_B5');
    var swir = image.select('SR_B6');
    
    // Calculating NDVI: (NIR - Red) / (NIR + Red)
    var ndvi = image.expression('(nir-red)/(nir+red)',
          {'nir':nir.multiply(0.0000275).add(-0.2),
           'red':red.multiply(0.0000275).add(-0.2)}).rename('NDVI');
    
    // Calculating EVI: 2.5 * ((NIR - Red) / (NIR + 6 * Red - 7.5 * Blue + 1))
    var evi = image.expression('2.5*(nir-red)/(nir+6*red-7.5*blue+1)',
          { 'nir':nir.multiply(0.0000275).add(-0.2),
            'red':red.multiply(0.0000275).add(-0.2),
            'blue': blue.multiply(0.0000275).add(-0.2)}).rename('EVI');
    
    // Calculating SAVI: ((1+0.5) * (NIR - Red)) / (NIR + Red + 0.5)
    var savi = image.expression('(1+0.5)*(nir-red)/(nir+red+0.5)',
          { 'nir':nir.multiply(0.0000275).add(-0.2),
            'red':red.multiply(0.0000275).add(-0.2)}).rename('SAVI');
    
    // Calculating MNDWI: (Green - SWIR) / (Green + SWIR)
    var mndwi = image.expression('(green-swir)/(green+swir)',
          {'green':green.multiply(0.0000275).add(-0.2),
           'swir':swir.multiply(0.0000275).add(-0.2)}).rename('MNDWI');
    
    // Calculating NDMI: (NIR - SWIR) / (NIR + SWIR)
    var ndmi = image.expression('(nir-swir)/(nir+swir)',
          {'nir':nir.multiply(0.0000275).add(-0.2),
           'swir':swir.multiply(0.0000275).add(-0.2)}).rename('NDMI');
    
    // Calculating BSI: ((Red + SWIR) - (NIR + Blue)) / ((Red + SWIR) + (NIR + Blue))
    var bsi = image.expression('((red+swir) - (nir+blue))/((red+swir) + (nir+blue))',
          { 'nir':nir.multiply(0.0000275).add(-0.2),
            'swir':swir.multiply(0.0000275).add(-0.2),
            'red':red.multiply(0.0000275).add(-0.2),
            'blue':blue.multiply(0.0000275).add(-0.2)}).rename('BSI');
    
    // Calculating NDBI: (SWIR - NIR) / (SWIR + NIR)
    var ndbi = image.expression('(swir-nir)/(swir+nir)',
          {'swir':swir.multiply(0.0000275).add(-0.2),
           'nir':nir.multiply(0.0000275).add(-0.2)}).rename('NDBI');
    
    // Calculating the minimum and maximum NDVI values within the image
    var ndvi_min = ee.Number(ndvi.reduceRegion({
        reducer: ee.Reducer.min(),
        scale: 30,
        maxPixels: 1e11}).values().get(0));
  
    var ndvi_max = ee.Number(ndvi.reduceRegion({
        reducer: ee.Reducer.max(),
        scale: 30,
        maxPixels: 1e11}).values().get(0));
    
    // Fraction of Vegetation Index
    // Calculating FVC: ((NDVI - NDVI_min) / (NDVI_max - NDVI_min))^2
    var fvc = image.expression('((ndvi-ndvi_min)/(ndvi_max-ndvi_min))**2',
               {'ndvi': ndvi, 'ndvi_min': ndvi_min, 'ndvi_max': ndvi_max}).rename('FVC');
    
    // Clipping FVC values to the range [0, 1]
    fvc = fvc.where(fvc.lt(0.0), 0.0);
    fvc = fvc.where(fvc.gt(1.0), 1.0);
    
    // Calculating emissivity: (c1 * FVC + c2)
    var emissivity = image.expression('(c1*fvc + c2)', {'fvc': fvc, 'c1': 0.004, 'c2': 0.986}).rename('Emissivity');
           
    // Adding all the calculated indices as bands to the original image
    return image.addBands(ndvi).addBands(ndbi).addBands(bsi).addBands(ndmi).addBands(mndwi).addBands(evi).addBands(savi).addBands(fvc).addBands(emissivity);
  }
  
// Function to mask clouds and cloud shadows in Landsat TOA data
function maskL8toa(image) {
    // Define the bit masks for cloud shadow and clouds
    var cloudShadowBitMask = (1 << 3);
    var cloudsBitMask = (1 << 4);
    
    // Select the quality assessment (QA) band
    var qa = image.select('QA_PIXEL');
    
    // Create a mask by checking if the cloud shadow and clouds bits are set to 0
    var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
      .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
    
    // Apply the mask to the image and return the masked image
    return image.updateMask(mask);
  }
  
  // Function to mask clouds in Landsat SR data
  function maskL8sr(image) {
    // Select the quality assessment (QA) band
    var cfmask = image.select('QA_PIXEL');
    
    // Create a mask by checking specific bits in the QA band
    var clear = cfmask.bitwiseAnd(1 << 2).eq(0)
      .and(cfmask.bitwiseAnd(1 << 3).eq(0))
      .and(cfmask.bitwiseAnd(1 << 5).eq(0));
    
    // Apply the mask to the image and return the masked image
    return image.updateMask(clear);
  }
  
  // Load the Landsat indices image collection for the specified date range
  var landsat_indices = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
    .filterDate(startDate, endDate);
  
  // Load the Landsat TOA image collection for the specified date range
  var landsat_TOA = ee.ImageCollection('LANDSAT/LC08/C02/T1_TOA')
    .filterDate(startDate, endDate);

// For loop will iterate on all the grid points and calculate the indices
for(var i=1;i<=241;i++) {
  // Select the specific grid using its Grid_ID
  var geo = south.filter(ee.Filter.eq('Grid_ID', i));
  
  // Get the first image within the selected grid
  var image = geo.first();
  
  // Extract the geometry of the selected grid
  var roi = geo.geometry();
  
  // Landsat indices image collection
  var indices_collection = landsat_indices
    .filterBounds(roi)
    .map(maskL8sr)
    .map(calculate_indices)
    .mean()
    .clip(roi);
  
  // Landsat TOA image collection
  var lst_collection = landsat_TOA
    .filterBounds(roi)
    .map(maskL8toa)
    .mean()
    .clip(roi);
  
  // Compute statistics for NDVI
  var ndvi_stats = indices_collection
    .select('NDVI')
    .reduceRegion({
      reducer: ee.Reducer.min()
        .combine(ee.Reducer.mean(), '', true)
        .combine(ee.Reducer.max(), '', true),
      geometry: roi,
      scale: 30,
      maxPixels: 1e9
    });
  
  // Compute statistics for NDBI
  var ndbi_stats = indices_collection
    .select('NDBI')
    .reduceRegion({
      reducer: ee.Reducer.min()
        .combine(ee.Reducer.mean(), '', true)
        .combine(ee.Reducer.max(), '', true),
      geometry: roi,
      scale: 30,
      maxPixels: 1e9
    });
  
  // Compute statistics for EVI
  var evi_stats = indices_collection
    .select('EVI')
    .reduceRegion({
      reducer: ee.Reducer.min()
        .combine(ee.Reducer.mean(), '', true)
        .combine(ee.Reducer.max(), '', true),
      geometry: roi,
      scale: 30,
      maxPixels: 1e9
    });
  
  // Compute statistics for SAVI
  var savi_stats = indices_collection
    .select('SAVI')
    .reduceRegion({
      reducer: ee.Reducer.min()
        .combine(ee.Reducer.mean(), '', true)
        .combine(ee.Reducer.max(), '', true),
      geometry: roi,
      scale: 30,
      maxPixels: 1e9
    });
  
  // Compute statistics for MNDWI
  var mndwi_stats = indices_collection
    .select('MNDWI')
    .reduceRegion({
      reducer: ee.Reducer.min()
        .combine(ee.Reducer.mean(), '', true)
        .combine(ee.Reducer.max(), '', true),
      geometry: roi,
      scale: 30,
      maxPixels: 1e9
    });
  
  // Compute statistics for NDMI
  var ndmi_stats = indices_collection
    .select('NDMI')
    .reduceRegion({
      reducer: ee.Reducer.min()
        .combine(ee.Reducer.mean(), '', true)
        .combine(ee.Reducer.max(), '', true),
      geometry: roi,
      scale: 30,
      maxPixels: 1e9
    });
  
  // Compute statistics for BSI
  var bsi_stats = indices_collection
    .select('BSI')
    .reduceRegion({
      reducer: ee.Reducer.min()
        .combine(ee.Reducer.mean(), '', true)
        .combine(ee.Reducer.max(), '', true),
      geometry: roi,
      scale: 30,
      maxPixels: 1e9
    });

                      
  // Select the thermal band from the lst_collection
  var thermal = lst_collection.select('B11');
  
  // Select the emissivity band from the indices_collection
  var emissivity = indices_collection.select('Emissivity');
  
  // Calculate Land Surface Temperature (LST) using the formula and assign it to the variable "LST"
  var LST = lst_collection.expression('(Tb / (1 + (11.5 * (Tb / 14388)) * log(Ep))) - 273.15',
    {'Tb': thermal, 'Ep': emissivity}).rename('LST');
  
  // Compute statistics (minimum, mean, maximum) for LST within the selected grid's geometry
  var lst_stats = LST.reduceRegion({
    reducer: ee.Reducer.min().combine(ee.Reducer.mean(), '', true).combine(ee.Reducer.max(), '', true),
    geometry: roi,
    scale: 30,
    maxPixels: 1e9
  });
  
  // Create a dictionary containing the grid and image information, as well as the computed statistics
  var dict = {
    Grid_num: image.get('GridNumber'),
    State_id: image.get('stateID'),
    State_name: image.get('admin_name'),
    City_name: image.get('city_ascii'),

    NDVI_min: ndvi_stats.get('NDVI_min'),
    NDVI_mean: ndvi_stats.get('NDVI_mean'),
    NDVI_max: ndvi_stats.get('NDVI_max'),

    NDBI_min: ndbi_stats.get('NDBI_min'),
    NDBI_mean: ndbi_stats.get('NDBI_mean'),
    NDBI_max: ndbi_stats.get('NDBI_max'),
    
    EVI_min: evi_stats.get('EVI_min'),
    EVI_mean: evi_stats.get('EVI_mean'),
    EVI_max: evi_stats.get('EVI_max'),

    SAVI_min: savi_stats.get('SAVI_min'),
    SAVI_mean: savi_stats.get('SAVI_mean'),
    SAVI_max: savi_stats.get('SAVI_max'),

    MNDWI_min: mndwi_stats.get('MNDWI_min'),
    MNDWI_mean: mndwi_stats.get('MNDWI_mean'),
    MNDWI_max: mndwi_stats.get('MNDWI_max'),

    NDMI_min: ndmi_stats.get('NDMI_min'),
    NDMI_mean: ndmi_stats.get('NDMI_mean'),
    NDMI_max: ndmi_stats.get('NDMI_max'),

    BSI_min: bsi_stats.get('BSI_min'),
    BSI_mean: bsi_stats.get('BSI_mean'),
    BSI_max: bsi_stats.get('BSI_max'),

    LST_min: lst_stats.get('LST_min'),
    LST_mean: lst_stats.get('LST_mean'),
    LST_max: lst_stats.get('LST_max')
  };
  
  // Add the dictionary to the list
  list.push(dict);

    
  // Create a new feature using the dictionary "dict" and assign it to the variable "feature"
  var feature = ee.Feature(null, dict);
  
  // Merge the newly created feature into the finalFeatureCollection
  finalFeatureCollection = finalFeatureCollection.merge(ee.FeatureCollection([feature]));
  
  // Convert the grid number to a string
  var num = ee.String(image.get('GridNumber'));
  
  // Define the export parameters for BSI image
//   var exportParams_bsi = {
//     image: indices_collection.select('BSI'),
//     description: 'bsi_' + num.getInfo(),
//     folder: 'GEE_South_bsi',
//     scale: 30,
//     region: image.geometry().bounds(),
//     maxPixels: 1e9,
//     formatOptions: { cloudOptimized: true }
//   };
  
  // Define the export parameters for LST image
  var exportParams_lst = {
    image: LST,
    description: 'lst_' + num.getInfo(),
    folder: 'GEE_South_lst',
    scale: 30,
    region: image.geometry().bounds(),
    maxPixels: 1e9,
    formatOptions: { cloudOptimized: true }
  };
  
  // Export the image to Google Drive
  Export.image.toDrive(exportParams_bsi);
}

print(list);