// // DELIMITACIÓN DE COBERTURA GLACIAR A PARTIR DE NDSI MEDIANTE IMÁGENES LANDSAT
// // Por: Hector Hugo Añamuro Luque -- Email: hugo.aluque@gmail.com
/** ######################################################################################################
###########################################################################################################*/
/** Establecer el periodo de tiempo para los promedios de imagenes las imagenes 
            === TENER EN CONSIDERACION QUE EL RANGO DEL PRODUCTO DE CUERPOS DE AGUA (1985 - 2019)==*/


var año = 2019;  //Para Producto de Cuerpos de agua
var nubes = 5; //Colocar el porcentaje de nubes deseado para filtrar la imagen
var U_NDSI = 0.45; // Colocar el Umbral de NDSI para un mejor delimitación
// var inicio = '1985-03-01';
// var fin = '1985-10-30';

/** Establecer la plataforma con la que trabajar ya sea Landsat 4-5 (L5) o Landsat 8 (L8)
    -Teniendo en cuenco el tiempo en la que cada uno de las misiones estuvieron activas
                Para Landsat 4 y 5 usar  == "L5"
                Para Lansat 8 usar == "L8"*/
                
var platforma = 'L8'; 

/** ######################################################################################################
###########################################################################################################*/

// FUNSION PARA CLACULO DE NDSI
var ndsi = function (image){
    var NDSI = image.expression('(green-swir1)/(green+swir1)', {
      'green':image.select('green'),
      'swir1':image.select('swir1')
    }).rename('ndsi')
    return image.addBands(NDSI)
}

// FUNSION PARA MASCARA DE NUBES L5 Y L7
var cloudMaskL457 = function(image) {
  var qa = image.select('pixel_qa');
  // If the cloud bit (5) is set and the cloud confidence (7) is high
  // or the cloud shadow bit is set (3), then it's a bad pixel.
  var cloud = qa.bitwiseAnd(1 << 5)
                  .and(qa.bitwiseAnd(1 << 7))
                  .or(qa.bitwiseAnd(1 << 3));
  // Remove edge pixels that don't occur in all bands
  var mask2 = image.mask().reduce(ee.Reducer.min());
  return image.updateMask(cloud.not()).updateMask(mask2);
};

// FUNSIONES PARA MASCARA DE NUBES L8
function maskL8sr(image) {
  // Bits 3 and 5 are cloud shadow and cloud, respectively.
  var cloudShadowBitMask = (1 << 3);
  var cloudsBitMask = (1 << 5);
  // Get the pixel QA band.
  var qa = image.select('pixel_qa');
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
                 .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return image.updateMask(mask);
}

// selecion de la plataforma para realizar el analisis 

if (platforma == 'L5' | platforma == 'L5') {
  var ImCol = 'LANDSAT/LT05/C01/T1_SR';
  var pl = 'Landsat 5';
} else {
  var ImCol = 'LANDSAT/LC08/C01/T1_SR';
  var pl = 'Landsat 8';
}                         
// MASCARA DE AGUA A PARTIR DE PRODUCTO: "JRC Yearly Water Classification History, v1.2"
// CARGANDO PRODUCTO DE CUERPOS DE AGUA ANUAL
var C_agua = ee.ImageCollection('JRC/GSW1_2/YearlyHistory')
            .filterBounds(Area_de_Interes)
            // .filterDate('2020-10-01','2020-10-30')
            .filter(ee.Filter.calendarRange(año,año, 'year'))
            // .select('waterClass');

var CA_mask = C_agua.mosaic().clip(Area_de_Interes)
var Mascara = ee.Image(0).where(CA_mask.lt(0),0)
                        .where(CA_mask.gte(0),0)
                        .where(CA_mask.gte(2),1);
// Map.addLayer(Mascara,{},'Mascara de Agua')
// // CARGAR COLECCIÓN LANDSAT 
// var Lcoll = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR")
//             .filterBounds(Area_de_Interes)
//             // .filterDate('2020-05-01','2020-09-30')
//             .filter(ee.Filter.calendarRange(2019, 2019, 'year'))
//             .filterMetadata('CLOUD_COVER', 'less_than', 10);
           

// Cargando la coleccion de imagenes de la plataforma seleccionada
var Landsat_Coll = ee.ImageCollection(ImCol) // aqui se establece la coleccion segun la plataforma 
                  .filterBounds(Area_de_Interes)
                  // .filterDate(inicio,fin)
                  .filter(ee.Filter.calendarRange(año,año, 'year'))
                  .filterMetadata('CLOUD_COVER', 'less_than', nubes)


// Aplicando la mascara de nubes conforme a la plataforma seleccionada
if (platforma == 'L5' | platforma == 'L5') {
  var L_coll_m = Landsat_Coll.map(cloudMaskL457);
  } else {
  var L_coll_m = Landsat_Coll.map(maskL8sr);
  }

// Aplicando un mosaico a la coleccion  ?==========================================================================
var Land_coll_Media = Landsat_Coll.median().clip(Area_de_Interes)


// SELECCIONAR LA BANDA CON A CUAL REALIZAR LA MASCARA
var datamask = Mascara.select('waterClass');

// CREAR UNA MASCARA BINARIA
var mask = Mascara.eq(0);

// ACTUALIZAR LA COLECCION CON LA MASRCA DE AGUA
var LColl_M = Land_coll_Media.updateMask(mask);



// print(Land_coll_Media)


// Calculando NDSI para la plataforma seleccionada
if (platforma == 'L5' | platforma == 'L5') {
  var NDSI = LColl_M.normalizedDifference(['B2', 'B5']);
} else {
  var NDSI = LColl_M.normalizedDifference(['B3', 'B6']);
}


// Aplicando umbrales para determinar zonas con Glaciar y sin Glaciar

var Cob_Glaciar = NDSI.gte(U_NDSI);


Map.addLayer(Land_coll_Media,{bands:["B6","B5","B4"],gamma: 1.5,max:7000, min:150},'Landsat 8',false)
Map.addLayer(Land_coll_Media,{bands:["B5","B4","B3"],gamma: 1.5,max:7000, min:150},'Landsat 5',false)
Map.addLayer(NDSI,{min:0,max:1,palette:['#FFFFFF','#11CDEE']},'NDSI',false)
Map.addLayer(Cob_Glaciar,{min:0,max:1,palette:['#FFFFFF','#11CDEE']},'Cob_Glaciar',true)
Map.centerObject(Area_de_Interes)

// Exportando Resultado

var link = Cob_Glaciar.getDownloadURL({
          name:'Cob_Glaciar',
          region:Area_de_Interes,
          scale:30,
          crs:'EPSG:32719',
          maxPixels:1e13,
          fileFormat: 'GeoTIFF'
});

print('Link para la Descarga:')
print(link)