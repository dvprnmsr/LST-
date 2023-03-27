# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 08:50:59 2023

@author: purnamas
"""

import ee
import math

try:
    ee.Initialize()
    
except Exception as error:
    ee.Authenticate()
    ee.Initialize()

class GEEImage:
    # conservative values for wet and dry soil albedo. Those values are adapted from AWRA (van Dijk, 2013)
    A_WET = ee.Image.constant(0.16) # Soil albedo for wet conditions [-]
    A_DRY = ee.Image.constant(0.26) # Soil albedo for dry conditions [-]
    VEG_CAP = ee.Image.constant(0.65) # Vegetation capacity [-]
    W_0_REF = ee.Image.constant(0.30) # Reference soil water content [-]
    SB_CONST = ee.Image.constant(5.67e-8) # Stefan-Boltzmann constant [W m-2 K-4]
    RHO_AIR = ee.Image.constant(1.293) # Air density [kg/m3]
    C_P = ee.Image.constant(1013) # Specific heat capacity of air [J/kg.K]
    w0 = ee.Image.constant(0.5) # Soil moisture [-], this should be from wflow_sbm
    h = ee.Image.constant(1) # Vegetation canopy which can derived from other map

    
    def __init__(self, image):
        self.image = image
 
    
    def from_bands(self):
        """Create a new GEEImage instance from selected bands"""
        bands = [
            'surface_pressure',
            'temperature_2m', 
            'u_component_of_wind_10m',
            'v_component_of_wind_10m',
            'surface_solar_radiation_downwards_sum', #without '***_sum'if it is ERA5 land hourly
            'surface_thermal_radiation_downwards_sum',
            'surface_net_thermal_radiation_sum',
            'surface_net_solar_radiation_sum',
            'total_evaporation_sum',
            'leaf_area_index_high_vegetation',
            'leaf_area_index_low_vegetation',
        ]
        
        selected_bands = self.image.select(bands)
        return GEEImage(selected_bands)


    def preprocessing(self):
        """Apply preprocessing functions to the image."""
        def wind_speed(image):
            """Add 2m wind speed band to the image."""
            u = image.select('u_component_of_wind_10m')
            v = image.select('v_component_of_wind_10m')
            wind_10m = (u.pow(2).add(v.pow(2))).pow(0.5).rename('wind_10m')
            wind_2m = wind_10m.multiply(math.pow((2 / 10), 0.143)).rename('wind_2m')
            return image.addBands(wind_2m)
        
        def conversion(image): 
            """Convert solar radiation from J/m2/day to W/m2."""
            Rs_in = image.select('surface_solar_radiation_downwards_sum').divide(86400).rename('Rs_in')
            return image.addBands(Rs_in)
        
        def fveg(image):
            """Add fractional vegetation cover [-] band to the image."""
            lai_high = image.select('leaf_area_index_high_vegetation')
            lai_low = image.select('leaf_area_index_low_vegetation')
            lai_tot = lai_high.add(lai_low).rename('lai')
            lai_ref = ee.Image(4) #find an alternative to this 
            fveg = ee.Image(1).subtract(lai_tot.divide(lai_ref).multiply(ee.Image(-1).exp())).rename('fveg')
            return image.addBands(fveg)
        
        def surf_albedo(image):
            """Surface albedo is calculated from vegetation fraction (fveg)
            and relative soil moisture content of the top soil layer"""
            fveg = image.select('fveg')
            fsol = ee.Image(1).subtract(fveg)
            aveg = ee.Image(0.452).multiply(self.VEG_CAP)
            asol = self.A_WET.add(self.A_DRY.subtract(self.A_WET).multiply(self.w0.divide(self.W_0_REF).multiply(ee.Image(-1)).exp()))
            surf_albedo = fveg.multiply(aveg).add(fsol.multiply(asol)).rename('surf_albedo')
            return image.addBands(surf_albedo)   
        
        def Rs_out(image): #Outgoing shortwave radiation [W m-2]
            """Daytime outgoing/upwelling shortwave radiation which is prescribed 
            by the ratio of reflected radiation by surface to incoming radiation 
            as prescribed in albedo values """
            Rs_out = image.select('surf_albedo').multiply(image.select('Rs_in')).rename('Rs_out')
            return image.addBands(Rs_out)
        
        def Rl_out(image):#Outgoing longwave radiation [W m-2]
            """Outgoing/upward longwave radiation calculates radiation emitted by earth 
            surface. The magnitude is governed by the Stefan-Boltzmann law with e = 1"""
            Rl_out = image.select('temperature_2m').pow(4).multiply(self.SB_CONST).rename('Rl_out')
            return image.addBands(Rl_out)
        
        def Rl_in(image):
            """Incoming longwave radiation calculates radiation received by earth 
            surface"""
            surface_pressure = image.select('surface_pressure')
            temp = image.select('temperature_2m')
            f_atm = ee.Image.constant(0.65).multiply((surface_pressure.divide(temp)).pow(0.14)) #Brutsaert (1975)
            Rl_in = f_atm.multiply(self.SB_CONST).multiply(temp.pow(4)).rename('Rl_in')
            return image.addBands(Rl_in)
        
        def aero_cond(image): 
            """Aerodynamic conductance"""
            f_h = ee.Image.constant(813).divide(self.h).subtract(5.45).log() #function of vegetation canopy [-]
            k = ee.Image.constant(0.305).divide(f_h.multiply(f_h.add(ee.Image.constant(2.3)))) #roughness coefficient for wind measurements at 2m [-] 
            ga = k.multiply(image.select('wind_2m')).rename('ga')
            return image.addBands(ga)
                
        def combined_function(image):
            result = wind_speed(image)
            result = conversion(result)
            result = fveg(result)
            result = surf_albedo(result)
            result = Rs_out(result)
            result = Rl_out(result)
            result = Rl_in(result)
            result = aero_cond(result)
            return result
            
        #Apply the processing to all images in the collection using the map metho
        return self.image.map(combined_function)
    
    def postprocessing(self):
        """Apply preprocessing functions to the image."""
        def latent_heat(image):
            LE = image.select('actevap').multiply(1000).multiply(2.5e6).divide(1000).divide(3600).rename('LE') #convert from mm/day to W/m2
            Rn = image.select('Rs_in').subtract(image.select('Rs_out')).add(image.select('Rl_in')).subtract(image.select('Rl_out')).rename('Rn')
            # Rn_resampled = Rn('bicubic').reproject(wflow_input.projection, None, res_LE).rename('Rn_resampled')
            H = Rn.subtract(LE).rename('H')
            # Ts = image.select('temperature_2m').add(H.divide(self.RHO_AIR.multiply(self.C_P).multiply(image.select('ga')))).rename('Ts')
            Ts = image.select('temperature_2m').add(H.divide(self.RHO_AIR.multiply(self.C_P).multiply(0.02))).rename('Ts')
            return image.addBands(Ts)
        
        return self.image.map(latent_heat)
    
    
class ImageResampler:
    def __init__(self, reference_image):
        self.reference_image = reference_image
        
        
    def resample(self, image, method = 'bicubic'):
        scale = self.reference_image.projection().nominalScale() 
        resampled = image.resample(method).reproject(self.reference_image.projection(), None, scale)
        return resampled
        
    
    # def EB(self):
    #     def Ts(image):
    #         """Latent heat flux is computed from evapotranspiration"""
    #         LE = image.select('actevap').multiply(1000).multiply(2.5e6).divide(1000).divide(3600).rename('LE') #convert from mm/day to W/m2
    #         Rn = image.select('Rs_in').subtract(image.select('Rs_out')).add(image.select('Rl_in')).subtract(image.select('Rl_out')).rename('Rn')
    #         Rn_resampled = Rn('bicubic').reproject(wflow_input.projection, null, res_LE).rename('Rn_resampled')
    #         H = Rn_resampled.subtract(LE).rename('H')
    #         Ts = image.select('temperature_2m').add(H.divide(self.RHO_AIR.multiply(self.C_P).multiply(image.select('ga')))).rename('Ts')
    #         return image.addBands(Ts)
        
    #     return self.image.map(Ts)
        
    

    
    
        # def add_wind_speed
        # u = ee.Image(self.image).select('u_component_of_wind_10m')
        # v = ee.Image(self.image).select('v_component_of_wind_10m') 
        # wind_10m = (u.pow(2).add(v.pow(2))).pow(0.5).rename('wind_10m') 
        # wind_2m = wind_10m.multiply(math.pow((2 / 10), 0.143)).rename('wind_2m')
        # preprocessing = self.image.addBands(wind_2m)
        # return ee.Image(self.image.addBands(wind_2m))
        
        
        # self.bands = ee.List([
        #     'surface_pressure',
        #     'temperature_2m',
        #     'u_component_of_wind_10m',
        #     'v_component_of_wind_10m',
        #     'surface_solar_radiation_downwards_hourly',
        #     'surface_thermal_radiation_downwards_hourly',
        #     'leaf_area_index_high_vegetation',
        #     'leaf_area_index_low_vegetation',
        # ])
        
    # def preprocessing(self):
    #    def add_wind_speed(image):
    #        u = image.select('u_component_of_wind_10m')
    #        v = image.select('v_component_of_wind_10m')
    #        wind_10m = (u.pow(2).add(v.pow(2))).pow(0.5).rename('wind_10m')
    #        wind_2m = wind_10m.multiply(math.pow((2 / 10), 0.143)).rename('wind_2m')
    #        return image.addBands(wind_2m)

    #    # Create an ImageCollection from the image
    #    collection = ee.ImageCollection.fromImages(self.image)

    #    # Apply the processing to all images in the collection using the map method
    #    collection_with_wind_speed = collection.map(add_wind_speed)

    #    return collection_with_wind_speed.select(self.bands)
    
    
     
    # def __init__(self, w0, Rs_in, Ta, surface_pressure, wind_2m, h, lai, actevap):
    #     self.w0 = w0 #relative soil moisture content [-]
    #     self.Rs_in = Rs_in #incoming shortwave radiation [W m-2]
    #     self.Ta = Ta #effective average air temperature [degree Celcius]
    #     self.surface_pressure = surface_pressure #atmospheric vapor pressure [Pa] or in short we can call it pe
    #     self.wind_2m = wind_2m #effective wind speed at 2 m [m/s]
    #     self.h = h #height of vegetation canopy [m]
    #     self.lai = lai
    #     self.lai_ref = 1
    #     self.actevap = actevap #actual evapotranspiration from wflow_sbm (actevap [mm/timestep])
        #
#     def fveg(self):
#         """Fractional canopy cover [-]"""
#         self.fveg = 1 - math.exp(-(self.lai/ self.lai_ref))
#         return self.fveg
    
#     def surf_albedo(self):
#         """Surface albedo is calculated from vegetation fraction (fveg)
#         and relative soil moisture content of the top soil layer"""
#         self.fsol = 1- self.fveg
#         self.aveg = 0.452 * self.VEG_CAP
#         self.asol = self.A_WET + (self.A_DRY - self.A_WET) * math.exp(-(self.w0 / self.W_0_REF))
#         self.surf_albedo = self.fveg * self.aveg + self.fsol * self.asol
#         return self.surf_albedo
    
#     def Rs_out(self): #Outgoing shortwave radiation [W m-2]
#         """Daytime outgoing/upwelling shortwave radiation which is prescribed 
#         by the ratio of reflected radiation by surface to incoming radiation 
#         as prescribed in albedo values """
#         self.Rs_out = self.surf_albedo * self.Rs_in
#         return self.Rs_out
        
#     def Rl_in(self):
#         """Incoming longwave radiation calculates radiation received by earth 
#         surface"""
#         self.f_atm = 0.65 * (self.surface_pressure / (self.Ta + 273.15))**0.14 #Brutsaert (1975)
#         self.Rl_in = self.f_atm * self.SB_CONST * (self.Ta + 273.15)**4
#         return self.Rl_in 
    
#     def Rl_out(self):#Outgoing longwave radiation [W m-2]
#         """Outgoing/upward longwave radiation calculates radiation emitted by earth 
#         surface. The magnitude is governed by the Stefan-Boltzmann law with e = 1"""
#         self.Rl_out = self.SB_CONST * (self.Ta + 273.15)**4
#         return self.Rl_out 
    
#     def aero_cond(self): 
#         """Aerodynamic conductance"""
#         self.f_h = math.log((813 / self.h) - 5.45) #function of vegetation canopy [-]
#         self.k = 0.305 / (self.f_h * (self.f_h + 2.3)) #roughness coefficient for wind measurements at 2m [-] 
#         self.ga = self.k * self.wind_2m
#         return self.ga
    
#     def Ts(self):
#         """Latent heat flux is computed from evapotranspiration"""
#         self.LE = self.actevap * 1000 * 2.5e6 / 1000 /86400
#         self.Rn = self.Rs_in - self.Rs_out  + self.Rl_in - self.Rl_out
#         self.H = self.Rn - self.LE
#         self.Ts = self.Ta + (self.H / (self.RHO_AIR * self.C_P * self.ga))
#         return self.Ts
        

        
# #     def G(self, fveg, Rs_in, Fgr_max, Fs_ref): #ground heat flux
# #         """Ground heat flux [W/m2]"""
# #         self.fveg = fveg
# #         self.fsol = 1 - self.fveg
# #         self.Rs_in = Rs_in
# #         self.Fgr_max = Fgr_max
# #         self.Fs_ref = Fs_ref
# #         self.Fgr = self.Fgr_max * (1 - math.exp(- self.fsol/self.Fs_ref))
# #         self.g_a = self.Fgr * Rs_in 
        

         
# # SEB = trial.SEB()
# SEB =  EnergyBalance(0.8, 167.8, 25, 25, 2, 2, 0.9, 3)                                                                                                                  
# fveg = SEB.fveg()
# alb = SEB.surf_albedo()
# Rs_out = SEB.Rs_out()
# Rl_in = SEB.Rl_in()
# Rl_out = SEB.Rl_out()
# ga = SEB.aero_cond()
# Ts = SEB.Ts()
      