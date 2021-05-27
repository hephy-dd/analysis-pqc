#!/usr/bin/env python3


import analysis_pqc as fp
import pqc_analysis_tools as tools
import pqc_analysis_json as pqc
import os
import sys
import argparse
import pandas as pd
import yaml
import math
import matplotlib.pyplot as plt


def df_to_csv(df, flute):

    f = '~/{}.csv'.format(flute)
    hdr = False  if os.path.isfile(f) else True
    

    df.to_csv(f, mode='a', header=hdr, index=False)



def read_config():
    
    with open('PQC_specs.yml', 'r') as f:
         conf = yaml.load(f, Loader=yaml.FullLoader)
    
    
    return conf



def append_dictionary(current_dict, dict_to_append, flute, test):
    
    conf = read_config()
    
    for key, values in dict_to_append.items():
        if key in current_dict.keys():
              if not math.isnan(float(values)) and abs(values-conf[flute][test]) < abs(current_dict[key]-conf[flute][test]):
                   current_dict.update({'{}'.format(key): values})
                
              else:
                  pass
        else:
            current_dict.update({'{}'.format(key): values})

    return current_dict  


def append_parameter(data):

    parameter_list =[]
    parameter_list.append(data)
    return parameter_list

##############################################################################################################################################################


class Flute1:

    
    def __init__(self, path, label):
        
        self.path = path
        self.label = label
        self.flute = 'Flute 1'
        
        self.functions = {'mos': self.mos,
              'fet': self.fet,
              'van-der-pauw': self.vdp,
              'capacitor': self.capacitor}
    

 



    def mos(self, datafile, data_dict):

       try:
           v_fb1, v_fb2, t_ox, n_ox, Q_ox, c_acc = pqc.analyse_mos_data(datafile, False)
        
       except Exception as err:
           v_fb2, c_acc, n_ox = float('NaN')

       par_dict = ({'Vfb [V]': round(v_fb2, 2) if not math.isnan(v_fb2) else v_fb2,
              'Cacc [pF]': round(c_acc*1e12, 2) if not math.isnan(c_acc) else c_acc,
              'n_ox [cm^-2]': n_ox})

       data_dict = append_dictionary(data_dict, par_dict, self.flute, 'mos_vfb')
       
       return data_dict




    def fet(self, datafile, data_dict):

       try:
           v_th = pqc.analyse_fet_data(datafile, False)
       except:
           v_th = float('NaN') 
    
    
       par_dict = ({'V_th [V]' : round(v_th, 2) if not math.isnan(v_th) else v_th})
  
       data_dict = append_dictionary(data_dict, par_dict, self.flute, 'fet_vth')
    
       return data_dict




    def vdp(self, datafile, data_dict, linewidth=False):
 
       try:
           r_sheet, lbl = pqc.analyse_van_der_pauw_data(datafile, False)
           

       except Exception as err:
           r_sheet = float('NaN')

       if "N" in datafile:
           test = 'N_vdp'
           if "Reverse" in datafile:
              par_dict = ({'VdP-N (r) [Ohm\sq]': round(r_sheet, 2) if not math.isnan(r_sheet) else r_sheet})
           else:
              par_dict = ({'VdP-N (s) [Ohm\sq]': round(r_sheet, 2) if not math.isnan(r_sheet) else r_sheet})    
   
       elif "P-stop" in datafile:
              test = 'P-stop_vdp'
              if "Reverse" in datafile:
                   par_dict = ({'VdP-Pstop (r) [kOhm\sq]' : round((r_sheet*1e-3) ,2) if not math.isnan(r_sheet) else r_sheet})
              else:
                   par_dict = ({'VdP-Pstop (s) [kOhm\sq]' : round((r_sheet*1e-3),2) if not math.isnan(r_sheet) else r_sheet}) 
   
       else:
             test = 'Rpoly_vdp'
             if "Reverse" in datafile:
                  par_dict = ({'VdP-Rpoly (r) [kOhm\sq]' : round((r_sheet*1e-3) ,2) if not math.isnan(r_sheet) else r_sheet}) 
             else:
                  par_dict = ({'VdP-Rpoly (s) [kOhm\sq]' : round((r_sheet*1e-3),2) if not math.isnan(r_sheet) else r_sheet}) 

       
          
       if linewidth:
       
             return r_sheet

       else: 
             
            data_dict = append_dictionary(data_dict, par_dict, self.flute, test) 
            return data_dict




    def capacitor(self, datafile, data_dict):

    
       try:  
           c_mean, c_median = pqc.analyse_capacitor_data(datafile, False)


       except Exception as err:
           c_median = float('NaN')


       if "_Left_" in datafile:
           par_dict = ({ 'CapW [pF]' : round(c_median*1e12, 2) if not math.isnan(c_median) else c_median})
       else:
             par_dict = ({'CapE [pF]' : round(c_median*1e12, 2) if not math.isnan(c_median) else c_median}) 

       data_dict = append_dictionary(data_dict, par_dict, self.flute, 'cap') 


       return data_dict

    
   
   
   
    def run(self):

      
      data_dict = {'Halfmoon': self.label} 
      for  key, f in self.functions.items():

        directory = tools.find_all_files_from_path(self.path, '{}'.format(key))
        for i in directory:
           if 'Flute_1' in i and 'eft' in i:
               data_dict = f(i, data_dict)
    
      if os.path.isfile('~/{}.csv'.format(self.flute)):
           df = pd.DataFrame(data=data_dict, columns= None, index=[0])
      else:
           df = pd.DataFrame(data=data_dict, columns = data_dict.keys(), index=[0])
    
      ####re-arrange the order of the parameters
      cols_order = ['Halfmoon', 'Vfb [V]', 'Cacc [pF]', 'n_ox [cm^-2]', 'V_th [V]', 'VdP-N (s) [Ohm\sq]', 'VdP-N (r) [Ohm\sq]', 'VdP-Pstop (s) [kOhm\sq]', 'VdP-Pstop (r) [kOhm\sq]', 'VdP-Rpoly (s) [kOhm\sq]', 'VdP-Rpoly (r) [kOhm\sq]', 'CapE [pF]']

      cols = df.columns.tolist()
      cols = cols_order
      df = df[cols]
      df_to_csv(df, self.flute)
      
          

#####################################################################################################################


class Flute2:
    
    def __init__(self, path, label):
        
        self.path = path
        self.label = label
        self.flute = 'Flute 2'
        
        self.functions = {'gcd': self.gcd,
              'meander': self.poly_meander, 
              'linewidth': self.linewidth,
              'breakdown': self.breakdown}
    
    
        
    def gcd(self, datafile, data_dict):
      
        try:
            i_surf, i_bulk, vfb_acc, vfb_inv, s0 = pqc.analyse_gcd_data(datafile, False)
        except Exception as err:
            i_surf=i_bulk=vfb_acc=vfb_inv=s0 = float("NaN")
      

        par_dict = ({'I_surf [pA]': round(i_surf,2) if not math.isnan(i_surf) else i_surf,
              'S_0 [cm/s] ':  round(s0,2) if not math.isnan(s0) else s0,
              'Vth_acc [V]': round(vfb_acc+5,2)  if not math.isnan(vfb_acc) else vfb_acc,
              'Vth_inv [V]': round(vfb_inv, 2) if not math.isnan(vfb_inv) else vfb_inv})

        data_dict = append_dictionary(data_dict, par_dict, self.flute, 'gcd_vfb')
       
        return data_dict


   
    def poly_meander(self, datafile, data_dict):

        try:
            rho_sq = pqc.analyse_meander_data(datafile)
            rho_sq = round(rho_sq*1e-6, 2)
        except:
            rho_sq = float('NaN')
            
        par_dict = ({'Poly_meander [MOhm]': rho_sq})
         
        data_dict = append_dictionary(data_dict, par_dict, self.flute, 'Poly_meander')

        return data_dict

  
  
   
    def linewidth(self, datafile, data_dict):
       

            vdp_list = tools.find_all_files_from_path(self.path, 'van-der-pauw') 
            
            if "_N_" in datafile:
              
              test = 'N_linewidth'
              try:
                for i in vdp_list:
                  if 'N' in i and 'Reverse' not in i: 

                    r_sheet = Flute1.vdp(self, i, data_dict, True)
                    
                    t_line = pqc.analyse_linewidth_data(datafile, r_sheet, False)
                    t_line = round(t_line,2)

              except:
                   t_line = float("NaN")
             
              par_dict = ({'N-line [um]': t_line})    
   
            elif "P-stop" in datafile:
               test = 'P-stop linewidth'
               try:
                for i in vdp_list:
                  if 'P-stop' in i and 'Reverse' not in i:
                    r_sheet = Flute1.vdp(self, i, data_dict, True)
                    t_line = pqc.analyse_linewidth_data(datafile, r_sheet, False)
                    t_line = round(t_line, 2)                   

               except:
                    t_line = float("NaN")
             
               if "2-wire" in datafile:      
                  par_dict = ({'Pstop-line [um] (2wire)' : t_line})     
               elif "4-wire" in datafile:
                  par_dict = ({'Pstop-line [um] (4wire)' : t_line})
            
            if par_dict: 
              data_dict = append_dictionary(data_dict, par_dict, self.flute, test)

          
              return data_dict



    
    def breakdown(self, datafile, data_dict):

        try:
            v_bd = pqc.analyse_breakdown_data(datafile, False)
            v_bd = round(v_bd, 2)

        except Exception as err:
            v_bd = float('NaN')

        par_dict = ({'V_bd [V]' : v_bd})

        data_dict = append_dictionary(data_dict, par_dict, self.flute, 'breakdown')

        return data_dict




    def run(self):

      
      data_dict = {'Halfmoon': self.label} 

      for  key, f in self.functions.items():

        directory = tools.find_all_files_from_path(self.path, '{}'.format(key))
        for i in directory:
           if 'Flute_2' in i and 'eft' in i:
               
               data_dict = f(i, data_dict)
    
      if os.path.isfile(os.getcwd() + '~/{}.csv'.format(self.flute) ):
           df = pd.DataFrame(data=data_dict, columns = None, index=[0])
      else:
           df = pd.DataFrame(data=data_dict, columns = data_dict.keys(), index=[0])
    
      df_to_csv(df, self.flute)




#####################################################################################################################

class Flute3:    


 def __init__(self, path, label):
        
        self.path = path
        self.label = label
        self.flute = 'Flute 3'
        
        self.functions = {'bulk': self.bulk_cross,
              'clover': self.clover_metal,  
              'meander': self.metal_meander,
              'van-der-pauw': self.pPlus_vdp,
              'linewidth': self.pPlus_linewidth,
              'cv': self.cv,
              'iv': self.iv}
    


 def bulk_cross(self, datafile, data_dict):

      try: 
          r_sheet, _ = pqc.analyse_cross_data(datafile)
          r_sheet = round(r_sheet*1e-3, 2)
      except:
          r_sheet = float('NaN')
      
      if 'Reverse' in datafile:
            par_dict = ({'vdP_Bulk (r) [kOhm/sq]': r_sheet})
      else:
            par_dict = ({'vdP_Bulk (s) [kOhm/sq]': r_sheet})

      data_dict = append_dictionary(data_dict, par_dict, self.flute, 'bulk_cross')

      return data_dict




 def clover_metal(self, datafile, data_dict):

      try: 
          rho_sq = pqc.analyse_meander_data(datafile)
          rho_sq = round(rho_sq*1e3, 2)
      except:
          rho_sq = float('NaN')
      
      if 'Reverse' in datafile:
             par_dict = ({'CloverMetal (r) [mOhm/sq]': rho_sq})
      else: 
             par_dict = ({'CloverMetal (s) [mOhm/sq]': rho_sq})
  
      data_dict = append_dictionary(data_dict, par_dict, self.flute, 'clover_metal')

      return data_dict

 
 
 def metal_meander(self, datafile, data_dict):

        try:
            rho_sq = pqc.analyse_meander_data(datafile)
            rho_sq - round(rho_sq, 2)
        except:
            rho_sq = float('NaN')

        par_dict = ({'R_meander [Ohm]': rho_sq})

        data_dict = append_dictionary(data_dict, par_dict, self.flute, 'metal_meander')

        return data_dict

 
  
  
 def pPlus_vdp(self, datafile, data_dict, line = False):

   if 'P_cross' in datafile:
    
      try: 
          r_sheet, _ = pqc.analyse_van_der_pauw_data(datafile, False)
          r_sheet = round(r_sheet*1e-3, 2)
      except:
          r_sheet = float('NaN')
     
    
      if 'Reverse' in datafile:
            par_dict = ({'pPlus_vdp (r) [kOhm/sq]': r_sheet})
      else:
            par_dict = ({'pPlus_vdp (s) [kOhm/sq]': r_sheet})

    
      if line:
          return r_sheet
      else:
          data_dict = append_dictionary(data_dict, par_dict, self.flute, 'pPlus_vdp')
          return data_dict
   
   else:
        return data_dict


 def pPlus_linewidth(self, datafile, data_dict):

    try:
         vdp_list = tools.find_all_files_from_path(self.path, 'van-der-pauw') 
         for i in vdp_list:
             if 'P_cross' in i: 
                r_sheet = self.pPlus_vdp(i, data_dict, line=True) 
            
                t_line = pqc.analyse_linewidth_data(datafile, r_sheet, False)
      
    except Exception as err:
                 t_line = float('NaN')
    
    par_dict = ({'Linewidth pPlus [um]': t_line})

    data_dict = append_dictionary(data_dict, par_dict, self.flute, 'pPlus_linewidth')

    return data_dict
    
       

 def cv(self, datafile, data_dict):

      try: 
          vfd = pqc.analyse_cv_data(datafile, False)
      except Exception as err:
          vfd = float('NaN')

      par_dict = ({' V_fd [V]': vfd})

      data_dict = append_dictionary(data_dict, par_dict, self.flute, 'cv')

      return data_dict


 def iv(self, datafile, data_dict):

      try: 
          i600, _ = pqc.analyse_iv_data(datafile, False)
      except Exception as err:
          i600 = float('NaN')

      par_dict = ({' I600 [A]': i600})

      data_dict = append_dictionary(data_dict, par_dict, self.flute, 'iv')

      return data_dict




 def run(self):

      
      data_dict = {'Halfmoon': self.label} 

      for  key, f in self.functions.items():
 
        directory = tools.find_all_files_from_path(self.path, '{}'.format(key))
        for i in directory:
           if 'Flute_3' in i: 
              data_dict = f(i, data_dict)
      if os.path.isfile('~/{}.csv'.format(self.flute) ):
           df = pd.DataFrame(data=data_dict, columns = None, index=[0])
      else:
           df = pd.DataFrame(data=data_dict, columns = data_dict.keys(), index=[0])
    
      df_to_csv(df, self.flute)





#####################################################################################################################
class Flute4:


 def __init__(self, path, label):
        
        self.path = path
        self.label = label
        self.flute = 'Flute 4'
        
        self.functions = {'gcd05': self.gcd,
              'cbkr': self.cbkr,  
              'contact': self.contact,
              }
 

 

 def gcd(self, datafile, data_dict):
      
        try:
            i_surf, i_bulk, vfb_acc, vfb_inv, s0 = pqc.analyse_gcd_data(datafile, False)
            
        except Exception as err:
            i_surf=i_bulk=vfb_acc=vfb_inv=s0 = float('NaN')

        par_dict = ({'I_surf [pA]': round(i_surf,2) if not math.isnan(i_surf) else i_surf,
              'S_0 [cm/s] ': round(s0,2) if not math.isnan(s0) else s0,
              'Vfb_acc [V]': round(vfb_acc+5,2) if not math.isnan(vfb_acc) else vfb_acc,
              'Vfb_inv [V]': round(vfb_inv,2) if not math.isnan(vfb_inv) else vfb_inv})

        data_dict = append_dictionary(data_dict, par_dict, self.flute, 'gcd_vfb')
       
        return data_dict




 def cbkr(self, datafile, data_dict):

        try:
            vdp_list = tools.find_all_files_from_path(self.path, 'van-der-pauw') 
            for i in vdp_list:
              if 'N_' in i and not 'Reverse' in i: 
                if 'N_' in datafile: 
                 r_sheet = Flute1.vdp(self, i, data_dict, True)
                 r = pqc.analyse_cbkr_data(datafile, r_sheet, False)
              elif 'Polysilicon' in i and not 'Reverse' in i:
                if 'Polysilicon' in datafile:
                 r_sheet = Flute1.vdp(self, i, data_dict, True)
                 r = pqc.analyse_cbkr_data(datafile, r_sheet, False)

        except Exception as err:
            r = float('NaN')

        if 'N' in datafile:
            par_dict = ({'CBKR_strip [Ohm/sq]' : round(r, 2) })
            test = 'cbkr_N'
        else:
            par_dict = ({ 'CBKR_poly [kOhm/sq]' : round(r*1e-3, 2)})
            test = 'cbkr_poly'

        data_dict = append_dictionary(data_dict, par_dict, self.flute, test)

        return data_dict


 def contact(self, datafile, data_dict):
         
        try:
            r_contact = pqc.analyse_contact_data(datafile)
             
        except Exception as err:
            r_contact = float('NaN')

        if '_P_' in datafile:
            par_dict = ({'CC_Edge [kOhm]' : round(r_contact*1e-3, 2) if not math.isnan(r_contact) else r_contact})
            test = 'cc_pPlus'
        elif '_N_' in datafile:
             par_dict = ({'CC_Strip [kOhm]' : round(r_contact*1e-3, 2) if not math.isnan(r_contact) else r_contact })
             test = 'cc_N'
        else:
            par_dict = ({ 'CC_poly [MOhm]' : round(r_contact*1e-6, 2) if not math.isnan(r_contact) else r_contact})
            test = 'cc_poly'

        data_dict = append_dictionary(data_dict, par_dict, self.flute, test)
        return data_dict




 def run(self):

      
      data_dict = {'Halfmoon': self.label} 

      for  key, f in self.functions.items():
 
        directory = tools.find_all_files_from_path(self.path, '{}'.format(key))
        
        for i in directory:
           
           if 'Flute_4' in i:
              data_dict = f(i, data_dict)

      if os.path.isfile('~/{}.csv'.format(self.flute) ):
           df = pd.DataFrame(data=data_dict, columns = None, index=[0])
      else:
           df = pd.DataFrame(data=data_dict, columns = data_dict.keys(), index=[0])
    
      df_to_csv(df, self.flute)






#####################################################################################################################


def parse_args():
        
        parser = argparse.ArgumentParser()
        parser.add_argument('path')
       #parser.add_argument('test')
        return parser.parse_args()

 


def main():
      
      
   args = parse_args()
    
     
   for subdir, dirs, files in os.walk(args.path):  
    for directory in dirs:

  
      label = '_'.join(directory.split('_')[1:])
      directory = args.path+ directory+'/'
      flu1 = Flute1(directory, label) 
      #name= directory.path.split('/')[-2]
      flu1.run()
      flu2 = Flute2(directory, label)
      flu2.run()
      flu3 = Flute3(directory, label)
      flu3.run()
      flu4 = Flute4(directory, label)
      flu4.run()
          



if __name__=="__main__":
         main()

