import os

class PQC_RawData:
    '''
    '''

    def __init__(self,path,test,meta,series):
        self.data={}
        
        self.path=path
        self.test=test

        self.location='Hephy'
        
        self.sample_name=meta.get('sample_name')
        self.sample_position=meta.get('sample_position')
        self.sample_comment=meta.get('sample_comment')
        self.contact_name=meta.get('contact_name')
        self.measurement_name=meta.get('measurement_name')
        self.start_timestamp=meta.get('start_timestamp')
        self.operator=meta.get('operator')
        self.waiting_time=meta.get('waiting_time')
        
        self.FILE_NAME=os.path.basename(self.path)
        self.WAITING_TIME_S=self.waiting_time.split(' ')[0]
        self.RUN_BEGIN_TIMESTAMP=self.start_timestamp.replace('T',' ')
        man,batch,wafer,sensortype,hm,location=self.sample_name.split('_')
        self.NAME_LABEL=self.sample_name[7:]#removing 'HPK_VPX'
        self.KIND_OF_PART='{} Halfmoon {}'.format(sensortype.replace('-',''),location[0])
        self.KIND_OF_HM_SET_ID={'L':'Left','R':'Right'}[location[1]]
        
        self.KIND_OF_HM_FLUTE_ID,self.KIND_OF_HM_STRUCT_ID,self.KIND_OF_HM_CONFIG_ID=self.get_structure(self.contact_name,self.measurement_name)
        
    def add_data(self,data_dict):
        self.data={**self.data,**data_dict}

    def get_structure(self,contact_name,measurement_name):
        lookup={
            'FET':['PQC1','FET_PSS','Not Used'],
            'MOS capacitor (HV Source)':['PQC1','MOS_QUARTER','Not Used'],
            'Capacitor test structure Left 10kHz 250mV (HV Source)':['PQC1','CAP_W','Not Used'],
            'Capacitor test structure Right 10kHz 250mV (HV Source)':['PQC1','CAP_E','Not Used'],
            'Polysilicon Van-der-Pauw cross':['PQC1','VDP_POLY','STANDARD'],
            'Reverse Polysilicon Van-der-Pauw cross':['PQC1','VDP_POLY','ROTATED'],
            'N+ Van-der-Pauw cross':['PQC1','VDP_STRIP','STANDARD'],
            'Reverse N+ Van-der-Pauw cross':['PQC1','VDP_STRIP','ROTATED'],
            'P-stop Van-der-Pauw cross':['PQC1','VDP_STOP','STANDARD'],
            'Reverse P-stop Van-der-Pauw cross':['PQC1','VDP_STOP','ROTATED'],
            'GCD':['PQC2','GCD','Not Used'],
            'N+ linewidth structure':['PQC2','LINEWIDTH_STRIP','Not Used'],
            'P-stop linewidth structure (2-wire)':['PQC2','LINEWIDTH_STOP','Not Used'],
            'P-stop linewidth structure (4-wire)':['PQC2','LINEWIDTH_STOP','Not Used'],
            'Polysilicon meander':['PQC2','R_POLY','Not Used'],
            'Dielectric Breakdown 1':['PQC2','DIEL_SW','Not Used'],
            'Diode IV':['PQC3','DIODE_HALF','Not Used'],
            'Diode CV':['PQC3','DIODE_HALF','Not Used'],
            'Metal clover leaf Van-der-Pauw':['PQC3','CLOVER_METAL','STANDARD'],
            'Reverse Metal clover leaf Van-der-Pauw':['PQC3','CLOVER_METAL','ROTATED'],
            'P+ cross-bridge Van-der-Pauw':['PQC3','VDP_EDGE','STANDARD'],
            'Reverse P+ cross-bridge Van-der-Pauw':['PQC3','VDP_EDGE','ROTATED'],
            'P+ cross-bridge linewidth':['PQC3','VDP_EDGE','LINEWIDTH'],
            'Bulk cross':['PQC3','VDP_BULK','STANDARD'],
            'Reverse bulk cross':['PQC3','VDP_BULK','ROTATED'],
            'Metal meander':['PQC3','MEANDER_METAL','Not Used'],
            'GCD05':['PQC4','GCD05','Not Used'],
            'N+ CBKR':['PQC4','CBKR_STRIP','STANDARD'],
            'Polysilicon CBKR':['PQC4','CBKR_POLY','STANDARD'],
            'Polysilicon contact chain':['PQC4','CC_POLY','Not Used'],
            'P+ contact chain':['PQC4','CC_EDGE','Not Used'],
            'N+ contact chain':['PQC4','CC_STRIP','Not Used']}
        FLUTE_ID,STRUCT_ID,CONFIG_ID=lookup[measurement_name]

        if FLUTE_ID[-1] != contact_name[-1]:
            raise ValueError('Flute nr. mismatch',self.sample_name,self.test)

        return FLUTE_ID,STRUCT_ID,CONFIG_ID
    

        
    #def get_param() #needed for xml file

        
